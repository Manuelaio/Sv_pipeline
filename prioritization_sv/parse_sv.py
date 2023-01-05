#TO DO: - Add file for taking vep anootation column DONE
#       - add DP
#       - add output with header DONE
# running : python test_rab/parse_sv.py -i testing_vep/chr7_SV_test_ann_VEPann_11-18-22.vcf.gz -ann test_rab/ann_sv.txt -ped /expanse/lustre/projects/ddp195/eiovino/ita_cohort/SNVs/model.txt  -o test1.tsv


import argparse
import numpy as np
import pysam
import sys

argparser = argparse.ArgumentParser(description = 'Extracts SV annotations from VEP-annotated VCF')
argparser.add_argument('-i', '--in-vcf', metavar = 'file', dest = 'in_VCF', required = True, help = 'VEP-annotated VCF')
argparser.add_argument('-o', '--out', metavar = 'file', required = True, help = 'Output tab-delimited annotation file')
argparser.add_argument('-ann', '--vep_annotation', metavar = 'file', dest = 'in_ann', required = True)
argparser.add_argument('-ped', '--ped_file', metavar = 'file', dest = 'in_ped', required = True)
argparser.add_argument('-fid', type = str)
args = argparser.parse_args()

fout = open(args.out, "w")

#------- Open VCF file -------#
vcf_iterator = pysam.VariantFile(args.in_VCF)

#------- Open pedigree file -------#
samples = []
with open(args.in_ped, "r") as f:
    for line in f: 
        linesplit = line.rstrip().split("\t")
        if args.fid:
            if linesplit[0] != args.fid: continue
        samples.append(linesplit[1])
index_sample_dict = dict(enumerate(vcf_iterator.header.samples))
sample_index_dict = {v: k for k, v in index_sample_dict.items()}



#------- Define VEP annotation -------#
vep_annotations_subset = [] # The vep annotations I want to get
with open(args.in_ann, "r") as f:
    for line in f: vep_annotations_subset.append(line.rstrip())

#------- Index all the VEP annotations -------#
vep_meta_name = "ANN"
vep_meta = vcf_iterator.header.info.get("ANN", None) # Check header to make sure VEP was used
if not vep_meta: 
    vep_meta_name = "CSQ"
    vep_meta = vcf_iterator.header.info.get("CSQ", None)
if not vep_meta: raise Exception('The input file is NOT annotated with VEP, no ANN header found.')

vep_annotations_list = vep_meta.description.split(":", 1)[1].strip().split("|")
def index_annotations(vep_annotations_list):
    vep_index_dict = {}
    for i in range(len(vep_annotations_list)): vep_index_dict[vep_annotations_list[i]] = i
    return vep_index_dict
vep_index_dict = index_annotations(vep_annotations_list)
number_of_annotations = len(vep_annotations_list)

consequence_levels_dict = {
"transcript_ablation" : 0,
"splice_acceptor_variant" : 1,
"splice_donor_variant" : 2,
"stop_gained" : 3,
"frameshift_variant" : 4,
"stop_lost" : 5,
"start_lost" : 6,
"transcript_amplification" : 7,
"inframe_insertion" : 8,
"inframe_deletion" : 9,
"missense_variant" : 10,
"protein_altering_variant" : 11,
"splice_region_variant" :12,
"splice_donor_5th_base_variant" : 13,
"splice_donor_region_variant" : 14 ,
"splice_polypyrimidine_tract_variant" : 15,
"incomplete_terminal_codon_variant" : 16,
"start_retained_variant" : 17,
"stop_retained_variant" : 18,
"synonymous_variant" : 19,
"coding_sequence_variant" : 20,
"mature_miRNA_variant" : 21,
"5_prime_UTR_variant" : 22,
"3_prime_UTR_variant" : 23,
"non_coding_transcript_exon_variant" : 24,
"intron_variant" : 25,
"NMD_transcript_variant" : 26 ,
"non_coding_transcript_variant" : 27,
"upstream_gene_variant" : 28,
"downstream_gene_variant" : 29,
"TFBS_ablation" : 30,
"TFBS_amplification" : 31,
"TF_binding_site_variant" : 32,
"regulatory_region_ablation" : 33,
"regulatory_region_amplification" : 34,
"feature_elongation" : 35,
"regulatory_region_variant" : 36,
"feature_truncation" : 37,
"intergenic_variant" : 38,

}

def rank_transcripts(vep_annotations_per_transcript):
    consequence_levels = np.zeros(len(vep_annotations_per_transcript))
    symbols = []
    # Loop through each transcript
    for i, transcript_annotations in enumerate(vep_annotations_per_transcript):
        consequence = transcript_annotations.split("|")[vep_index_dict["Consequence"]].split("&")[0]
        if consequence not in consequence_levels_dict: consequence_levels[i] = 100
        else: consequence_levels[i] = consequence_levels_dict[consequence]
    return vep_annotations_per_transcript[np.argsort(consequence_levels)[0]], consequence_levels[0]

# Header
print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format("chrom", "pos", "ref", "alt", "AC","SVLEN","SVTYPE", "AMCG"), end = "\t", file = fout)
for annotation in vep_annotations_subset:
    print(annotation)
    print(annotation, end = "\t", file = fout)
print("SYMBOLS_ALL", end = "\t", file = fout)
for sample in samples:
    print("{}_GT".format(sample), end = "\t", file = fout)
#for sample in samples:
#    print("{}_DP".format(sample), end = "\t", file = fout)
for sample in samples:
    print("{}_SCORE".format(sample), end = "\t", file = fout)
print("FID", file = fout)

for record in vcf_iterator:
    vep_annotations_per_transcript = record.info[vep_meta_name]
    worst_transcript, transcript_level = rank_transcripts(vep_annotations_per_transcript)
    if transcript_level == 100: continue

    symbols = set()
    for transcript_annotations in vep_annotations_per_transcript: 
        symbols.add(transcript_annotations.split("|")[vep_index_dict["SYMBOL"]])

    sample_genotypes = []
    num_alternate = 0
    for sample in samples:
        sample_genotype = record.samples[sample_index_dict[sample]]["GT"]
        sample_genotypes.append(sample_genotype)
        if 1 in sample_genotype: num_alternate += 1 
    if num_alternate == 0: continue

    #sample_depths = []
    #for sample in samples:
    #    sample_depth = record.samples[sample_index_dict[sample]]["DP"]
    #    sample_depths.append(sample_depth)

    sample_DRs = []
    for sample in samples:
        sample_DR = record.samples[sample_index_dict[sample]]["SCORE"][0]
        sample_DRs.append(sample_DR)

    worst_transcript_list = worst_transcript.split("|")

    vep_annotations_subset_dict = {}
    for annotation in vep_annotations_subset:
        vep_annotations_subset_dict[annotation] = worst_transcript_list[vep_index_dict[annotation]]

    # Variant record information...
    # record.chrom, record.pos...
    alt = record.alts[0]
    filters = []
    for filter_ in record.filter: filters.append(filter_) 
    #AC = record.info["AC"][0]
    formatted_genotypes = []
    for genotype in sample_genotypes:
        #formatted_genotype = [-1, -1]
        formatted_genotype = ['.', '.']
        if genotype[0] != None: formatted_genotype[0] = genotype[0]
        if genotype[1] != None: formatted_genotype[1] = genotype[1]
        formatted_genotypes.append(formatted_genotype)
    AC = record.info['SUPP']
    svlen = record.info['SVLEN']
    svtype = record.info['SVTYPE']
    amcg = record.info['ACMG_Classification']
    #print(AC)
    # Printing...
    # Variant record data...
    print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(record.chrom, record.pos, record.ref, alt, AC, svlen, svtype, amcg), end = "\t", file = fout)
    #print("{}\t{}\t{}\t{}\t{}\t{}".format(record.chrom, record.pos, record.ref, alt, svlen, svtype), end = "\t", file = fout)
    # Variant annotation data
    for key, val in vep_annotations_subset_dict.items():
        print(val, end = "\t", file = fout)
    print(",".join(list(symbols)), end = "\t", file = fout)
    # Genotype data 
    for genotype in formatted_genotypes: 
        print("/".join(map(str, genotype)), end = "\t", file = fout)
    ## Depth data
    for dr in sample_DRs:
        print(dr, end = "\t", file = fout)

    print(args.fid, file = fout)
