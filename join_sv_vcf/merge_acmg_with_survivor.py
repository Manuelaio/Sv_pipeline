import pysam
import sys

scoresheet_dict = {}
with open("output_amcg/Scoresheet.txt", "r") as f:
    header = f.readline()
    for line in f:
        linesplit = line.rstrip().split("\t")
        VariantID, Chromosome, Start, End, Type, Classification = linesplit[0:6]
        scoresheet_dict[VariantID] = Classification

vcf_path = sys.argv[1]
vcf_iterator = pysam.VariantFile(vcf_path, mode = "r")
vcf_out = pysam.VariantFile("ita_cohort_survivor_intersect_merged_ACMG.vcf.gz", mode = "w", header = vcf_iterator.header)
vcf_out.header.add_meta(key = "INFO", items = [ ("ID", "ACMG_Classification"), ("Number", 1), ("Type", "Integer"), ("Description", "ACMG Classification (Coding is this: Benign: 0, Likely benign: 1, Uncertain significance: 2, Likely pathogenic: 3, Pathogenic: 4)") ] )

classification_dict = {"Benign": 0, "Likely benign": 1, "Uncertain significance": 2, "Likely pathogenic": 3, "Pathogenic": 4}

for record in vcf_iterator:
    VariantID = "{}_{}_{}_{}".format(record.chrom, record.start, record.stop, record.info["SVTYPE"])

    if VariantID not in scoresheet_dict:
        vcf_out.write(record)
        continue

    new_record = vcf_out.header.new_record(contig = record.chrom, start = record.start, stop = record.stop, alleles = record.alleles, id = record.id, qual = record.qual, filter = record.filter, info = record.info)
    new_record.info["ACMG_Classification"] = classification_dict[scoresheet_dict[VariantID]]
    for i in range(len(record.samples)):
        for format_ in record.samples[i]: new_record.samples[i][format_] = record.samples[i][format_]
    vcf_out.write(new_record)
