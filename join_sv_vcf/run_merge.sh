#!/bin/bash
#
#	bash bash run_merge.sh list_78_sample.list SV_cohort.vcf
#
#Input is SURVIVOR vcfs lists and the output vcf 

merged_vcf=$2
list_vcf=$1

SURVIVOR merge $list_vcf 1000 1 1 1 1 50 $merged_vcf

# # bcftools reheader takes in a pre-generated "new_header.txt" file that fixes some of the lines in the output SURVIVOR merged VCF, as well as re-ordering the contigs in the header (using the same order as the GRCh38 FASTA file). We also rename the samples such that they share the same names as the ones in the SNV VCF, and re-order them as well.
# # bcftools view --exclude "(END-POS)>1000000" gives the same number of variants as bcftools view --exclude "abs(SVLEN)>1000000"
bcftools reheader --header files/new_header.txt --samples files/renamed_samples.txt $merged_vcf | bcftools view --samples-file files/snv_samples.txt | bcftools view --exclude "abs(SVLEN)>1000000" | bcftools sort --output-type z --output-file $merged_vcf.gz

tabix $merged_vcf.gz
# merged_bed=ita_CNV_ls1MB.bed
bcftools query -f "%CHROM\t%POS0\t%END\t%SVTYPE\n" $merged_vcf.gz | grep --perl "DEL|DUP" > $merged_vcf.bed
python /work/emanuela.iovino/ClassifyCNV/ClassifyCNV.py --infile $merged_vcf.bed --GenomeBuild hg38 --outdir output_amcg
python merge_acmg_with_survivor.py $merged_vcf.gz 
python make_scores_for_survivor.py
