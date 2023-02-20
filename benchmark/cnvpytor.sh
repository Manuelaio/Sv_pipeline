#1 Download ad conda

bam=$1
name=$2
fasta=$3

#samtools idxstats $bam | cut -f 1 > chr.txt
#awk '{if($1 !~ "_" && $1 !~ "HLA-")print $0}' test.bed > no_decoy.bed
value=`cat no_decoychr.txt`
cnvpytor -root $name.pytor -chrom $value -rd $bam -T $fasta
cnvpytor -root $name.pytor -his 100

cnvpytor -root $name.pytor -partition 100
cnvpytor -root $name.pytor  -call 100 > $name.calls.tsv
awk '{ print $2 } END { print "exit" }' $name.calls.tsv | cnvpytor -root $name.pytor  -genotype 100
#perl cnvnator2VCF.pl -prefix HG002 -reference GRCh38 $name.calls.tsv /archive/ngsbo/db/trioCEU_1KGP_resources/GRCh38_full_analysis_set_plus_decoy_hla.fa > CNVcall.1000.vcf 
