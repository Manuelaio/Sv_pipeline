delly call -g $reference_fasta -x $delly_exclude_regions_bed  -o tmp.bcf $alignment_file -q 20 -s 15 -z 5
delly call -g $reference_fasta -v tmp.bcf -x $delly_exclude_regions_bed -o $sample_id-delly.bcf  $alignment_file -q 20 -s 15 -z 5

#### use a homemade script 
python /expanse/lustre/projects/ddp195/eiovino/sv_pipeline_nf/data/recode_delly.py NA12878.bcf NA12878.recode.vcf
#### remove small calls 
bcftools view --include "SVLEN>=50 || SVLEN<=-50" NA12878.recode.vcf | bcftools view --include "FILTER='PASS'" > NA12878_vcf-filtered.vcf
