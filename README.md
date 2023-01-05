# Strucural Varinats  pipeline
Nextflow pipeline for short read WGS data

```
nextflow run sv_main.nf --reference_fasta GRCh38_full_analysis_set_plus_decoy_hla.fa --sample_alignments_tsv sample_sv2.csv --outdir /. --delly_exclude_regions_bed Hg38/exlude.regions.delly.human.hg38.excl.tsv --smoove_exclude_regions_bed Hg38/exclude.smoove.cnvnator_100bp.GRCh38.20170403.bed --expansion_hunter_variant_catalog_json Hg38/ExpansioHunter_variant_catalog.json --recode_delly_python_script Hg38/recode_delly.py --mosdepth_segmental_duplications_bed Hg38/Segmental_dups_hg38_frt_srt.bed --bind_path . -profile slurm -resume

```
