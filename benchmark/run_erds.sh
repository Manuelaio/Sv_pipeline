#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH -t 2-00:00 # Runtime in D-HH:MM
#SBATCH -o /home/eiovino/sv_callers/ERDS/erds.out # File to which STDOUT will be written 
#SBATCH -e /home/eiovino/sv_callers/ERDS/erds.err # File to which STDERR will be written 
#SBATCH --job-name=sbatch_erds
#SBATCH --mail-type=END # Type of email notification- BEGIN,END,FAIL,ALL 
#SBATCH --mail-user=emanuelaiovino.87@gmail.com # Email to which


perl /home/eiovino/sv_callers/ERDS/erds1.1/erds_pipeline.pl -b /expanse/projects/sebat1/genomicsdataanalysis/pilot_data/ASH_trio/processed_data/hg002_ill/sorted_indexed/HG002.GRCh38.300x.bam  -v /expanse/lustre/projects/ddp195/j3guevar/sv_testing/resources/snv_indel_vcfs_for_cerds/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz -o /expanse/lustre/projects/ddp195/eiovino/erds/ -r /expanse/lustre/projects/ddp195/j3guevar/sv_testing/resources/GRCh38_full_plus_hs38d1_analysis_set_minus_alts.fa 
