#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH -t 2-00:00 # Runtime in D-HH:MM
#SBATCH -o /expanse/lustre/projects/ddp195/eiovino/svaba/HG002_bwa/logs/svaba.out # File to which STDOUT will be written 
#SBATCH -e /expanse/lustre/projects/ddp195/eiovino/svaba/HG002_bwa/logs/svaba2.err # File to which STDERR will be written 
#SBATCH --job-name=sbatch_svaba
#SBATCH --mail-type=END # Type of email notification- BEGIN,END,FAIL,ALL 
#SBATCH --mail-user=emanuelaiovino.87@gmail.com # Email to which notifications wil

source /home/eiovino/miniconda3/bin/activate
conda activate sv_caller

date
svaba run -t /expanse/lustre/projects/ddp195/eiovino/cram_HG002/HG002.100_reads.markdup.recal.bam  -I -a //expanse/lustre/projects/ddp195/eiovino/svaba/HG002_bwa/germline_HG002 -p 64  -G /expanse/lustre/projects/ddp195/j3guevar/resources/GRCh38_reference_genome/testing/GRCh38_full_analysis_set_plus_decoy_hla.fa

date
