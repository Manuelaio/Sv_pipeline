#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH -t 2-00:00 # Runtime in D-HH:MM
#SBATCH -o /expanse/lustre/projects/ddp195/eiovino/lumpy/HG002_bwa/logs/lumpy.out # File to which STDOUT will be written 
#SBATCH -e /expanse/lustre/projects/ddp195/eiovino/lumpy/HG002_bwa/logs/lumpy.err # File to which STDERR will be written 
#SBATCH --job-name=sbatch_lumpy_NA12878
#SBATCH --mail-type=END # Type of email notification- BEGIN,END,FAIL,ALL 
#SBATCH --mail-user=emanuelaiovino.87@gmail.com # Email to which notifications wil

source /home/eiovino/miniconda3/bin/activate
conda activate sv_caller

date


smoove call --name NA12878 --exclude $smoove_exclude_regions_bed --fasta $reference_fasta  -p $task.cpus --genotype $alignment_file


date
