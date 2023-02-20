#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH -t 2-00:00 # Runtime in D-HH:MM
#SBATCH -o /expanse/lustre/projects/ddp195/eiovino/melt_HG002/melt.out # File to which STDOUT will be written 
#SBATCH -e /expanse/lustre/projects/ddp195/eiovino/melt_HG002/melt.err # File to which STDERR will be written 
#SBATCH --job-name=HG002_melt_preprocessing
#SBATCH --mail-type=END # Type of email notification- BEGIN,END,FAIL,ALL 
#SBATCH --mail-user=emanuelaiovino.87@gmail.com # Email to which notifications wil

source /home/eiovino/miniconda3/bin/activate
conda activate melt
date
#java -jar /expanse/lustre/projects/ddp195/eiovino/MELTv2.2.2/MELT.jar Preprocess -bamfile /expanse/lustre/projects/ddp195/eiovino/cram_REACH000236/REACH000236.markdup.recal.sort.bam -h /expanse/lustre/projects/ddp195/eiovino/fasta/Homo_sapiens_assembly38.fasta 
#Preprocess generates three files, all as suffixes to the current bam file: sorted.bam.disc – discordant pairs from the current BAM file,sorted.bam.disc.fq,index of these discordant pairs, fastq version of all discordant pairs for MEI alignment
bash /expanse/lustre/projects/ddp195/eiovino/melt_HG002/NA12878/individualAnalysis.sh
date
