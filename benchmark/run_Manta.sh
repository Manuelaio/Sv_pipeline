#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH -t 2-00:00 # Runtime in D-HH:MM
#SBATCH -o Manta.out # File to which STDOUT will be written 
#SBATCH -e Manta.err # File to which STDERR will be written 
#SBATCH --job-name=sbatch_Manta_NA12878
#SBATCH --mail-type=END # Type of email notification- BEGIN,END,FAIL,ALL 
#SBATCH --mail-user=emanuelaiovino.87@gmail.com # Email to which notifications wil

source /home/eiovino/miniconda3/bin/activate
conda activate sv_caller

date
configManta.py --bam NA12878.alt_bwamem_GRCh38DH.20150706.CEU.illumina_platinum_ped.cram --referenceFasta GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa --runDir  NA12878/
python NA12878/runWorkflow.py
date
