#!/bin/bash
#SBATCH -J launcher_collect_error_events_in_alignment
#SBATCH -c 1
#SBATCH --mem=2GB
#SBATCH -t 03-00
#SBATCH --output=slurm-%x.%j.out
#SBATCH -p workq


source /usr/local/bioinfo/src/Anaconda/Anaconda3-5.2.0/etc/profile.d/conda.sh
conda activate compute_error_rate_env

python collect_error_events_in_alignment.py --bam_file shotgun_reads_to_plasmide.bam --debug
 
 
