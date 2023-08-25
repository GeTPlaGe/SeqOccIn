#!/bin/bash
#SBATCH -J frogs_preprocess
#SBATCH -c 6
#SBATCH --mem=20GB
#SBATCH -t 03-00 #Acceptable time formats include "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds".
#SBATCH --output=slurm-%x.%j.out
#SBATCH -p workq


source /usr/local/bioinfo/src/Anaconda/Anaconda3-5.2.0/etc/profile.d/conda.sh
conda activate frogs@4.1.0

five_prime_primer=GCAGTCGAACATGTAGCTGACTCAGGTCACAGCTTTAATGCGGTAGTTTATC
#three_prim_primer=TGGATCACTTGTGCAAGCATCACATCGTAGATCGATGATAAGCTGTCAAAC  need to be reverse complement to work with frogs
three_prim_primer=GTTTGACAGCTTATCATCGATCTACGATGTGATGCTTGCACAAGTGATCCA


preprocess.py longreads --min-amplicon-size 4000  --max-amplicon-size 5000 \
                          --five-prim-primer $five_prime_primer --three-prim-primer $three_prim_primer -p 6 \
                          --input-R1 fastq/*


