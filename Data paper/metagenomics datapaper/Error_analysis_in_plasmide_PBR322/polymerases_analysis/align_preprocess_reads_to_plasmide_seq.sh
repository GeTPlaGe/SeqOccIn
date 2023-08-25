#!/bin/bash
#SBATCH -J minimap_plasmide
#SBATCH -c 8
#SBATCH --mem=20GB
#SBATCH -t 03-00
#SBATCH --output=slurm-%x.%j.out
#SBATCH -p workq


module load bioinfo/samtools-1.16.1
module load bioinfo/minimap2-2.24

# reference sequence downloaded manually here: 
# https://international.neb.com/-/media/nebus/page-images/tools-and-resources/interactive-tools/dna-sequences-and-maps/text-documents/pbr322fsa.txt?rev=b1ef8762bbe9431886127f019c09df5b&hash=355F5D1EE26BAC8271AA68193BF15F36
# and cut on HindIII restriction site to match amplicon expected sequence 
ref='plasmide_Pbr322_neb_HindIII_cut.fasta'

minimap2 -ax map-hifi $ref preprocess.fasta --cs  -t 8 | samtools sort -o preprocess_reads_to_plasmide_seq.bam --threads 8

samtools index preprocess_reads_to_plasmide_seq.bam -@ 8

samtools flagstat preprocess_reads_to_plasmide_seq.bam  --threads 8 > preprocess_reads_to_plasmide_seq.flagstat 
 
 
