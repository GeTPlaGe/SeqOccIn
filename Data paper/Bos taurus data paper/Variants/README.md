# SeqOccIn - A Bos taurus data paper 

## Structural variants

The variants described in the paper were detected using the following approaches:

- The long reads based variation detection was performed using a dedicated pipeline available at https://github.com/SeqOccin-SV/SeqOccinVariants
- The short reads based variant detection was performed using the manta sofwtare (https://github.com/Illumina/manta) v1.6.0, with default parameters, more precisely the following commands were used

```
bwa mem -M -T 30 reference.fa reads_R1.fastq.gz reads_R2.fastq.gz | \
                  samtools view -bS -  | 
                  samtools sort -o alignments.bam
                  
configManta.py  --bam alignments.bam \
                --referenceFasta reference.fa \
                --runDir  RUNDIR \
                --callRegions regions.bed.gz
python2 RUNDIR/runWorkflow.py
```
- The assembly based variant detection was permformed using svim-asm (https://github.com/eldariont/svim-asm) according to recommandations, with minimap2 (v2-2.11) and svim-asm (v1.0.2)

```
minimap2 -a -x asm5 --cs -r2k reference.fa haplotype1.fasta > haplotype1.sam
minimap2 -a -x asm5 --cs -r2k reference.fa haplotype1.fasta > haplotype2.sam
samtools sort -o happlotype1.sorted.bam  happlotype1.sam
samtools sort -o happlotype2.sorted.bam  happlotype2.sam

svim-asm haploid OUTDIR happlotype1.sorted.bam happlotype2.sorted.bam reference.fa
```

All the vcf files are aviable here:

- Variants detected using PacBio HiFi reads (40X) :   
  https://entrepot.recherche.data.gouv.fr/file.xhtml?persistentId=doi:10.57745/7UB49W&version=1.0
- Variants detected using ONTeads (42X) :   
  https://entrepot.recherche.data.gouv.fr/file.xhtml?persistentId=doi:10.57745/QVUJCA&version=1.0


