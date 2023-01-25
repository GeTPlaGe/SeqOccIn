# all the software packages used to run this script should be installed before hand and available in the $PATH

# for bioinfo genotoul users they can be uploaded by module in order to do this you have to decomment the following lines
#module load bioinfo/jellyfish-2.2.10
#module load bioinfo/wtdbg2-2.3
#module load bioinfo/supernova-2.1.1
#module load bioinfo/racon-v1.4.10
#module load bioinfo/pilon-v1.2
#module load bioinfo/hifiasm-v0.15.5
#module load bioinfo/juicer-1.5.6
#module load bioinfo/minimap2-2.5
#module load bioinfo/busco-5.1.2
#module load bioinfo/Inspector-v1.2
#module load bioinfo/seqkit-v0.16.0
#module load bioinfo/RepeatModeler-open-1.0.11
#module load bioinfo/RepeatMasker-4-0-7
#module load bioinfo/yak-0.1
#module loadbioinfo/dsk-v2.2.0


# generating histogram file for genomescope2 
zcat -f ERR1037805[45678].fastq.gz | jellyfish count /dev/fd/0 -C -m 21 -o reads.jf
jellyfish histo -h 10000000 reads.jf > reads.histo


# ONT wtdbg2 v2.3 assembly
zcat ERR10386215.fastq.gz ERR10386216.fastq.gz ERR10386217.fastq.gz ERR10386218.fastq.gz ERR10386219.fastq.gz > ONT.fastq
wtdbg2 -g 2.9g -x ont -i  -fo assembly_ONT_raw 
wtpoa-cns -i assembly_ONT_raw.ctg.lay.gz -fo assembly_ONT_raw.ctg.fa

# CLR wtdbg2 v2.3 assembly
wtdbg2 -g 2.9g -x sq -i ERR10378053.fastq.gz -fo assembly_CLR_raw 
wtpoa-cns -i assembly_CLR_raw.ctg.lay.gz -fo assembly_CLR_raw.ctg.fa

# Supernova} v2.1.1 assembly
mv ERR10310247_1.fastq.gz ERR10310247_2.fastq.gz ERR10310248_1.fastq.gz ERR10310248_2.fastq.gz ERR10310249_1.fastq.gz ERR10310249_2.fastq.gz ERR10310250_1.fastq.gz ERR10310250_2.fastq.gz ERR10310251_1.fastq.gz ERR10310251_2.fastq.gz ERR10310252_1.fastq.gz ERR10310252_2.fastq.gz  10X_Chromium_Directory/
supernova run --id=assembly --maxreads=1082666666 --fastqs=10X_Chromium_Directory/
supernova mkoutput --style=pseudohap --asmdir=assemblyDir --outprefix=assembly


# Racon v1.4.10 polishing
minimap2 -a -x map-ont assembly_ONT_raw.ctg.fa ONT.fastq > alignment.sam
racon ONT.fastq alignement.sam assembly_ONT_raw.cgt.fa > assembly_ONT_racon.fa


# Pilon v1.22 polishing
longranger mkref assembly_racon.fa longranger align \
  --reference=assembly_ONT_racon.fa --id=ALIGN \
  --fastqs=10X_Chromium_Directory 
java -jar pilon-1.22.jar --genome assembly_racon.fa \
  --bam longranger_alignement.bam --fix bases,gaps --output assembly_pilon.fa
# racon and pilon polishing can also be performed in the same way on the wtdbg2 assembly

# Hifiasm v0.15.5 assembly
zcat ERR1037805[45678].fastq.gz > CCS_reads.fastq
hifiasm -o assembly_HiFi.asm CCS_reads.fastq
awk '/^S/{print ">"$2"\n"$3}' assembly_HiFi.asm.bp.p_ctg.gfa | fold > assembly_HiFi.asm.bp.p_ctg.gfa.fa


# Hifiasm v0.15.5 diploid Hi-C assembly
zcat ERR10310241_1.fastq.gz ERR10310242_1.fastq.gz ERR10310243_1.fastq.gz ERR10310244_1.fastq.gz > HiC_1.fastq
zcat ERR10310241_2.fastq.gz ERR10310242_2.fastq.gz ERR10310243_2.fastq.gz ERR10310244_2.fastq.gz > HiC_2.fastq
hifiasm -o assembly_HiFi_HiC.asm --h1 HiC_1.fastq \
    --h2 HiC_2.fastq CCS_reads.fastq
awk '/^S/{print ">"$2"\n"$3}' assembly_HiFi_HiC.asm.bp.hap1.p_ctg.gfa | fold > assembly_HiFi_HiC.asm.bp.hap1.p_ctg.gfa.fa
awk '/^S/{print ">"$2"\n"$3}' assembly_HiFi_HiC.asm.bp.hap2.p_ctg.gfa | fold > assembly_HiFi_HiC.asm.bp.hap2.p_ctg.gfa.fa

# Hifiasm v0.15.5 diploid parental assembly
yak count -b37 -o mother.yak illumina_mother_reads1.fastq.gz \
  illumina_mother_reads2.fastq.gz
yak count -b37 -o father.yak illumina_father_reads1.fastq.gz \
  illumina_father_reads2.fastq.gz
hifiasm -o assembly_HiFi_Parental.asm -1 mother.yak \
  -2 Father.yak CCS_reads.fastq
awk '/^S/{print ">"$2"\n"$3}' assembly_HiFi_Parental.asm.bp.hap1.p_ctg.gfa | fold > assembly_HiFi_Parental.asm.bp.hap1.p_ctg.gfa.fa
awk '/^S/{print ">"$2"\n"$3}' assembly_HiFi_Parental.asm.bp.hap2.p_ctg.gfa | fold > assembly_HiFi_Parental.asm.bp.hap2.p_ctg.gfa.fa


# Hi-C scaffolding Juicer + 3D-Dna v180114 heatmap production and scaffold construction 
mkdir references
cp assembly_HiFi.asm.bp.p_ctg.gfa.fa references/assembly.fa
bwa index references/assembly.fa
python generate_site_positions.py HindIII assembly.fa \
  ../references/assembly.fa > restriction_sites/assembly_HindIII.txt
  
awk '{print $1,$NF}' restriction_sites/assembly_HindIII.txt > \
  references/chrom_size.tsv
  
scripts/juicer.sh -g Bos_taurus -s HindIII -d . \
  -z references/assembly.fa -p references/chrom_size.tsv \
  -y restriction_sites/assembly_HindIII.txt -S early -D .

run-asm-pipeline.sh -r 0 assembly.fa juicer_merged_nodups.txt
run-asm-pipeline-post-review.sh --sort-output -r genome.assembly \
  assembly.fa juicer_merged_nodups.txt


# 10X: Scripts and github repository to produce 3DDna Hi-C-like heatmap, with 10X data.
# http://genoweb.toulouse.inra.fr/~sigenae/10X_and_Juicebox.html#HiC_asm
# https://forgemia.inra.fr/patrice.dehais/tag_compare

# PhasingValidation, For assemblies, with dsk v2.2:
dsk -file assembly_HiFi_Parental.asm.bp.hap1.p_ctg.gfa.fa -kmer-size 21 -out dsk_hap1 
dsk -file assembly_HiFi_Parental.asm.bp.hap2.p_ctg.gfa.fa -kmer-size 21 -out dsk_hap2
dsk2ascii -file dsk_hap1.h5 -out hap1.tsv
dsk2ascii -file dsk_hap2.h5 -out hap2.tsv
sort hap1.tsv > hap1_sorted.tsv
sort hap2.tsv > hap2_sorted.tsv
join -j 1 -v 2 hap1_sorted.tsv hap2_sorted.tsv > joinedHapv2.tsv
join -j 1 -v 1 hap1_sorted.tsv hap2_sorted.tsv > joinedHapv1.tsv
join hap1_sorted.tsv hap2_sorted.tsv > joinedBothHap.tsv


# For reads, with dsk v2.2:
dsk -file reads1.fa -kmer-size 21 -out dsk_reads1 
dsk -file reads2.fa -kmer-size 21 -out dsk_reads2
dsk2ascii -file dsk_reads1.h5 -out reads1.tsv
dsk2ascii -file dsk_reads2.h5 -out reads2.tsv
sort reads1.tsv > reads1_sorted.tsv
sort reads2.tsv > reads2_sorted.tsv
join -j 1 -v 2 reads1_sorted.tsv reads2_sorted.tsv > joinedReadsv2.tsv
join -j 1 -v 1 reads1_sorted.tsv reads2_sorted.tsv > joinedReadsv1.tsv
join reads1_sorted.tsv reads2_sorted.tsv > joinedBothReads.tsv

# Evaluate phasing quality:
join joinedReadsv1.tsv joinedHapv1.tsv > Shared_Reads1_Hap1.tsv
join joinedReadsv1.tsv joinedHapv2.tsv > Shared_Reads1_Hap2.tsv
join joinedReadsv2.tsv joinedHapv1.tsv > Shared_Reads2_Hap1.tsv
join joinedReadsv2.tsv joinedHapv2.tsv > Shared_Reads2_Hap2.tsv
# To reproduce analysis on each contigs, use the contigs sequence as Assembly.hap1.fa or Assembly.hap2.fa


# minimap v2.5 alignment against reference for DGenies visualization
minimap2 -cx asm5 reference.fa assembly.fa > map.paf


# BUSCO v5.1
busco -i assembly.fa -o Busco_assembly -l mammalia_odb10 -m geno

# Inspector v1.2
inspector.py -c assembly.fasta -r reads.fastq.gz \
   --ref ARS-UCD1.2_genomic.fa.gz -o inspectorOut 

#seqkit v0.16.0
seqkit fx2tab -q -l -n reads.fastq.gz > stats.tsv


# RepeatModeler / repeatMasker  v1.0.11 / v4-0-7, for HiFi and ONT assemblies
BuildDatabase -name Off2_HiFi_noMatch assembly_HiFi_Parental.asm.bp.hap1.p_ctg.gfa.fa
RepeatModeler -database Off2_HiFi_noMatch -pa 10
RepeatMasker -pa 10 -lib consensi.fa.masked -dir RMasker_Hifidb assembly_HiFi.fa


