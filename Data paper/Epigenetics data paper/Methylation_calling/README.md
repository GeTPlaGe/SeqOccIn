# SeqOccIn - CpG methylation calling from short and long read sequencing technologies

Methylation calling for pig and quail was performed according to ONT, PacBio, EM-seq and WGBS methods used in the paper.

## EM-Seq and WGBS methylation calling with nfcore/methylseq

Methylation extractions from EM-seq and WGBS reads were performed using the [nfcore/methylseq](https://github.com/nf-core/methylseq) pipeline v1.6.1 and [Nextflow](https://github.com/nextflow-io/nextflow) v21.04.1, following the bismark workflow as follows:

```
nextflow run nf-core/methylseq \
--fasta /path/to/reference/genome \
--input '/path/to/files/*R{1,2}.fastq.gz' \
-profile genotoul \ # cluster configuation 
--zymo \ # if WGBS data
--maxins 800 \ # if WGBS data
--em_seq \ # if EM-Seq data
--cytosine_report
```

The output from bismark are stranded CpG report files in tab-delimited text format and contain information about the methylation status of all cytosines :

|chrm | position| strand |count methylated|count unmethylated|C-context|trinucleotide context|
| ----| ----| ----| ----| ----| ----| ----|
KQ966912.1|236|+| 0| 0| CG|CGC
KQ966912.1|237|-| 0| 0| CG|CGC
KQ966912.1|903|+| 0| 0| CG|CGG
KQ966912.1|904|-| 0| 0| CG|CGT
KQ966912.1|1118|+| 0| 0| CG|CGG

## ONT methylation calling with Megalodon

CpG methylation calling from ONT reads was performed using [Megalodon](https://github.com/nanoporetech/megalodon) v2.5.0 along with guppy v5.0.17 using res_dna_r941_prom_modbases_5mC_CpG_v001.cfg as calling model with default parameters as follows:

```
megalodon /path/to/fast5/ \
    --outputs basecalls mods \
    --guppy-params "-d /path/to/rerio/basecall_models/" \
    --guppy-config res_dna_r941_min_modbases_5mC_CpG_v001.cfg \
    --reference /path/to/reference/genome \
    --mod-motif m CG 0 --devices 0 1 2 3 --processes 30 \ # 4 GPU
    --verbose-read-progress 3 \
    --guppy-server-path /path/to/guppy_basecall_server \
    --overwrite
```

The output of the command is a `megalodon_results` folder holding the basecalls in fastq format and the modified bases in bedmethyl format :

|chrm | start| stop|score|coverage|strand|start2|stop2|RGB|count|fraction|
| ----| ----| ----| ----| ----| ----| ----| ----| ----| ----| ----|
|CM003781.1|54|55|.|5|+|54|55|0,0,0|5|80.0
CM003781.1|55|56|.|5|-|55|56|0,0,0|5|100.0
CM003781.1|106|107|.|4|+|106|107|0,0,0|4|100.0
CM003781.1|107|108|.|4|-|107|108|0,0,0|4|100.0
CM003781.1|166|167|.|5|+|166|167|0,0,0|5|60.0

## PacBio methylation calling with primrose

[Primrose](https://github.com/mattoslmp/primrose) software package in [SMRTLink v11.0](https://www.pacb.com/wp-content/uploads/SMRT_Tools_Reference_Guide_v11.0.pdf) was used to call CpG methylation from HiFi reads extracted from the subreads bam file using the ccs command line tools with `--hifi-kinetics` option as follows:

```
ccs movie.subreads.bam movie.ccs.kinetics.bam --hifi-kinetics
primrose --log-level INFO movie.ccs.kinetics.bam movie.ccs.5mc.bam
```

The output is adhering to the [SAM tag specification from 9. Dec 2021](https://samtools.github.io/hts-specs/SAMtags.pdf),
using `Mm` and `Ml` tags. It's also described in the [PacBio BAM file formats](https://pacbiofileformats.readthedocs.io/en/latest/BAM.html#use-of-read-tags-for-per-read-base-base-modifications) as :

| Tag  | Type  |           Description            |
| ---- | ----- | -------------------------------- |
| `MM` | `Z`   | Base modifications / methylation |
| `ML` | `B,C` | Base modification probabilities  |

This bam file was next mapped to the reference using [pbmm2](https://github.com/PacificBiosciences/pbmm2) and converted in bedmethyl format with [modbam2bed](https://github.com/epi2me-labs/modbam2bed) script from ONT :

```
pbmm2 align --sort reference.mmi movie.ccs.5mc.bam movie.ccs.5mc.align.bam
modbam2bed -e -m 5mC -t8 --cpg reference.fa movie.ccs.5mc.align.bam > movie.ccs.bedmethyl
```

The output obtained in bedmethyl format looks like this :

|chrm | start| stop|score|coverage|strand|start2|stop2|RGB|count|fraction|
| ----| ----| ----| ----| ----| ----| ----| ----| ----| ----| ----|
CM003781.1|25|26|5mC|500|+|25|26|0,0,0|2|100.00|0|1|1
CM003781.1|54|55|5mC|1000|+|54|55|0,0,0|2|100.00|0|2|0
CM003781.1|55|56|5mC|1000|-|55|56|0,0,0|1|100.00|0|1|0
CM003781.1|106|107|5mC|1000|+|106|107|0,0,0|2|100.00|0|2|0
CM003781.1|107|108|5mC|1000|-|107|108|0,0,0|1|100.00|0|1|0
