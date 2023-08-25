# Get the fastq files 

To fill when the data will be uploaded


# Plasmid reference sequence 

The plasmid reference sequence is needed for detecting errors in amplicon sequences. You can locate the plasmid sequence on the NEB website via the following link:[NEB Plasmid Reference Sequence.](https://international.neb.com/-/media/nebus/page-images/tools-and-resources/interactive-tools/dna-sequences-and-maps/text-documents/pbr322fsa.txt?rev=b1ef8762bbe9431886127f019c09df5b&hash=355F5D1EE26BAC8271AA68193BF15F36)


# Error analysis of the polymerases

## Preprocessing of the sequencing reads of the plasmid Pbr322 amplified Accustart and Hotstart polymerases

The amplicon sequences from the plasmid Pbr322 will undergo preprocessing using the FROGS v4.1.0 preprocess tool. This step includes removing PacBio adapters and PCR primers, along with dereplicating identical sequences.


### Install FROGS 4.1.0

1. Get FROGS 4.1.0 repository
```bash
wget https://github.com/geraldinepascal/FROGS/archive/refs/tags/v4.1.0.tar.gz

tar -xf v4.1.0.tar.gz
rm v4.1.0.tar.gz
```
2. Install with conda

```bash
conda env create --name frogs@4.1.0 --file FROGS-4.1.0/frogs-conda-requirements.yaml
# to use FROGS, first you need to activate your environment
conda activate frogs@4.1.0
```

### Launch FROGS preprocess

Pacbio adapters have not been removed by pacbio postprocess  so we will remove them along with PCR primers.
We need to give the sequences of the adapter and primers we want the tool to remove from the sequence.

Below is the list of sequences to remove:
PacBio Adapters:

- Forward Adapter: `GCAGTCGAACATGTAGCTGACTCAGGTCAC`
- Reverse Adapter: `TGGATCACTTGTGCAAGCATCACATCGTAG`

Primer + PacBio Adapters:

- Forward Primer: `GCAGTCGAACATGTAGCTGACTCAGGTCACAGCTTTAATGCGGTAGTTTATC`
- Reverse Primer: `TGGATCACTTGTGCAAGCATCACATCGTAGATCGATGATAAGCTGTCAAAC`

Now, proceed with the following command to launch FROGS preprocess:

```bash

conda activate frogs@4.1.0

five_prime_primer=GCAGTCGAACATGTAGCTGACTCAGGTCACAGCTTTAATGCGGTAGTTTATC
three_prim_primer=GTTTGACAGCTTATCATCGATCTACGATGTGATGCTTGCACAAGTGATCCA


preprocess.py longreads --min-amplicon-size 4000  --max-amplicon-size 5000 \
                          --five-prim-primer $five_prime_primer --three-prim-primer $three_prim_primer -p 6 \
                          --input-R1 fastq/*
```

The tool produce two files:
- `preprocess.fasta` containing the dereplicated reads in fasta format
- `preprocess_counts.tsv` cointaining the count of each reads in each samples. 

## Error analysis 

Tools and package needed to perform this part are listed in the conda_env.yaml file and can be installed within a Conda environment using the following commands:

```bash
conda env create -f conda_env.yaml -n compute_error_rate_env
conda activate compute_error_rate_env 

```

### Aligning the preprocessed sequence to the plasmid reference sequence

For alignment, the plasmid has been linearized at the HindIII restriction site. To align the sequences with the expected amplicon sequence, the fasta file has been manually rearranged to start at the HindIII restriction site : [plasmide_Pbr322_neb_HindIII_cut.fasta](plasmide_Pbr322_neb_HindIII_cut.fasta).

Dereplicated sequences have been aligned to the reference sequence using minimap2 and samtools:

```bash
conda activate compute_error_rate_env 

# reference sequence downloaded manually here: 
ref='plasmide_Pbr322_neb_HindIII_cut.fasta'

minimap2 -ax map-hifi $ref preprocess.fasta --cs  -t 8 | samtools sort -o preprocess_reads_to_plasmide_seq.bam --threads 8

samtools index preprocess_reads_to_plasmide_seq.bam -@ 8

samtools flagstat preprocess_reads_to_plasmide_seq.bam  --threads 8 > preprocess_reads_to_plasmide_seq.flagstat 

```


### Collecting error events from the alignment

Putative errors have been collected using the python script [collect_error_events_in_alignment.py](scripts/collect_error_events_in_alignment.py). This script requires the following Python modules:
- numpy
- pandas
- bamnostic

This script takes the previously computed BAM file as input:

```bash
conda activate compute_error_rate_env

python scripts/collect_error_events_in_alignment.py --bam_file preprocess_reads_to_plasmide_seq.bam --debug

```

The script generates a TSV file named `reads_events_count.tsv`, which describes the type and number of errors observed for each read. The initial lines of the file will resemble this structure:

```tsv
read_id	match	substitution	deletion	insertion	read_len	aln_len	identity	query_coverage	all_error
m64122_210520_233236/109314552/ccs	4339	0	0	0	4421	4339	100.0	98.14521601447636	0
m64122_210520_233236/113509532/ccs	4339	0	0	0	4707	4339	100.0	92.18185680900785	0
m64122_210520_233236/88279866/ccs	4339	0	0	0	4708	4339	100.0	92.16227697536108	0
m64122_210520_233236/51710607/ccs	4338	1	0	1	4734	4340	99.95391705069125	91.67722855935784	2
m64122_210520_233236/73205827/ccs	4338	1	0	1	4422	4340	99.95391705069125	98.14563545906829	2
m64122_210520_233236/96340812/ccs	4338	1	0	2	4699	4341	99.93089149965446	92.38135773568844	3
m64122_210520_233236/108201842/ccs	4338	1	0	3	4423	4342	99.90787655458314	98.16866380284874	4
```
Please be aware that the reads have undergone dereplication using the frogs preprocess tool. This dereplication process needs to be considered when calculating the final error rates. The count of each dereplicated read is provided in the `preprocess_counts.tsv` file, which is generated as output by the frogs preprocess tool.

# Error analysis of plasmid shotgun sequencing

The plasmid has undergone shotgun sequencing without any polymerase amplification, aiming to assess errors introduced solely by the sequencing process. In a manner similar to the amplicon sequences, these reads can be aligned against the plasmid's reference sequence to identify and extract errors.

The Conda environment used for the plasmid's amplicon sequences can also be employed for the following tasks.


## Aligning the Shotgun Reads to the Plasmid Reference Sequence

The same plasmid reference used for the polymerase analysis can also be employed as a reference sequence for this step. The shotgun reads are aligned to the plasmid reference using minimap2, and the resulting alignment is then sorted using samtools.

```bash
conda activate compute_error_rate_env

minimap2 -ax map-hifi $ref fastq/PBR322-proto-shotgun-ligation.fastq.gz --cs  -t 8 | samtools sort -o shotgun_reads_to_plasmide.bam --threads 8

samtools index shotgun_reads_to_plasmide.bam -@ 8

samtools flagstat shotgun_reads_to_plasmide.bam  --threads 8 > shotgun_reads_to_plasmide.flagstat
```

An sbatch script performing the same operations as this code snippet can be found here:  [shotgun_analysis/align_shotgun_reads_to_plasmide.sh](shotgun_analysis/align_shotgun_reads_to_plasmide.sh)


### Collecting error events from the alignment

Similar to the process with the plasmid amplicon reads, errors from the alignments can be extracted using the Python script [collect_error_events_in_alignment.py](scripts/collect_error_events_in_alignment.py). 

```bash
conda activate compute_error_rate_env

python collect_error_events_in_alignment.py --bam_file shotgun_reads_to_plasmide.bam --debug

```

# Compute error rates and plot them 

The error rates for both the plasmids with polymerase amplification and the shotgun sequencing have been computed within this Jupyter notebook: [plot_error_rates.ipynb](plot_error_rates.ipynb).

The notebook produce this final summary table `error_rate_accustart_hotstart_shotgun.tsv`: 


| polymerase | substitution_mean_rate | deletion_mean_rate     | deletion_std_rate      | insertion_mean_rate   | insertion_std_rate     | all_error_mean_rate   | all_error_std_rate     |
|------------|------------------------|------------------------|------------------------|-----------------------|------------------------|-----------------------|------------------------|
| Accustart  | 0.0017034985061008922  | 0.00037539763877669774 | 1.2342878715049305e-05 | 0.0010065167566872494 | 4.5255706875587496e-05 | 0.0030854129015648393 | 4.5295355303927724e-05 |
| Hotstart   | 0.0002623945622201094  | 0.0002971164935496605  | 3.340757540331315e-06  | 0.0008982427957152317 | 2.006612294022626e-05  | 0.0014577538514850014 | 1.7921511685849816e-05 |
| shotgun    | 8.509540236488563e-05  | 0.00027217389668779683 |                        | 0.0008751111282520822 |                        | 0.0012323804273047646 |                        |

