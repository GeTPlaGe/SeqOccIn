# SeqOccIn - Benchmarking of sequencing technologies for CpG methylation calling in quail (Coturnix japonica) and pig (Sus scrofa)

The publication about these data is available here : [doi:](https://www.nature.com/articles/)

All raw reads files produced from this project are available at EMBL/EBI ENA under the BioProject PRJEB55065 for quail (https://www.ebi.ac.uk/ena/browser/view/PRJEB55065) and PRJEB55066 for pig (https://www.ebi.ac.uk/ena/browser/view/PRJEB55066).

The result files from methylation calling per technology can be accessed at Recherche Data Gouv using [DOI:10.57745/PIHOIA] (https://doi.org/10.57745/PIHOIA) for quail and [DOI:10.57745/EMMOSV] (https://doi.org/10.57745/EMMOSV) for pig.

Used reference assemblies were Coturnix japonica version 2.0 available on https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/577/835/GCA_001577835.1_Coturnix_japonica_2.0/ and Sus scrofa version 11.1 available on https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/003/025/GCA_000003025.6_Sscrofa11.1/.

## CODE USAGE

The following scripts were run to produce the plots found in [seqoccin DNA methylation datapaper](https://www.nature.com/articles/).
For both samples (Coturnix Japonica, Sus scrofa), we generated plots in order to describe the quality of our CpG methylation detection with WGBS, EM-seq, ONT and PacBio.

All scripts should be adapted to run in your environment. Paths to result files and reference assemblies should be modified.
