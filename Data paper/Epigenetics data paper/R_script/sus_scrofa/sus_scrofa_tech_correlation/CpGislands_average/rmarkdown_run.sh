#!/bin/bash
#SBATCH -p workq
#SBATCH -J sscrofadtpaper_wholecpg
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
#SBATCH --mem=50G
#SBATCH --cpus-per-task=2


module purge
module load system/pandoc-2.1.3
module load system/Miniconda3-4.7.10

source activate r-env

Rscript -e 'rmarkdown::render("Compare_tech_heatmap_datapaper_no_filter_cpgisland_average.Rmd")'
