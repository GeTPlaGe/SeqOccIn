#!/bin/bash
#SBATCH -p testq
#SBATCH -J rmarkdown_generator
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
#SBATCH --mem=60G
#SBATCH --cpus-per-task=8


module purge
module load system/pandoc-2.1.3
module load system/Miniconda3-4.7.10

source activate r-env

Rscript -e 'rmarkdown::render("dtpaper_Sscrofa_cumulative_coverageplot_alltech.Rmd")'
