# SeqOccIn - Variants detection in the 29 maize lines

The structural variants and SNPs were detected using the following pipeline:  
   https://github.com/SeqOccin-SV/SeqOccinVariants

## PCA plot

The PCA plot in the paper was made using the following R code.

```R
library(SNPRelate)
library(ggplot2)
library(tidyverse)

# snps.vcf.gz is a vcf file made of 1M randomly selected SNPs among variants
# for which an alternative allele appears at least in two samples
vcf.fn <- "snps.vcf.gz"
gds.fn <- "data.gds"
snpgdsVCF2GDS(vcf.fn, gds.fn, method="biallelic.only")

genofile <- snpgdsOpen(gds.fn)
# Compute PCA using snpgdsPCA from SNPRelate library
pca <- snpgdsPCA(genofile, autosome.only=TRUE)
# Extract percent of variance explained for each component
pve <- pca$varprop

pc.df <- data.frame(sample = pca$sample.id, PC1 = pca$eigenvect[,1], PC2 = pca$eigenvect[,2])

# Load metadata provinding population associated to each sample
metadata <- read.csv("../Populations.csv", sep="\t")

# Merge PCA results with metadata (based on SampleID)
pca_data <- merge(pc.df, metadata, by.x = "sample", by.y = "sample")

cbPallete <- c("Stiff Stalk Synthetic"= "#3399ff", 
               "Iodent"="#339900", 
               "Lancaster"="#00cc33", 
               "Lancaster / Iodent"="#999933", 
               "European Flint"="#ffcc33", 
               "Italian Flint" ="#ff9900", 
               "Northern Flint" = "#cc6633", 
               "Tropical Highland" ="#ff3366", 
               "Tropical Spanish" ="#cc0066" )

var_explained <- round(pca$varprop[1:2] * 100, 2)
x_lab <- paste0("PC1 (", var_explained[1], "%)")
y_lab <- paste0("PC2 (", var_explained[2], "%)")

# Plot the PCA
options(repr.plot.width = 8, repr.plot.height = 6)
pca_data %>% rename(Subgroup = subgroup) %>%
  ggplot(aes(x = PC1, y = PC2, color = Subgroup)) +
  geom_point(size = 3) +
  geom_text(aes(label = sample), vjust = -1, show.legend = FALSE) +
  xlab(x_lab) +
  ylab(y_lab) +
  theme_minimal() +
  scale_color_manual(values = cbPallete) 


```
