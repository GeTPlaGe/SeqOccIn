# SeqOccIn - Variants detection in the 29 maize lines

The structural variants and SNPs were detected using the following pipeline:  
   https://github.com/SeqOccin-SV/SeqOccinVariants



The PCA plot was made using the following R code.


```R
library(SNPRelate)
library(ggplot2)
library(tidyverse)

vcf.fn <- "snps.vcf.gz"
gds.fn <- "data-2strains.gds"
snpgdsVCF2GDS(vcf.fn, gds.fn, method="biallelic.only")

genofile <- snpgdsOpen(gds.fn)
pca <- snpgdsPCA(genofile, autosome.only=TRUE)
pve <- pca$varprop

pc.df <- data.frame(sample = pca$sample.id, PC1 = pca$eigenvect[,1], PC2 = pca$eigenvect[,2])

# Load metadata (make sure the SampleID matches the sample names in your GDS file)
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
