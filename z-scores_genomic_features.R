## ---------------------------
##
## Script name: z-scores_genomic_features.R
##
## Purpose of script: Compute the z-score for the enrichment of DMRs in genomic regions (permutation test) 
##
## Author: Verena Dietrich
##
## Date Created: 2023-05-10
## Updated: 2023-05-23
## Last Update: 2023-09-29
##
## Email: verena.dietrich@uni-muenster.de
##
## R version 4.3.0 (2022-06-23)
## ---------------------------

##############################
# 0 - Load libraries
##############################
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomicRanges")
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("annotatr")
BiocManager::install("regioneR")

library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(annotatr)
library(regioneR)
library(ggplot2)

set.seed(123)

##############################
# 1 - Load peaks of interest (DMRs)
##############################

files = c("/path/to/DMRs/peaks_4N_1N.txt",
          "/path/to/DMRs/peaks_Diff_1N.txt",
          "/path/to/DMRs/peaks_Diff_4N.txt",
          "/path/to/DMRs/peaks_Undiff_1N.txt",
          "/path/to/DMRs/peaks_Undiff_4N.txt",
          "/path/to/DMRs/peaks_Undiff_diff.txt",
          "/path/to/DMRs/Final_DMRs_4N_CZvsCTR_15range.txt",
          "/path/to/DMRs/Final_DMRs_Diff_CZvsCTR_15range.txt",
          "/path/to/DMRs/Final_DMRs_Undiff_CZvsCTR_15range.txt")

##############################
# 2 - Set regions of interest
##############################

annots <- c('hg38_cpgs', 'hg38_basicgenes', 'hg38_genes_intergenic','hg38_genes_intronexonboundaries',
           'hg38_genes_promoters',
           'hg38_genes_5UTRs',
           'hg38_genes_exons',
           'hg38_genes_introns',
           'hg38_genes_3UTRs', 
           'hg38_genes_1to5kb', 
           'hg38_genes_cds',
           'hg38_genes_firstexons', 
           'hg38_genes_exonintronboundaries')

# Annotations
annotations <- build_annotations(genome = 'hg38', annotations = annots)
annotations_unique <- unique(annotations$type)

##############################
# 3 - Perform permutation test with peaks
##############################

df = data.frame() # Empty data frame for results

# Loop over all files with peaks of interest
for (file in files){
  # Load peaks of interest as GRange object
  comparison <- read.table(file , header = TRUE)
  comaprison_gRanges <- GRanges(seqnames=comparison$seqnames,ranges=IRanges(comparison$start, comparison$stop))
  # Loop over all features in the annotation
  for (feature in annotations_unique){
    feature_object = annotations[annotations$type == feature,]
    feature_object_unique <- unique(feature_object)
    # Permutation test
    pt <- overlapPermTest(A=comaprison_gRanges, B=feature_object_unique , ntimes=10000, genome = "hg38")
    print(paste0(file, ": z_score for ", feature, ": ", pt$numOverlaps$zscore, "p_value: ", pt$numOverlaps$pval))
    df = rbind(c(file, feature, pt$numOverlaps$zscore, pt$numOverlaps$pval),df )
  }
}

##############################
# 4 - Save results
##############################

colnames(df) <- c("file", "feature", "z_score", "p_value")
write.table(df, "/path/to/results/results_z-score_features_cryptoControl_seed123.csv", sep = ";", row.names = FALSE)

