## ---------------------------
##
## Script name: z-scores_repeats.R
##
## Purpose of script: Compute the z-score for the enrichment of DMRs in repeats (permutation test) 
##
## Author: Verena Dietrich
##
## Date Created: 2023-05-04
## Updated: 2023-05-25
## Last Update: 2023-10-11
##
## Email: verena.dietrich@uni-muenster.de
##
## R version 4.3.0 (2023-04-21)
## ---------------------------

##############################
# 0 - Load libraries
##############################
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomicRanges")
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
BiocManager::install("annotatr")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("regioneR")

library(GenomicRanges)
library(regioneR)
library(ggplot2)

set.seed(123)

##############################
# 1 - Load peaks of interest (DMRs)
##############################

files = c(#"/path/to/DMRs/peaks_4N_1N.txt",
          #"/path/to/DMRs/peaks_Diff_1N.txt",
          #"/path/to/DMRs/peaks_Diff_4N.txt",
          #"/path/to/DMRs/peaks_Undiff_1N.txt",
          #"/path/to/DMRs/peaks_Undiff_4N.txt",
          #"/path/to/DMRs/peaks_Undiff_diff.txt",
          #"/path/to/DMRs/Final_DMRs_4N_CZvsCTR_15range.txt",
          "/path/to/DMRs/Final_DMRs_Diff_CZvsCTR_15range.txt",
          "/path/to/DMRs/Final_DMRs_Undiff_CZvsCTR_15range.txt")

##############################
# 2 - Set regions of interest
##############################

repeatmasker <- read.csv('/path/to/spermatogenesis-methylome/annotation/repeat-masker-hg38.csv', header=TRUE, sep = "\t" )
repeatmasker_sub <- subset(repeatmasker, select = c("genoName", "genoStart", "genoEnd", "repClass")) #, "repFamily" 
annotations_unique <- unique(repeatmasker_sub$repClass)

##############################
# 3 - Perform permutation test with peaks
##############################

df = data.frame() # Empty data frame for results

# Loop over all files with peaks of interest
for (file in files){
  # Load peaks of interest as GRange object
  comparison <- read.table(file , header = TRUE)
  comaprison_gRanges <- GRanges(seqnames=comparison$seqnames,ranges=IRanges(comparison$start, comparison$end))
  # Loop over all features in the annotation
  for (feature in annotations_unique){
    feature_object = repeatmasker_sub[repeatmasker_sub$repClass == feature,]
    feature_object_gRanges <- GRanges(seqnames=feature_object$genoName,ranges=IRanges(feature_object$genoStart, feature_object$genoEnd))
    # Permutation test
    pt <- overlapPermTest(A=comaprison_gRanges, B=feature_object_gRanges , ntimes=10000, genome = "hg38")
    print(paste0(file, ": z_score for ", feature, ": ", pt$numOverlaps$zscore, "p_value: ", pt$numOverlaps$pval))
    df = rbind(c(file, feature, pt$numOverlaps$zscore, pt$numOverlaps$pval),df )
  }
}

##############################
# 4 - Save results
##############################
colnames(df) <- c("file", "feature", "z_score", "p_value")
write.table(df, "/path/to/results/results_z-score_repeats_cryptoControl_seed123.txt", sep = ";", row.names = FALSE)






