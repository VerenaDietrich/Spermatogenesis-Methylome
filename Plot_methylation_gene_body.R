## ---------------------------
##
## Script name: Plot_methylation_gene_body.R
##
## Purpose of script: Plot the methylation levels over the gene body
##
## Author: Verena Dietrich
##
## Date Created: 2023-07-11
## Last Updated: 2023-10-11
##
## Email: verena.dietrich@uni-muenster.de
##
## R version 4.2.1 (2022-06-23)
## ---------------------------

##############################
# 0 - Load libraries
##############################
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install("rtracklayer")
library(rtracklayer)
library(stringr)
library(ggplot2)

##############################
# 1 - Import gtf-annotation-file (provided by wg-blimp)
##############################
annotation_df <- as.data.frame(rtracklayer::import('/path/to/spermatogenesis-methylome/annotation/genes.hg38.gtf/gencode.v34.annotation.gtf'))
# prepare annotation-data frame for further computation
transcript_sites <- annotation_df[annotation_df$type == "transcript",][c("seqnames", "start", "end")]
transcript_sites$seqnames <-gsub("chr","", transcript_sites$seqnames)

##############################
# 2 - Define number of bins and base pairs up and downstream
##############################

bins_before <- 5
bins_after <- 5
intervall_surr <- seq(0,5000, length.out = bins_before+1) # 5000 bp for and after the transcript divided into the number of bins that should be displayed before and after
number_poitions_plot <- 50
## For every transcript has its bin size, so that every transcript has the same number of bins
transcript_sites$bin_width <- (transcript_sites$end - transcript_sites$start)/number_poitions_plot
boundaries_bin = data.frame(matrix(nrow = nrow(transcript_sites), ncol = 0))

##############################
# 3 - Define boundaries of bins according to transcript length
##############################

## Get for each transcript the start and end of the bins
#Intervals before
for (i in length(intervall_surr):1){
  boundaries_bin <- cbind(boundaries_bin, transcript_sites$start - intervall_surr[i])
}
#Intervals in between
for (i in 1:(number_poitions_plot-1)){
  boundaries_bin <- cbind(boundaries_bin, transcript_sites$start + transcript_sites$bin_width * i)
}
boundaries_bin <- round(boundaries_bin, digits = 0)
#Intervals after
for (i in 1:length(intervall_surr)){
  boundaries_bin <- cbind(boundaries_bin, transcript_sites$end + intervall_surr[i])
}

##############################
# 4 - Import raw data
##############################

file_path <- "/path/to/bed-graphs/"
files <- c("179937_CpG.bedGraph", "179946_CpG.bedGraph", "179954_CpG.bedGraph", "179940_CpG.bedGraph", "179948_CpG.bedGraph", "179956_CpG.bedGraph",
           "179942_CpG.bedGraph", "179950_CpG.bedGraph", "179958_CpG.bedGraph", "179944_CpG.bedGraph", "179952_CpG.bedGraph", "179960_CpG.bedGraph",
           "196372_CpG.bedGraph", "196374_CpG.bedGraph", "196376_CpG.bedGraph", "196378_CpG.bedGraph", "196380_CpG.bedGraph", "196382_CpG.bedGraph"
)

##############################
# 5 - define data frame for results
##############################

methylation_combined <- data.frame(
  Bin = integer(),
  Methylation = double(),
  Sample = character()
) 

##############################
# 6 - Compute mean methylation level for each bin for each sample
##############################

for (file in files){
  print(file)
  methylation <- import.bedGraph(paste(file_path, file, sep = "/"))
  for (bin in 1:(ncol(boundaries_bin)-1)){
    print(paste0("bin: ", bin)) 
    ### add for gene in transcript_sites -> min coverage of 5, mean for each gene -> add to methylaton_level list  
    GRange_bin <- GRanges(seqname = transcript_sites$seqnames, ranges = IRanges(boundaries_bin[,bin], boundaries_bin[,bin+1]))
    overlap <- subsetByOverlaps(methylation, GRange_bin)
    overlap_mincov <- overlap[(overlap$NA. + overlap$NA..1) >=5]
    methylation_level <- overlap_mincov$score/100 # for one bin with all genes for one sample
    
    # subtable
    methylation_sub <- data.frame(
      Bin = bin,
      Methylation = mean(methylation_level),
      Sample = str_split(file, "_")[[1]][1]
    )
    methylation_combined <- rbind(methylation_combined, methylation_sub)
  }
}

##############################
# 7 - Add sample information 
##############################
# Modify according to your needs
methylation_combined['Group'] <- NA
methylation_combined[methylation_combined$Sample %in% c("179937", "179946", "179954"),]$Group <- "Control_Undiff"
methylation_combined[methylation_combined$Sample %in% c("179940", "179948", "179956"),]$Group <- "Control_Diff"
methylation_combined[methylation_combined$Sample %in% c("179942", "179950", "179958"),]$Group <- "Control_4N"
methylation_combined[methylation_combined$Sample %in% c("179944", "179952", "179960"),]$Group <- "Control_1N"
methylation_combined[methylation_combined$Sample %in% c("196372", "196378"),]$Group <- "Crypto_Undiff"
methylation_combined[methylation_combined$Sample %in% c("196374", "196380"),]$Group <- "Crypto_Diff"
methylation_combined[methylation_combined$Sample %in% c("196376", "196382"),]$Group <- "Crypto_4N"

##############################
# 4 - Save data for further plots  
##############################

#write.csv(methylation_combined, "/path/to/results/result_TSS-TES_minCov.csv", row.names=FALSE)
#methylation_combined <- read.csv('/path/to/results/result_TSS-TES_minCov.csv', header=TRUE)

##############################
# 5 - Plots 
##############################

methylation_combined$Sample <- factor(methylation_combined$Sample)
methylation_combined$Group <- factor(methylation_combined$Group, levels = c("Control_Undiff", "Control_Diff", "Control_4N", "Control_1N", "Crypto_Undiff", "Crypto_Diff", "Crypto_4N"))
methylation_combined$Patient  <- factor(methylation_combined$Patient, levels = c("Control", "Crypto"))
methylation_combined$State <- factor(methylation_combined$State, levels = c("Undiff", "Diff", "4N", "1N"))
methylation_combined$Group <- factor(methylation_combined$Group, levels = c( "Control_Undiff", "Control_Diff", "Control_4N", "Control_1N", "Crypto_Undiff", "Crypto_Diff", "Crypto_4N" ), ordered = TRUE)

## Plot all Samples 
ggplot(data  = methylation_combined, mapping = aes( y=Methylation, x=Bin, colour=Sample))+
  geom_line(linewidth = 0.5) +
  scale_color_manual(values=c("#80146E", "#478EC1", "#7ED5B8", "#F5F2D8", "#7E4E90ff" , "#5a829e", "#69B3A2", "#80146E", "#478EC1", "#7ED5B8", "#F5F2D8", "#7E4E90ff" , "#5a829e", "#69B3A2", "#80146E", "#478EC1", "#7ED5B8", "#F5F2D8")) +
  geom_vline(xintercept = bins_before + 0.5, linetype="dotted", 
                  color = "blue", linewidth=0.5) +
  geom_vline(xintercept = number_poitions_plot + bins_before + 0.5, linetype="dotted", 
             color = "blue", linewidth=0.5) +
  ylim(0, 1)+
  theme_classic()

## Plot groupwise
ggplot(data = methylation_combined, mapping = aes( y=Methylation, x=Bin))+ 
  stat_summary(aes(group=Group, colour=Group), fun=mean, geom="line")+
  scale_color_manual(values=c("#80146E", "#478EC1", "#7ED5B8", "#F5F2D8", "#7E4E90ff" , "#5a829e", "#69B3A2")) +
  geom_vline(xintercept = bins_before + 0.5, linetype="dotted", 
             color = "blue", linewidth=0.5) +
  geom_vline(xintercept = number_poitions_plot + bins_before + 0.5, linetype="dotted", 
             color = "blue", linewidth=0.5) +
  ylim(0, 1)+
  theme_classic()

## Plot groupwise only Control 
ggplot(data = methylation_combined[methylation_combined$Patient == "Control",], mapping = aes( y=Methylation, x=Bin))+ 
  stat_summary(aes(group=Group, colour=Group), fun=mean, geom="line")+
  scale_color_manual(values=c("#80146E", "#478EC1", "#7ED5B8", "#F5F2D8", "#7E4E90ff" , "#5a829e", "#69B3A2")) +
  geom_vline(xintercept = bins_before + 0.5, linetype="dotted", 
             color = "blue", linewidth=0.5) +
  geom_vline(xintercept = number_poitions_plot + bins_before + 0.5, linetype="dotted", 
             color = "blue", linewidth=0.5) +
  ylim(0, 1)+
  theme_classic()

#ggsave("/path/to/results/TSS-TES-Plot/crypto/Plot_sample_minCov.png", width = 7.5, height = 6, dpi = 300)


