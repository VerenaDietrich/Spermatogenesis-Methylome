## Date Created: 2023-10-09
## Updated: 2023-10-09
## Last Update: 2023-10-09
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

BiocManager::install("rtracklayer")
BiocManager::install("reshape")
library(rtracklayer)
library(stringr)
library(ggplot2)
library(reshape)
library(GenomicRanges)
# for split plot
install.packages("devtools")
devtools::install_github("psyteachr/introdataviz")
library(introdataviz)

##############################
# 1 - Load files
##############################

# import annotation file 
gtf <- rtracklayer::import('/path/to/spermatogenesis-methylome/annotation/genes.hg38.gtf/gencode.v34.annotation.gtf')
annotation_df <- as.data.frame(gtf)
types <- c("gene", "transcript", "exon", "CDS", "UTR") # "start_codon", "stop_codon", "Selenocysteine"

# import input-files
file_path <- "/path/to/bed-graphs/"
# Modify to your needs
files <- c("179937_CpG.bedGraph", "179946_CpG.bedGraph", "179954_CpG.bedGraph", "179940_CpG.bedGraph", "179948_CpG.bedGraph", "179956_CpG.bedGraph",
           "179942_CpG.bedGraph", "179950_CpG.bedGraph", "179958_CpG.bedGraph", "179944_CpG.bedGraph", "179952_CpG.bedGraph", "179960_CpG.bedGraph",
           "196372_CpG.bedGraph", "196374_CpG.bedGraph", "196376_CpG.bedGraph", "196378_CpG.bedGraph", "196380_CpG.bedGraph", "196382_CpG.bedGraph"
)

##############################
# 2 - Create data frame for result
##############################

methylation_combined <- data.frame(
  Sample = character(), 
  Methylation = double(), 
  Type = character()
) 

for (file in files){
  print(file)
  methylation <- import.bedGraph(paste(file_path, file, sep = "/"))
  for (type in types){
    
    print(type)
    type_sites <- annotation_df[annotation_df$type == type,][c("seqnames", "start", "end")]
    type_sites$seqnames <-gsub("chr","", type_sites$seqnames)
    type_sites_object <- GRanges(seqnames=type_sites$seqnames,
                                 ranges=IRanges(start = type_sites$start, end = type_sites$end))
    overlap <- subsetByOverlaps(methylation, type_sites_object)
    # min coverage of 5
    overlap_mincov <- overlap[(overlap$NA. + overlap$NA..1) >=5]
    methylation_level <- overlap_mincov$score/100
    
    # subtable
    methylation_sub <- data.frame(
      Methylation = methylation_level,
      Sample = str_split(file, "_")[[1]][1],
      Type = type
    )
    methylation_combined <- rbind(methylation_combined, methylation_sub)
  }
}

##############################
# 3 - Add sample information 
##############################

# Modify to your needs
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

#write.csv(methylation_combined, "/path/to/results/result_gtf_annotation_cryptoANDcontrol_minCov.csv", row.names=FALSE)

#methylation_combined <- read.csv('/path/to/results/result_gtf_annotation_cryptoANDcontrol_minCov.csv', header=TRUE)

##############################
# 5 - Plots 
##############################

methylation_combined$Group <- factor(methylation_combined$Group, levels = c("Control_Undiff", "Control_Diff", "Control_4N", "Control_1N", "Crypto_Undiff", "Crypto_Diff", "Crypto_4N"))
# Patient groups
methylation_combined$Patient <- gsub("_.*", "\\1", methylation_combined$Group)
methylation_combined$Patient  <- factor(methylation_combined$Patient, levels = c("Control", "Crypto"))
# States of cell development
methylation_combined$State <- gsub(".*\\_", "\\1", methylation_combined$Group)
methylation_combined$State <- factor(methylation_combined$State, levels = c("Undiff", "Diff", "4N", "1N"))

####### Plot all types
plot <- ggplot(data = methylation_combined, mapping = aes(y=Methylation, x=Group, fill = Group)) + 
  geom_violin(scale = "count")+ 
  scale_fill_manual(values=c("#80146E", "#478EC1","#7ED5B8",  "#F5F2D8", "#7E4E90ff", "#5a829e", "#69B3A2"))+
  facet_wrap(~Type,ncol = 2)+
  theme_classic() 
#ggsave("/path/to/results/gtf_Types.png", width = 10, height = 10, dpi = 300)


###### Violine Split Plot
plot <- ggplot(data = methylation_combined, aes(x = State, y = Methylation, fill = Patient)) +
  geom_split_violin(scale = "count") + # oder "area" oder "with"
  facet_wrap(~Type,ncol = 2)+
  scale_fill_manual(values=c("#80146E", "#478EC1","#7ED5B8",  "#F5F2D8", "#7E4E90ff"))

#ggsave("/path/to/results/gtf_Splitplot.png", width = 14, height = 6, dpi = 300)

###### All samples separate
methylation_combined$Sample <- as.character(methylation_combined$Sample) 
plot_gtf <- ggplot(methylation_combined, aes(x = Sample, y = Methylation)) + 
  geom_violin(fill = "blue", color = "black", alpha = 0.7) + 
  #geom_boxplot(width = 0.1) + 
  theme_minimal() + 
  labs(title = "Verteilung der Methylierung pro Sample", x = "Sample", y = "Methylierung") +
  theme(axis.text.x = element_text(angle = 25, vjust = 0.5)) +
  facet_wrap(~ Type, ncol = length(unique(methylation_combined$Type))) +
  ggtitle("GTF")

#ggsave("/path/to/results/gtf_all_samples.jpg", plot = plot_gtf, width = 15, height = 10, units = "cm", dpi = 300)

