## Date Created: 2023-10-09
## Updated: 2023-09-12
## Last Update: 2023-10-10
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
library(tidyverse)

##############################
# 1 - Load files
##############################

repeatmasker <- read.csv('/path/to/spermatogenesis-methylome/annotation/repeat-masker-hg38.csv', header=TRUE, sep = "\t" )
repeatmasker_sub <- subset(repeatmasker, select = c("genoName", "genoStart", "genoEnd", "repName")) #, "repFamily" 
annotation_df <- as.data.frame(repeatmasker_sub)
types <- c("L1HS", "L1PA2", "L1PA3", "L1PA4", "L1PA5", "SVA_D", "SVA_F", "HERVH-int", "L1M7")

file_path <- "/path/to/bed-graphs/"
files <- c("179937_CpG.bedGraph", "179946_CpG.bedGraph", "179954_CpG.bedGraph", "179940_CpG.bedGraph", "179948_CpG.bedGraph", "179956_CpG.bedGraph",
           "179942_CpG.bedGraph", "179950_CpG.bedGraph", "179958_CpG.bedGraph", "179944_CpG.bedGraph", "179952_CpG.bedGraph", "179960_CpG.bedGraph",
           "196372_CpG.bedGraph", "196374_CpG.bedGraph", "196376_CpG.bedGraph", "196378_CpG.bedGraph", "196380_CpG.bedGraph", "196382_CpG.bedGraph"
)

##############################
# 2 - Create data frame for result and statistics
##############################

methylation_combined <- data.frame(
  Sample = character(), 
  Methylation = double(), 
  Type = character()
) 

df_statsitc <- data.frame(
  sample = character(),
  type = character(),
  length = integer(),
  overlap = integer(),
  overlap_minCov= integer()
)

for (file in files){
  print(file)
  methylation <- import.bedGraph(paste(file_path, file, sep = "/"))
  for (type in types){
    print(type)
    type_sites <- annotation_df[annotation_df$repName == type,][c("genoName", "genoStart", "genoEnd")]
    type_sites$genoName <-gsub("chr","", type_sites$genoName)
    type_sites_object <- GRanges(seqnames=type_sites$genoName,
                                 ranges=IRanges(start = type_sites$genoStart, end = type_sites$genoEnd))
    overlap <- subsetByOverlaps(methylation, type_sites_object)
    overlap_mincov <- overlap[(overlap$NA. + overlap$NA..1) >=5]
    methylation_level <- overlap_mincov$score/100
    df_statsitc_new <- data.frame(
      sample = file,
      type = type,
      length = sum(lengths(type_sites_object)),
      overlap = sum(lengths(overlap)),
      overlap_minCov= sum(lengths(overlap_mincov)))
    df_statsitc <- rbind(df_statsitc,df_statsitc_new)
    
    
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
# 3 - Statistics
##############################


### Statistics
df_statsitc$overlap_ratio <- df_statsitc$overlap/df_statsitc$length
df_statsitc$overlap_minCov_ratio <- df_statsitc$overlap_minCov/df_statsitc$overlap
#write.csv(df_statsitc, "/path/to/results/statistics_retroposons/result_repeatMasker_transposons_minCov_statistic.csv", row.names=FALSE)

# Histogramm
png(file="/path/to/results/statistics_retroposons/transposons_overlap_hist.png",
    width=600, height=350)
hist <- hist(df_statsitc$overlap_ratio, main="Overlap der Regionen mit den Methylierungsdaten", xlab="Anteil des Overlaps", col="grey", border="black", breaks = 20)
dev.off()

# Boxplot
mincover_plot <- ggplot(data = df_statsitc, aes(x = type, y = overlap_minCov_ratio, fill = type)) +
  geom_boxplot() +
  labs(title="Ratio fullfilling min coverage of 5", x="Retrotransposons", y="Overlap ratio")
ggsave("/path/to/results/statistics_retroposons/transposons_min_cov_boxplot.png", plot = mincover_plot, width = 7.5, height = 6, dpi = 300)

# Jitterplot
mincover_plot_jitter <- ggplot(data = df_statsitc, aes(x = type, y = overlap_minCov_ratio, color = type)) +
  geom_jitter(width=0.2) +
  #geom_text(aes(label=sample), vjust=-1) +
  labs(title="Ratio fullfilling min coverage of 5", x="Retrotransposons", y="Overlap ratio")
ggsave("/path/to/results/statistics_retroposons/transposons_min_cov_jitterplot.png", plot = mincover_plot_jitter, width = 7.5, height = 6, dpi = 300)

##############################
# 4 - Add sample information 
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
# 5 - Save data for further plots  
##############################

#write.csv(methylation_combined, "/path/to/results/result_repeatMasker_transposons_minCov.csv", row.names=FALSE)
#methylation_combined <- read.csv('/path/to/results/result_repeatMasker_transposons_minCov.csv', header=TRUE)

##############################
# 6 - Plots 
##############################

methylation_combined$Group <- factor(methylation_combined$Group, levels = c("Control_Undiff", "Control_Diff", "Control_4N", "Control_1N", "Crypto_Undiff", "Crypto_Diff", "Crypto_4N"))
methylation_combined$Patient <- gsub("_.*", "\\1", methylation_combined$Group)
methylation_combined$State <- gsub(".*\\_", "\\1", methylation_combined$Group)
methylation_combined$Patient  <- factor(methylation_combined$Patient, levels = c("Control", "Crypto"))
methylation_combined$State <- factor(methylation_combined$State, levels = c("Undiff", "Diff", "4N", "1N"))
methylation_combined$Type <- factor(methylation_combined$Type, levels = c("SVA_D", "SVA_F", "L1HS", "L1PA2", "L1PA3", "L1PA4", "L1PA5", "HERVH-int", "L1M7"))

### Split violine plot
methylation_combined_splitplot <- methylation_combined[methylation_combined$State %in% c("Undiff", "Diff", "4N"),]
library(introdataviz)
plot_1 <- ggplot(data = methylation_combined_splitplot, aes(x = State, y = Methylation, fill = Patient)) +
  geom_split_violin(scale = "area", linewidth = 0.3) + # oder "area" oder "width" oder "count"
  facet_wrap(~Type,ncol = 3)+
  theme_classic() +
  theme(strip.text.x= element_text(color="black",size = 7), strip.background = element_rect(
    fill = "white", color = "black",linewidth = 0.3), 
    axis.title.y = element_text(colour = "black", size = 7),
    legend.position="none",
    axis.title.x =element_blank(),axis.text.x=element_text(colour = "black", size = 7, angle = 45, vjust = 1, hjust=1),axis.text.y=element_text(colour = "black", size = 7),axis.line = element_line(linewidth = 0.3, colour = "black"))+
  scale_fill_manual(values=c("aquamarine3", "darkorchid4"))

#ggsave("/path/to/results/Repeat_retroposons_Splitplot_minCov_white_linewidth_area.png", plot = plot_1, width = 8.7, height = 8.7, units = "cm", dpi = 300)


## Control
methylation_combined_interest <- methylation_combined[methylation_combined$Patient == "Control",]
methylation_combined_interest$State <- factor(methylation_combined_interest$State, levels = c( "Undiff", "Diff",   "4N" ,    "1N"  ))
methylation_combined_interest$Type <- factor(methylation_combined_interest$Type, levels = c("L1HS", "L1PA2", "L1PA3", "L1PA4", "L1PA5", "SVA_D", "SVA_F", "HERVH-int", "L1M7"))

g1 <- ggplot(data = methylation_combined_interest, mapping = aes(y=Methylation, x=State, fill = State)) + 
  geom_violin(scale = "count",linewidth = 0.1, colour = "black")+
  scale_fill_manual(values=c("#001889","#91008D", "#D24E71","#EDA200"))+facet_wrap(~factor(Type,levels=c("SVA_D", "SVA_F", "L1HS", "L1PA2", "L1PA3", "L1PA4", "L1PA5", "HERVH-int", "L1M7")),ncol=3, nrow=3)+theme_classic()+ 
  theme(strip.text.x= element_text(color="black",size = 7), strip.background = element_rect(
    fill = "grey", color = "black",linewidth = 0.3),axis.title.y = element_text(colour = "black", size = 7),
    #legend.text = element_text(size = 5), legend.title = element_text(size = 7), legend.key.size = unit(0.5, 'cm'),
    legend.position="none",
    axis.title.x =element_blank(),axis.text.x=element_text(colour = "black", size = 7, angle = 45, vjust = 1, hjust=1),axis.text.y=element_text(colour = "black", size = 7),axis.line = element_line(linewidth = 0.3, colour = "black"))

#ggsave("/path/to/results/Repeat_retroposons_Control_minCov_grey.png", plot = g1, height = 10, width =10, units = "cm", dpi = 300)
