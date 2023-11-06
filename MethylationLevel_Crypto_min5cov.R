# Function: extracting of the methylation level for specific regions
# Input: bedGraph for samples of which the methylation level should be computed
# Otput: .csv-file with methylation level of the specified regions

# Author: Verena Dietrich, Jun 2022

BiocManager::install("plyranges")
BiocManager::install("rtracklayer")
install.packages("stringr")
library(stringr)
library(rtracklayer)

# To be modified
source <- "imprinted genes"
#source <- "3000 genes"
file_path <- "/path/to/bed-graphs/"
files <- c("BloodSperm/60_Hmp01_blood_young_CpG.bedGraph", 
           "BloodSperm/60_Hmp02_blood_old_CpG.bedGraph", 
           "BloodSperm/60_Hmp01_sperm_young_CpG.bedGraph", 
           "BloodSperm/60_Hmp02_sperm_old_CpG.bedGraph",
           "179937_CpG.bedGraph", "179946_CpG.bedGraph", "179954_CpG.bedGraph", "179940_CpG.bedGraph", "179948_CpG.bedGraph", "179956_CpG.bedGraph",
           "179942_CpG.bedGraph", "179950_CpG.bedGraph", "179958_CpG.bedGraph", "179944_CpG.bedGraph", "179952_CpG.bedGraph", "179960_CpG.bedGraph",
           "196372_CpG.bedGraph", "196374_CpG.bedGraph", "196376_CpG.bedGraph", "196378_CpG.bedGraph", "196380_CpG.bedGraph", "196382_CpG.bedGraph",
           "contaminated/196365_CpG.bedGraph", "contaminated/196368_CpG.bedGraph", "contaminated/196370_CpG.bedGraph"
           )

# Sources:
# Chose between the collection of genes you want the methylation level for
# imprinted genes, 50 genes (https://static-content.springer.com/esm/art%3A10.1186%2Fs13148-020-00854-0/MediaObjects/13148_2020_854_MOESM2_ESM.pdf,
#                            https://clinicalepigeneticsjournal.biomedcentral.com/articles/10.1186/s13148-020-00854-0#MOESM1)
# 3000 genes (https://d-nb.info/1213726247/34)

## Looping over all regions

if (source == "imprinted genes") {
  regions_file <- read.csv(file = "/path/to/spermatogenesis-methylome/annotation/Imprinted_genes_location.csv")
  DMR_names <- regions_file$Imprinted.DMR.name
  regions <- regions_file$Location_.hg38.
  df <- data.frame(DMR_Names = DMR_names,
                   Region = regions,
                   Number_CpGs_publication = regions_file$Number.CpGs
  )
} else if (source == "3000 genes") {
  regions_file <- read.csv(file = "/path/to/spermatogenesis-methylome/annotation/3000_genes_locations.csv")
  gene_names <- regions_file$gene_name
  promotor_names <- regions_file$promoter_name
  regions <- regions_file$coordinates..hg38.
  df <- data.frame(Gene_Names = gene_names,
                   Promotor_Names = promotor_names,
                   Region = regions
  )

} else {
  print("choose source/file with positions for methylation level")
  break
}

regions_split <- str_split(regions, ":|-")

for (file in files){
  print(file)
  methylation <- import.bedGraph(paste(file_path, file, sep = "/"))
  metylation_levels <- c()
  number_CPGs<- c()
  coverage <- c()
  for (region in regions_split) {
    chr <- strtoi(str_remove(region[1], "chr"))
    if(is.na(chr)) {
      chr <- str_remove(region[1], "chr")
    }
    start <- strtoi(region[2])
    end <- strtoi(region[3])
    region_object <- GRanges(seqnames = chr, ranges = IRanges(start = start+1, end = end+1))
    overlap <- subsetByOverlaps(methylation, region_object)
    overlap_mincov <- overlap[(overlap$NA. + overlap$NA..1) >=5]
    metylation_levels <- c(metylation_levels, mean(overlap_mincov$score))
    number_CPGs <- c(number_CPGs, length(overlap_mincov$score))
  }
  sample <- paste("number_CPGs", tail(str_split(file, "/")[[1]], n= 1), sep = "_")
  df <- cbind(df, number_CPGs)
  names(df)[length(names(df))] <- sample
  df <- cbind(df, metylation_levels/100)
  names(df)[length(names(df))] <- paste("Methylation_Level", tail(str_split(file, "/")[[1]], n= 1), sep = "_")
}

#write.csv(df, file = "/path/to/results/methylation_level_imprinted_genes.csv")
