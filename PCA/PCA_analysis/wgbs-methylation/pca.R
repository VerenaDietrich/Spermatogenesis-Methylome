library("data.table")
#library("bsseq")
library("ggplot2")
library("RColorBrewer")
library("scatterplot3d")
wgbs.methylation.convertGSEToCov <- function (coveragePath, sampleNames) {
  
  fileNames <- paste(coveragePath, "/", groups$sample, ".bed.gz", sep = "")
  
  for (sampleName in sampleNames) {
    
    bedFile <- paste(coveragePath, "/", sampleName, ".bed.gz", sep = "")
    covFile <- paste(coveragePath, "/", sampleName, ".cov.gz", sep = "")
    
    rawData <- fread(paste('zcat', bedFile), data.table = FALSE)
    
    convertedData <- data.frame(
      chr = substr(rawData[,1], 4, nchar(rawData[,1])),
      pos1 = rawData$Pos,
      pos2 = rawData$Pos,
      rate = rawData$MetRate,
      met = rawData$Met,
      unmet = rawData$UnMet
    )
    
    gzippedFile <- gzfile(covFile, "w")
    write.table(convertedData, gzippedFile, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  }
  
}

wgbs.methylation.loadAllData <- function (coveragePath, groups, cores = 1) {

  fileNames <- paste(coveragePath, "/", groups$sample, ".cov.gz", sep = "")

  coverageData <- read.bismark(fileNames, groups$sample, strandCollapse = FALSE, fileType = "cov", mc.cores = cores)

  smoothedData <- BSmooth(coverageData, mc.cores = cores, parallelBy = "chromosome")

  return(smoothedData)
}

wgbs.methylation.performPCA <- function (smoothedData, groups, minSamplesPerLocus = nrow(groups), excludeChromosomes = "X", cores = 4) {
  
  print("Performing pca with min samples per locus of:")
  print(minSamplesPerLocus)
  
  covg <- getCoverage(smoothedData)
  covg[covg > 0] <- 1
  
  groupedRowSums <- sapply(unique(groups$group), function (g) {

    samplesInGroup <- subset(groups, group == g)$sample
    groupCovg <- covg[,samplesInGroup]

    return(rowSums(groupCovg > 0))
  })

  loci <- granges(smoothedData)
  
  rowsSufficientlyCovered <- rowSums(covg > 0) >= minSamplesPerLocus &
    rowSums(groupedRowSums) >= length(unique(groups$group))
    
  validLoci <- !(seqnames(loci) %in% excludeChromosomes)
  
  keepingRows <- as.logical(rowsSufficientlyCovered & validLoci) # why as.logical? because of 'R'easons
  
  print(paste("Excluding loci from chromosomes:", excludeChromosomes))
  print(paste("Total loci excluded:",sum(!validLoci)))
  
  print("Rows with sufficient coverage:")
  print(sum(rowsSufficientlyCovered))
  
  print("Keeping in total:")
  print(sum(keepingRows))
  
  print("Sample coverage distribution (colSums):")
  print(colSums(covg[keepingRows,] > 0))
  
  plot(ggplot(data.frame(chr = seqnames(loci[keepingRows])), aes(x = chr)) + geom_histogram(stat = "count") + ggtitle("Loci included in PCA analysis"))
  
  methylData <- getMeth(smoothedData, type = "smooth")

  pcaSource <- methylData[keepingRows,]
  pcaSource[is.na(pcaSource)] <- 0.5
  
  pca <- prcomp(t(pcaSource))
  
  return(pca)
}

wgbs.methylation.plotPCA <- function (pca, groups, type = "TWO_D") {

  data <- data.frame(
    PC1 = pca$x[,1],
    PC2 = pca$x[,2],
    PC3 = pca$x[,3],
    group = groups$group,
    sample = groups$sample
  )
  
  switch(type,
         TWO_D = {
           ggplot(data, aes(x = PC1, y = PC2, colour = group, label = sample)) + geom_text()
         },
         THREE_D = {
           pal <- brewer.pal(length(unique(groups$group)), "Set1")
           colors <- pal[as.numeric(as.factor(data$group))]
           scatterplot3d(data$PC1, data$PC2, data$PC3, type = "h", color = colors, xlab = "PC1", ylab = "PC2", zlab = "PC3")
           legend("bottomright", legend = levels(as.factor(data$group)), col = pal, pch = 1)
  })
}