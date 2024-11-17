# compare WT organoids to organoid markers described here:
#   - 2015 Paşca et al. (https://doi.org/10.1038/nmeth.3415) sup table 1
#   - 2019 Sloan et al. (https://doi.org/10.1038%2Fs41596-018-0032-7) table 1
#   - 2023 Mulder et al. (https://doi.org/10.1186/s13287-023-03302-x) sup file 5

# ---- Step 0: Set up and package management ----
library(readxl)
# plots
library(ggplot2)
library(patchwork)
library(gridExtra)
library(grid)
library(pheatmap)
library(ggpubr)
# data manipluation
library(stringr)
library(tibble)
library(dplyr)
library(tidyverse)
library(data.table)

#library(dendextend)
resultDir <- "Z:/mnt/Data1/PROJECTS/ZMYND8/Organoid_BulkRNA-seq/Analysis_Results_Kaya/"
inputDir <- "Z:/mnt/Data1/PROJECTS/ZMYND8/Organoid_BulkRNA-seq/Data/"

# ---- 1. Take in data from organoid paper and from WT DE -----
organoidMarkers <- read_excel(paste0(inputDir, "organoidMarkers.xlsx"))
WTDERes <- read.csv("Z:/mnt/Data1/PROJECTS/ZMYND8/Organoid_BulkRNA-seq/Analysis_results_Kaya/DE/WT/mergedDEAcrossTimepointsWT.csv")
allWTRes <- readRDS("Z:/mnt/Data1/PROJECTS/ZMYND8/Organoid_BulkRNA-seq/Analysis_results_Kaya/DE/WT/WT_DEResults.rds")
preprocessedData <- readRDS("Z:/mnt/Data1/PROJECTS/ZMYND8/Organoid_BulkRNA-seq/Analysis_Results_Kaya/preprocessing/preprocessedData.rds")
IDToName <- preprocessedData$IDToName

# ---- 2. Select genes found in both lists and merge -----

# merged contains DE data and information about the marker gene:
# - only for genes where both are present
# - only when the gene is significantly DE in at least 1 timepoint comparison
filteredMarkers <- organoidMarkers %>% filter(Gene %in% WTDERes$geneNames)
matchRes <- WTDERes[match(filteredMarkers$Gene, WTDERes$geneNames), ]
merged <- cbind(filteredMarkers, matchRes) %>%
  dplyr::select(geneIDs, everything())

allGeneNames <- IDToName[row.names(allWTRes$D20vsD60)]

# check gene names match up after merging
stopifnot(sum(merged$Gene != merged$geneNames) == 0)
merged <- merged %>% dplyr::select(-geneNames)

# ---- 3. Save csv for combined organoid marker and DE raw data -----
fwrite(merged, paste0(resultDir, "DE/WT/organoidMarkersWithWTRes.csv"))

# ---- 4. plot TPM for the organoid marker genes ----
mk <- merged$Gene
meta <- preprocessedData$metadata %>% filter(Genotype == "WT")
tpm <- preprocessedData$filteredGeneTPM[c("gene_id", "gene_name", meta$Sample)]
markerTPMPlots <- list()
pdf(paste0(resultDir, "markerGenes/TPMMarkerPlotsWT.pdf"), height = 8, width = 8)
for (gene in mk) {
  # main plot based off "Z:/mnt/Data1/PROJECTS/ZMYND8/Organoid_BulkRNA-seq/Scripts_Kaya/Olivia_ThesisPlots_IV.R")
  plotdata <- data.frame(colnames(tpm[, -c(1, 2)]),
                         t(tpm[which(tpm$gene_name %in% gene), -c(1, 2)]))
  colnames(plotdata) <- c("Sample", "TPM")
  plotdata <- merge(plotdata, meta, by = "Sample")

  markerTPMPlots[[gene]] <- ggplot(plotdata, aes(y = TPM, x = factor(Day, levels = c("D20", "D60", "D100")))) +
    geom_boxplot(width = 0.3, size = 0.4, staplewidth = 0.3,
                 position = position_dodge(0.8)) +
  
    geom_dotplot(
      aes(fill = Batch, color = Batch),
      binaxis = "y", stackdir = "center", dotsize = 0.8,
      position = position_dodge(0.2)
      ) +
    labs(title = gene, x = "Day")
  print(markerTPMPlots[[gene]])
}
dev.off()


# a way to look at everything together inlcuding the data about the organoid markers
pdf(paste0(resultDir, "markerGenes/TPMMarkerPlotsWTCombinedInfo.pdf"),
    height = 5, width = 7)
# in this case mk is a list of marker genes that are also significantly DE
matchDEs <- list()
geneTables <- list()
for (gene in mk) {
  geneData <- filteredMarkers[filteredMarkers$Gene == gene, ]
  geneMain <- geneData[1]

  # creates data frame of all the markers with a positive value for the current gene
  marker <- geneData[-c(1, 2)][colSums(geneData[-c(1, 2)]) > 0]
  # add number of occurences of each cell type to marker data
  colnames(marker) <- str_c(colnames(marker), " (", as.character(marker[1, ]), ")")
  # Convert to a string of all the markers for the current gene
  marker <- str_flatten(colnames(marker), collapse = ",\n")
  # Include the overall marker category of the gene
  geneTable <- cbind(geneMain, marker)

  # Create table of the DE results for the current gene
  matchDE <- matchRes[matchRes$geneNames == gene, ][c(3:8)] %>%
    signif(digits = 3)
  rownames(matchDE) <- "DE Results"
  colnames(matchDE) <- c("D20vsD60: log2FC", "padj",
                         "D20vsD100: log2FC", "padj",
                         "D60vsD100: log2FC", "padj")
  
  # set up where each section is plotted
  # - marker WT plot across the left side and the two tables on the right
  layout <- matrix(c(1, 2,
                     1, 3), ncol = 2, byrow = TRUE)
  t1 <- tableGrob(t(geneTable), theme = ttheme_default(base_size = 10))
  t2 <- tableGrob(t(matchDE), theme = ttheme_default(base_size = 10))

  print(grid.arrange(markerTPMPlots[[gene]], t2, t1,
                     layout_matrix = layout))

  geneTables[[gene]] <- geneTable
  matchDEs[[gene]] <- matchDE
}
dev.off()

# create dataframe ready to make tables for all DE genes, not just markers
matchDEsFull <-  WTDERes[c(3:8)] %>% signif(digits = 3)

# name unique genes by their name and duplicated ones by their ensemblID
DEGeneNames <- WTDERes$geneNames
isDuplicatedDE <- DEGeneNames %in% DEGeneNames[duplicated(DEGeneNames)]
DEGeneNames[isDuplicatedDE] <- WTDERes[isDuplicatedDE, ]$geneIDs
rownames(matchDEsFull) <- DEGeneNames
colnames(matchDEsFull) <- c("D20vsD60: log2FC", "padj",
                            "D20vsD100: log2FC", "padj",
                            "D60vsD100: log2FC", "padj")
matchDEs[[gene]] <- matchDE

# ---- save exact data needed for shiny app ----
# the marker list
# the merged data (should have the marker list)
# start with just the TPM data
shinyData <- list(metaData = preprocessedData$metadata,
                  TPM = preprocessedData$filteredGeneTPM)
# improve gene names because gene names directly from data has repeated values (e.g. U3)
isDuplicated <- shinyData$TPM$gene_name %in%
  shinyData$TPM$gene_name[duplicated(shinyData$TPM$gene_name)]
shinyData$TPM$gene_name[isDuplicated] <- shinyData$TPM[isDuplicated, c(1, 2)] %>%
  mutate(gene_name = paste0(gene_id, rep("(", 373), gene_name, rep(")", 373))) %>%
  .$gene_name
shinyData$markerTables <- geneTables
shinyData$sigDEs <- matchDEsFull

# An error occurs in the shiny app where newlines turn up as /n
# These are replaced with a newline character here
for (gene in mk) {
  shinyData$markerTables[[gene]]$marker <- str_replace_all(shinyData$markerTables[[gene]]$marker, "/n", "\n")
}

saveRDS(shinyData, "Z:/mnt/Data1/PROJECTS/ZMYND8/Organoid_BulkRNA-seq/Scripts_Kaya/shiny/WTMarkerGenes/shinyData.rds")


# ---- compare to pasca data ----
# Further conversions - not applied to data:
  # C14orf37 is ARMH4 (https://www.genecards.org/cgi-bin/carddisp.pl?gene=ARMH4)
  # > "ARMH4" %in% WTDERes$geneNames - FALSE
  # NAT6 is NAA80 (https://www.ncbi.nlm.nih.gov/gene/24142)
  # > "NAA80" %in% WTDERes$geneNames - FALSE
  # C1orf9 is https://www.genecards.org/cgi-bin/carddisp.pl?gene=SUCO
  # > "SUCO" %in% WTDERes$geneNames - TRUE but FC < 1.5
  # TROVE2 is TROVE2
  # C1orf124 is "DVC1" https://www.nature.com/articles/nsmb.2395
  # incomplete

validateDE <- read_excel("Z:/mnt/Data1/PROJECTS/ZMYND8/Organoid_BulkRNA-seq/Data/rawDETableFromSergiu2015.xlsx")
validateDE <- validateDE %>% rename_at("LogFC_3D    differentiation", ~"LogFC_3D") %>%
  rename_at("Pval_3D    Differentiation", ~"Pval_3D") %>%
  rename_at("PVal_Fetal    brain    stage    1vs6", ~"PVal_Fetal_brain_stage_1vs6") %>%
  rename_at("LogFC_Fetal    Brain    stage    1vs6", ~"LogFC_Fetal_brain_stage_1vs6")
validateDE$geneSymbol[validateDE$geneSymbol %in% allGeneNames]
valDEMatch <- validateDE[validateDE$geneSymbol %in% WTDERes$geneNames, ]
validateDE$geneSymbol[validateDE$geneSymbol %in% matchRes$geneNames]

valDEWT <- cbind(WTDERes[match(valDEMatch$geneSymbol, WTDERes$geneNames), ],
                 valDEMatch)

validationMerged <- valDEWT %>%
  dplyr::select(geneIDs, geneNames,
                D20vsD100.log2FC, D20vsD100.padj,
                LogFC_3D, Pval_3D,
                LogFC_Monolayer, PVal_Monolayer,
                LogFC_Fetal_brain_stage_1vs6, PVal_Fetal_brain_stage_1vs6)


fwrite(validationMerged, "Z:/mnt/Data1/PROJECTS/ZMYND8/Organoid_BulkRNA-seq/Analysis_results_Kaya/DE/WT/validationD20vsD100.csv")


# ---- Perform hypergeometric test to see if the overlap is significant ----
# number of overlapping genes = 23
overlap <- sum(validationMerged$D20vsD100.padj < 0.05 &
                 abs(validationMerged$D20vsD100.log2FC) > log2(1.5))
# number of our DE genes, D20 vs D100 = 4769
DE20vs100 <- allWTRes$D20vsD100[!is.na(allWTRes$D20vsD100$padj) &
                                  allWTRes$D20vsD100$padj < 0.05 &
                                  abs(allWTRes$D20vsD100$log2FoldChange) > log2(1.5), ]
DE_Ourdata <- nrow(DE20vs100)
# all genes in our data = 23116
AllGenes_Ourdata <- nrow(preprocessedData$filteredCounts)
# num DE genes in their data = 71
DE_TheirData <- nrow(validateDE)
hypergeomTest <- phyper(overlap, DE_Ourdata, AllGenes_Ourdata - DE_Ourdata,
                        DE_TheirData, lower.tail = FALSE)

# interesting, our organoids are more significantly correlated to fetal data than pasca's
lmVsOrganoid <- lm(D20vsD100.log2FC ~ LogFC_3D, data = validationMerged)
lmVsFetal <- lm(D20vsD100.log2FC ~ LogFC_Fetal_brain_stage_1vs6, data = validationMerged)
lmPascaOnly <- lm(LogFC_3D ~ LogFC_Fetal_brain_stage_1vs6, data = validationMerged)

# lm uses pearsons by default too, so this is identical to the value used for graphing
plotFCvsFC <- function(dataset, var1Name, var2Name,
                       upLabelCoord, downLabelCoord, pearsonLabelCoord,
                       labNames, relSize = 3) {
  pearsonRes <- cor.test(dataset[[var1Name]], dataset[[var2Name]], method = "pearson")
  pearsonText <- paste("Pearson's correlation coefficient:",
                       paste("p-val =", signif(pearsonRes$p.value, 3)),
                       paste("cor =", signif(pearsonRes$estimate[[1]], 3)),
                       sep = "\n")
  plot <- ggplot(dataset, aes(x = .data[[var1Name]], y = .data[[var2Name]])) +
    geom_point() +
    geom_smooth(method = lm, color = "black") +
    annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf,
             alpha = 0.2, fill = "red") +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
             alpha = 0.2, fill = "blue") +
    annotate("text", x = upLabelCoord$x, y = upLabelCoord$y,
             label = "Upregulated in both", color = "red", size = rel(relSize)) +
    annotate("text", x = downLabelCoord$x, y = downLabelCoord$y,
             label = "Downregulated in both", color = "blue", size = rel(relSize)) +
    annotate("label", x = pearsonLabelCoord$x,  y = pearsonLabelCoord$y,
             label = pearsonText, size = rel(relSize)) +
    labs(title = labNames$title,
         x = labNames$x, y = labNames$y)
  return(plot)
}

# visually show the correlation evidenced by the 0.006669962 hypergeometric test result
pdf(file = "Z:/mnt/Data1/PROJECTS/ZMYND8/Organoid_BulkRNA-seq/Analysis_Results_Kaya/markerGenes/scatterPlotsOurOrganoidsVsSloan.pdf")
# our organoids vs Pasca organoids:
upLabelCoord <- data.frame(x = 1, y = 2.5)
downLabelCoord <- data.frame(x = -2.5, y = -0.25)
pearsonLabelCoord <- data.frame(x = 2, y = -2)
labNamesOurVsPasca <- data.frame(title = "FC in our organoids D20 to D100 vs\noriginal Sloan organoids D52 to D76",
                                 x = "Log2FC D20 vs D100 Voineagu organoids",
                                 y = "LogFC D52 vs D76 Pasca organoids")
print(plotFCvsFC(validationMerged, "D20vsD100.log2FC", "LogFC_3D",
                 upLabelCoord, downLabelCoord, pearsonLabelCoord,
                 labNamesOurVsPasca))

# our organoids vs Pasca fetal stage 1-6:
upLabelCoord <- data.frame(x = 1, y = 4)
downLabelCoord <- data.frame(x = -2.5, y = -0.3)
pearsonLabelCoord <- data.frame(x = 2, y = -3.7)
labNamesOurVsFetal <- data.frame(title = "FC in our organoids D20 to D100\nvs fetal stage 1 to 6",
                                 x = "Log2FC D20 vs D100 organoids",
                                 y = "LogFC fetal stage 1 to 6")
print(plotFCvsFC(validationMerged, "D20vsD100.log2FC", "LogFC_Fetal_brain_stage_1vs6",
                 upLabelCoord, downLabelCoord, pearsonLabelCoord,
                 labNamesOurVsFetal))

# Pasca organoids vs Pasca fetal stage 1-6:
upLabelCoord <- data.frame(x = 0.5, y = 4)
downLabelCoord <- data.frame(x = -0.9, y = -0.3)
pearsonLabelCoord <- data.frame(x = 1.3,  y = -2.5)
labNamesPascaVsFetal <- data.frame(title = "FC in original Pasca organoids D52 to D76\nvs fetal stage 1 to 6",
                                   x = "LogFC D52 vs D76 Pasca organoids",
                                   y = "LogFC fetal stage 1 to 6")
print(plotFCvsFC(validationMerged, "LogFC_3D", "LogFC_Fetal_brain_stage_1vs6",
                 upLabelCoord, downLabelCoord, pearsonLabelCoord,
                 labNamesPascaVsFetal))
dev.off()

# ---- heatmap of TPM of main markers ----

# heamap row scaling function from https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/
cal_z_score <- function(x) {
  (x - mean(x)) / sd(x)
}

# load list of markers with cell type info and expression start time data
# From 2019 Sloan et al. https://doi.org/10.1038%2Fs41596-018-0032-7 TABLE 1
SloanMarkers <- read_excel("Z:/mnt/Data1/PROJECTS/ZMYND8/Organoid_BulkRNA-seq/Data/SloanMarkers.xlsx")
SloanMarkers <- SloanMarkers %>% mutate(expStart = str_replace(`Temporal expression`, "From ~ ", ""),
                                        markerOf = `Spatial expression and specificity`) %>%
  select(-`Spatial expression and specificity`)
# to categorise the markers, the start time should be viewed as numeric
SloanMarkers$numericExpStart <- SloanMarkers$expStart %>% gsub("day ", "", .) %>%
  gsub("Throughout differentiation", "0", .) %>%
  as.numeric()

WTMarkerTPMRaw <- tpm[tpm$gene_name %in% SloanMarkers$Marker, ]
# sloan markers not in TPM:
sum(!SloanMarkers$Marker %in% WTMarkerTPMRaw$gene_name) == 0
# none - excellent

# prep TPM for plotting + normalise
WTMarkerTPM <- WTMarkerTPMRaw %>% remove_rownames() %>%
  column_to_rownames("gene_name") %>% select(c(2:10))
colnames(WTMarkerTPM) <- str_replace(colnames(WTMarkerTPM), "_GOK.*", "")
WTMarkerTPM <- WTMarkerTPM %>%
  select("B2_D20_WT", "B7_D20_WT", "B4_D21_WT", "B3_D22_WT",
         "B2_D60_WT", "B3_D60_WT", "B4_D60_WT",
         "B4_D100_WT", "B5_D100_WT") %>%
  as.matrix()
WTMarkerTPMNorm <- t(apply(WTMarkerTPM, 1, cal_z_score))
# unannotated heatmap
pheatmap(WTMarkerTPMNorm)

# Group into progenitors and other markers:
rowAnn <- data.frame(maturity = SloanMarkers$markerOf)
maturityCategories <- c("Progenitors", "Deep layer neurons",
                        "Superficial layer neurons", "GABAergic neurons",
                        "Interneurons", "Astrocytes")
for (maturityTerm in maturityCategories) {
  rowAnn$maturity[grepl(maturityTerm, SloanMarkers$markerOf, ignore.case = TRUE)] <- maturityTerm
}
rownames(rowAnn) <- SloanMarkers$Marker

# Group into time expression begins:
rowAnn$startTime <- SloanMarkers$expStart
rowAnn$startTime[SloanMarkers$numericExpStart >= 20 &
                   SloanMarkers$numericExpStart <= 30] <- "day 20-30"
rowAnn$startTime["Throughout differentiation" == SloanMarkers$expStart] <- "all"
rowAnn$startTime <- factor(rowAnn$startTime, levels = c("all", "day 20-30", "day 50", "day 75", "day 100", "day 200"))

# simpleRA - simpler/clearer row annotations
# Group into progenitors and mature markers:
simpleRA <- data.frame(maturity = rep("Mature", 16))
simpleRA$maturity[grepl("progenitor", SloanMarkers$markerOf)] <- "Progenitor"
rownames(simpleRA) <- SloanMarkers$Marker

simpleRA$startTime <- SloanMarkers$expStart
simpleRA$startTime[SloanMarkers$numericExpStart <= 20] <- "by day 20"
simpleRA$startTime[SloanMarkers$numericExpStart > 20 &
                     SloanMarkers$numericExpStart <= 50] <- "by day 50"
simpleRA$startTime[SloanMarkers$numericExpStart > 50 &
                     SloanMarkers$numericExpStart <= 100] <- "by day 100"
simpleRA$startTime["day 200" == SloanMarkers$expStart] <- "by day 200"
simpleRA$startTime <- factor(simpleRA$startTime, levels = c("by day 20", "by day 50",
                                                            "by day 100", "by day 200"))

# column annotations of day and batch
meta$SampleShort <- str_replace(meta$Sample, "_GOK.*", "")
colAnn <- data.frame(day = meta$Day[match(colnames(WTMarkerTPM), meta$SampleShort)])
colAnn$batch <- meta$Batch[match(colnames(WTMarkerTPM), meta$SampleShort)]
colAnn$day <-  factor(colAnn$day, levels = c("D20", "D60", "D100"))

rownames(colAnn) <- colnames(WTMarkerTPM)

# Using the callback function from the pheatmap help page and the discussion here https://stackoverflow.com/questions/73872992/pheatmap-manually-re-order-leaves-in-dendogram:
# Reorder cluster while keeping structure intact
callback <- function(hc, mat) {
  sv <- svd(t(mat))$v[7, ]
  dend <- reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

pdf(file = "Z:/mnt/Data1/PROJECTS/ZMYND8/Organoid_BulkRNA-seq/Analysis_Results_Kaya/markerGenes/heatmapsWTTPM.pdf")
print(pheatmap(WTMarkerTPMNorm, annotation_row = simpleRA, annotation_col = colAnn, clustering_callback = callback))
dev.off()

# ---- create heatmap of larger group of marker genes ----
orgMarkerTPMRaw <- tpm[tpm$gene_name %in% organoidMarkers$Gene, ]

# prep TPM for plotting + normalise
orgMarkerTPM <- orgMarkerTPMRaw %>% remove_rownames() %>%
  column_to_rownames("gene_name") %>% select(c(2:10))
colnames(orgMarkerTPM) <- str_replace(colnames(orgMarkerTPM), "_GOK.*", "")
orgMarkerTPM <- orgMarkerTPM %>%
  select("B2_D20_WT", "B7_D20_WT", "B4_D21_WT", "B3_D22_WT",
         "B2_D60_WT", "B3_D60_WT", "B4_D60_WT",
         "B4_D100_WT", "B5_D100_WT") %>%
  as.matrix()
orgMarkerTPMNorm <- t(apply(orgMarkerTPM, 1, cal_z_score))

# unnanotated heatmap
pheatmap(orgMarkerTPMNorm, fontsize = 10, fontsize_row = 6)
# GABA is removed as it is not a gene
organoidMarkers <- organoidMarkers[organoidMarkers$Gene != "GABA", ]

# sort markers into progenitor or not based on the highest value in the organoid marker data
progMarkers <- colnames(organoidMarkers)[str_detect(colnames(organoidMarkers), "PRG")]
organoidMarkersHM <- organoidMarkers %>% column_to_rownames("Gene") %>% select(-c("Category"))
organoidMarkersHM$mainMarker <- colnames(organoidMarkersHM)[apply(organoidMarkersHM, 1, which.max)]
organoidMarkersHM$progenitor <- organoidMarkersHM$mainMarker %in% progMarkers
progGenes <- rownames(organoidMarkersHM[organoidMarkersHM$progenitor, ])
matureGenes <- rownames(organoidMarkersHM[!organoidMarkersHM$progenitor, ])


# put together row annotation
bigRowAnn <- data.frame(maturity = rep("Mature", 55))
bigRowAnn$maturity[organoidMarkersHM$progenitor] <- "Progenitor"     
rownames(bigRowAnn) <- rownames(organoidMarkersHM)

# annotated heatmap
pheatmap(orgMarkerTPMNorm, annotation_row = bigRowAnn, fontsize = 10, fontsize_row = 6)

# heatmap without dendogram but split into days and marker categories
orgMarkerTPMNormSort <- rbind(orgMarkerTPMNorm[progGenes, ], orgMarkerTPMNorm[matureGenes, ])
pdf(file = "Z:/mnt/Data1/PROJECTS/ZMYND8/Organoid_BulkRNA-seq/Analysis_Results_Kaya/markerGenes/heatmaps55MarkerWTTPM.pdf")
extendedMarkerHM <- pheatmap(orgMarkerTPMNormSort, annotation_row = bigRowAnn,
                             annotation_col = colAnn, fontsize = 10,
                             fontsize_row = 5, cluster_rows = FALSE,
                             cluster_cols = FALSE, gaps_row = 24, gaps_col = c(4, 7),
                             labels_col = str_replace(colnames(orgMarkerTPMNormSort), "_WT", ""),
                             main = "Expression of marker genes across WT organoids")
print(extendedMarkerHM)
dev.off()

# ---- set up plots for practice thesis ----
# plan - large heatmap at top, then three interesting genes below
# BCL11B: drop from D60 to D100 confirmed by literature. Though they also described high expression early on
# RELN: increases over time clearly from D20 - 60. Found on organoid surface at D62 (sergiu)
# NES: codes for the common progenitor marker, nestin
# PAX6: strange to see it absent early on, as it is a progenitor marker, and has been observed at D18. Though the later expression could be due to layers
# OLIG2: oligodendrocyte progenitor. The pattern implies development of oligo progenitors over time, and then perhaps differentiation?
markerTPMPlots$BCL11B
markerTPMPlots$PAX6
markerTPMPlots$OLIG2
pdf(file = "Z:/mnt/Data1/PROJECTS/ZMYND8/Organoid_BulkRNA-seq/Analysis_Results_Kaya/markerGenes/heatmapAndSelectedMarkers.pdf") 
row1 <- ggarrange(extendedMarkerHM[[4]], labels = "A")
row2 <- ggarrange(markerTPMPlots$NES,
                  markerTPMPlots$PAX6,
                  markerTPMPlots$OLIG2, labels = c("B", "C", "D"),
                  ncol = 3, nrow = 1, common.legend = TRUE, legend = "right")
ggarrange(row1, row2, ncol = 1, heights = c(2, 1))
dev.off()


# next, plot hypergeometric distribution
hypergeomTest <- phyper(overlap, DE_Ourdata, AllGenes_Ourdata - DE_Ourdata, DE_TheirData, lower.tail = FALSE)
hyperDist <- data.frame(overlapNum = rhyper(100000, DE_Ourdata, AllGenes_Ourdata - DE_Ourdata, DE_TheirData))
distPlot <- ggplot(hyperDist, aes(x = overlapNum)) +
  geom_histogram(aes(y = after_stat(count) / sum(after_stat(count))),
                 binwidth = 1) 

pdf("Z:/mnt/Data1/PROJECTS/ZMYND8/Organoid_BulkRNA-seq/Analysis_Results_Kaya/markerGenes/hypergeometricDist.pdf")
overlapDiagram <- distPlot +
  geom_vline(xintercept = 27.5, color = "red", linetype = 'dashed') +
  annotate("rect", xmin = 27.5, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.2) +
  annotate("label", x = 6, y = 0.11, label = paste0("Probability of \n≥ 23 gene overlap:\n P = ",
                                                    formatC(hypergeomTest, format = "e", digits = 2)),
           size = rel(2.5)) +
  annotate("text", x = 30.75, y = 0.1125, label = "Overlap ≥ 23", color = "red", size = rel(2)) +
  labs(title = "Hypergeometric test of overlap\nbetween known organoid DE genes\nand those identified in our organoids",
       x = "Number of overlapping DE genes",
       y = "Frequency of overlap")
print(overlapDiagram)
dev.off()

upLabelCoord <- data.frame(x = 1.1, y = 2.5)
downLabelCoord <- data.frame(x = -2.3, y = -0.25)
pearsonLabelCoord <- data.frame(x = 2,  y = -2)
plotOrgVsOrg <- plotFCvsFC(validationMerged, "D20vsD100.log2FC", "LogFC_3D",
                           upLabelCoord, downLabelCoord, pearsonLabelCoord,
                           labNamesOurVsPasca, 2.5)

upLabelCoord <- data.frame(x = 1.1, y = 4)
downLabelCoord <- data.frame(x = -2.3, y = -0.3)
pearsonLabelCoord <- data.frame(x = 2,  y = -3.7)
plotOrgVsFetal <- plotFCvsFC(validationMerged, "D20vsD100.log2FC", "LogFC_Fetal_brain_stage_1vs6",
                             upLabelCoord, downLabelCoord, pearsonLabelCoord,
                             labNamesOurVsFetal, 2.5)

upLabelCoord <- data.frame(x = 0.6, y = 4)
downLabelCoord <- data.frame(x = -0.9, y = -0.3)
pearsonLabelCoord <- data.frame(x = 1.3,  y = -2.5)
plotLitOrgVsFetal <- plotFCvsFC(validationMerged, "LogFC_3D", "LogFC_Fetal_brain_stage_1vs6",
                                upLabelCoord, downLabelCoord, pearsonLabelCoord,
                                labNamesPascaVsFetal, 2.5)


pdf("Z:/mnt/Data1/PROJECTS/ZMYND8/Organoid_BulkRNA-seq/Analysis_Results_Kaya/markerGenes/hypergeomAndOrganoidComparison.pdf")
ggarrange(overlapDiagram, plotOrgVsOrg,
          plotOrgVsFetal, plotLitOrgVsFetal,
          ncol = 2, nrow = 2, labels = "AUTO")
dev.off()


# Looking at D20 ventral marker variation

WTDERes <- read.csv("Z:/mnt/Data1/PROJECTS/ZMYND8/Organoid_BulkRNA-seq/Analysis_results_Kaya/DE/WT/mergedDEAcrossTimepointsWT.csv")
allWTRes <- readRDS("Z:/mnt/Data1/PROJECTS/ZMYND8/Organoid_BulkRNA-seq/Analysis_results_Kaya/DE/WT/WT_DEResults.rds")
preprocessedData <- readRDS("Z:/mnt/Data1/PROJECTS/ZMYND8/Organoid_BulkRNA-seq/Analysis_Results_Kaya/preprocessing/preprocessedData.rds")
IDToName <- preprocessedData$IDToName
meta <- preprocessedData$metadata

#ventral MGE-like progenitor markers: NKX2.1 (or NKX2-1) and LHX6.
MGEs <- c("NKX2-1", "LHX6")
#ventral LGE-like progenitor markers: DLX5, GSX2, BCL11B, MEIS2.
LGEs <- c("DLX5", "GSX2", "BCL11B", "MEIS2")
#dorsal progenitor markers: PAX6.
dor <- c("PAX6")
#A ventral LGE-like marker (for hStrS) that seems to be high in our D20 WT organoid in our scRNA-seq data: SP8.
matLGEs <- c("SP8")

pivot_longer(data.frame(MGEs, LGEs, dor, matLGEs), 1)

TPM <- preprocessedData$filteredGeneTPM
D20TPM <- TPM[, c("gene_id", "gene_name", meta$Sample[meta$Day == "D20" & meta$Genotype != "Het"])]
D20TPM$gene_name <- IDToName[D20TPM$gene_id]
rownames(D20TPM) <- D20TPM$gene_id

EMarkers <- list("LGEs" = LGEs, "MGEs" = MGEs, "dor" = dor, "matLGEs" = matLGEs)
allMarkers <- unlist(EMarkers)

geneToMarkerOf <- function(gene, markerList) {
  for (markerType in names(markerList)) {
    if (gene %in% markerList[[markerType]]) {
      return(markerType)
    }
  }
}

markerTPMRaw <- D20TPM[D20TPM$gene_name %in% allMarkers, ]
# markers not in TPM:
sum(!allMarkers %in% markerTPMRaw$gene_name) == 0
# none - excellent

markerTPMRaw$markerOf <- unlist(lapply(markerTPMRaw$gene_name, geneToMarkerOf, markerList = EMarkers))

longTPM <- pivot_longer(markerTPMRaw, c(4:ncol(markerTPMRaw) - 1))
longTPM$Batch <- as.factor(gsub("_D.*", "", longTPM$name))
longTPM$Genotype <- as.factor(gsub("_.*", "", gsub(".*D[0-9]*_", "", longTPM$name)))
# WT only - D20 expression between batches
WTPlot <- ggplot(longTPM[longTPM$Genotype == "WT",], aes(x = gene_name, y = value)) +
  geom_boxplot() + 
  geom_point(aes(color = Batch)) + 
  labs(title = "WT D20 marker expression",y = "TPM", x = "Marker genes") +
  facet_grid(~markerOf, scales='free', space="free")

# comparing D20 WT vs KO
WTKOPlot <- ggplot(longTPM, aes(x = gene_name, y = value)) +
  geom_boxplot(aes(fill = Genotype)) + 
  geom_point(aes(color = Batch, group = Genotype), 
             position=position_dodge(0.75)) + 
  labs(title = "WT and KO D20 marker expression", y = "TPM", x = "Marker genes") +
  facet_grid(~markerOf, scales='free', space="free")

# exporting graphs
pdf(paste0(resultDir, "markerGenes/ganglionicEminenceD20.pdf"))
print(WTPlot)
print(WTKOPlot)
ggarrange(WTPlot,WTKOPlot, nrow = 2)
dev.off()

# are any of these genes in the DE results?
D20DE <- read.csv("Z:/mnt/Data1/PROJECTS/ZMYND8/Organoid_BulkRNA-seq/Analysis_results_Kaya/DE/timepoints/DEGenesD20WTvsKO.csv")
# yes, all markers are significantly DE
allMarkers %in% D20DE$geneName
D20DE[D20DE$geneName %in% allMarkers,]