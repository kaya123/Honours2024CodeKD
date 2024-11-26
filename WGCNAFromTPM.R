# Overview from https://docs.google.com/document/d/1AJ1IW-3AZ9_yB2gNR8Ri7irl-GU-JHTVbN1nZKAPoV8/edit
# WGCNA - Use the Salmon TPM data, with log2 transformation


# Proposed program structure:
# 1. Take in preprocessed data from "Z:/mnt/Data1/PROJECTS/ZMYND8/Organoid_BulkRNA-seq/Analysis_Results_Kaya/preprocessing/preprocessedData.rds")
# 2. Perform WCGNA on log2 TPM
# 3. Get Module Eigengenes per cluster 

# ---- Step 0: Set up and package management ---
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("Rgraphviz", suppressUpdates = TRUE) 

#install.packages("dplyr")
#install.packages("WGCNA") # networking
#install.packages("tidyverse") # data manipulation

library(WGCNA)
library(openxlsx)

# data manipulation
library(data.table)
library(dplyr)
library(tidyr)
library(genefilter)
library(tibble)
library(stringr)
# for GO
library(clusterProfiler)
library(enrichplot)
library(biomaRt)
library(Rgraphviz)
# plotting
library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)
library(pheatmap)
library(colorspace)
# fct_inorder function
library(forcats)

resultDir <- "Z:/mnt/Data1/PROJECTS/ZMYND8/Organoid_BulkRNA-seq/Analysis_Results_Kaya/"
netDir <- paste0(resultDir, "WGCNA/TPMThresh1IV/")
if (!dir.exists(netDir)) { dir.create(netDir) }
infoCols <- c("gene_id", "gene_name")
thresh <- 1

# ---- 1. Take in preprocessed data -----
preprocessedData <- readRDS(paste0(resultDir, "preprocessing/preprocessedData.rds"))
IDToName <- preprocessedData$IDToName
# WGCNA will be performed on only WT and KO data
meta <- preprocessedData$metadata
metaNet <- data.frame(meta) %>% 
  dplyr::select(-Clone) %>%
  filter(Genotype != "Het") %>%
  mutate(numericDay = as.numeric(Day %>% gsub("D","", .)))
rownames(metaNet) <- metaNet$Sample

TPM <- preprocessedData$geneTPM
TPM <- TPM[, -grep("Het", colnames(TPM))]
TPM_filtered<- TPM[rowSums(TPM[, -c(1,2)] > thresh ) >= 2 , ]

log2TPM <- TPM_filtered 
log2TPM[,-c(1,2)] <- log2(TPM_filtered[, -c(1,2)] + 1 )


# ---- 2. Perform WCGNA on log2 TPM ----
# remove row names row and transpose ready for WGCNA processing
rownames(log2TPM) <- log2TPM$gene_id
WGCNA_TPM <- log2TPM %>% dplyr::select(-c(gene_id, gene_name)) %>% t()

# test range of powers to pick lowest soft threshold where R^2 > 0.8
allowWGCNAThreads() 
powers = c(c(1:10), seq(from = 12, to = 30, by = 2))
powers = seq(from = 2, to = 20, by = 2)
sft = pickSoftThreshold(
  WGCNA_TPM,
  powerVector = powers,
  verbose = 7
)

# produce and save scale free topology plot
pdf(file = paste0(netDir, "WTKOScaleFreeTopologyTPM2.pdf"), width = 8, height = 5)
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type="n",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.8, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")
dev.off()

#https://bioinformaticsworkbook.org/tutorials/wgcna.html#gsc.tab=0
# single block WGCNA networking
picked_power = 8
temp_cor <- cor
cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)
netwk <- blockwiseModules(WGCNA_TPM,
                          power = picked_power,                
                          networkType = "signed", # negative and positive correlations treated separately
                          
                          deepSplit = 2, # default
                          pamRespectsDendro = F, # tree formation
                          minModuleSize = 100,
                          maxBlockSize = 30000,
                          
                          reassignThreshold = 0, # skipped for manual reassignment
                          mergeCutHeight = 0.15, # default
                          
                          # Can save tom for faster run in future - but runs slower first time 
                          saveTOMs = FALSE,
                          
                          # output details
                          numericLabels = T,
                          verbose = 15,
                          nThreads = 24)
cor <- temp_cor     # Return cor function to original namespace

# save WGCNA results, if already run, can retrieve saved results
saveRDS(netwk, file = paste0(netDir, "network.rds"))
netwk <- readRDS(file = paste0(netDir, "network.rds"))

# ---- 3. Get Module Eigengenes per cluster ----
# set up module information table
modules <- as.data.frame(table(netwk$colors))
colnames(modules) <- c("Label", "nGenes")
modules <- modules %>% 
  mutate(Color = c("grey", labels2colors(Label[-1])),
         Label = paste0("ME", Label))
fwrite(modules, paste0(netDir, "modulesPreReassignment.csv"))

#---- Save Module Eigenegens
me <- data.frame(rownames(WGCNA_TPM), netwk$MEs)
colnames(me)[1] <- "Sample"
fwrite(me, paste0(netDir, "ME.csv"))

#---- Save KMEs
moduleLabel <- paste0("ME", netwk$colors)
moduleColor <- modules$Color[match(moduleLabel, modules$Label)]
modules$colLabel <- paste0("ME", modules$Color)

KMEs <- signedKME(WGCNA_TPM, netwk$MEs, outputColumnName = "M")
kme <- data.frame(log2TPM[,infoCols], moduleColor, moduleLabel, KMEs)
fwrite(kme, paste0(netDir, "kMEPreReassignment.csv"))

# ---- save dendorgram prereassignment ----
# Convert labels to colors for plotting
# dendrogram of modules and their connections to one another
allColors <- modules$Color[c(2:nrow(modules), 1)]
mergedColors <- labels2colors(netwk$colors, colorSeq = allColors)
                  
pdf(file =  paste0(netDir, "dendrogramAllModulesPreReassignment.pdf"), height = 5, width = 8)
par(mfrow = c(1,1))
plotDendroAndColors(netwk$dendrograms[[1]],
                    mergedColors,
                    "Module colors",
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05,
                    main = "Cluster dendrogram of WGCNA modules")

dev.off()


# ---- reassignment of modules based on kMEs ---- 
kme <- read.csv(paste0(netDir, "kMEPreReassignment.csv"))
rownames(kme) <- kme$gene_id
kme$moduleLabel <- gsub("ME", "M", kme$moduleLabel)
# Re-assign genes to modules by kME. Any not passing the criteria of correlation pvalue <0.05, and kME>0.1, will be assigned to the grey/M0 module
kmeInfoCols <- c(1:4)
kmedata <- kme[,-kmeInfoCols]
pvalBH <- kmedata
pvalBH[,] <- NA

for (i in c(1:ncol(pvalBH)))
{
  p <- corPvalueStudent(kmedata[,i], nSamples=18)
  pvalBH[,i] <- p.adjust(p, method = "BH")
}
kme$newModule <- NA
minChangedKME <- 1
for (j in c(1:nrow(kmedata)))
{
  if (j==1) print("Working on genes 1:10000")
  if (j==10000) print(paste("Working on genes 10000:", nrow(kmedata)))
  m <- which(kmedata[j,] == max(kmedata[j,]))
  if ((pvalBH[j,m]<0.05) & (kmedata[j,m]>0.1)) {
    if (kmedata[j,m]) {
      minChangedKME <- min(minChangedKME, kmedata[j,m])
      kme$newModule[j] <- as.character(colnames(kmedata)[m])
    }
  }
}
# Assign genes not associated to any module to M0
kme$newModule[is.na(kme$newModule)] <-"M0"
t <- as.data.frame(table(old = kme$moduleLabel,new = kme$newModule))
fwrite(t, paste0(netDir, "kMEReassignmentCounts.csv"))

leftGreyMod <- kme[kme$newModule != "M0" & kme$moduleLabel == "M0",]
fwrite(leftGreyMod, paste0(netDir, "reassignedFromGreyKMEs.csv"))


# Replacing the old moduleLabel and moduleColor columns with the updated newModule and newColor columns respectively
kme$newColor=kme$moduleColor[match(kme$newModule, kme$moduleLabel)]
kme$moduleLabel=kme$newModule; kme$moduleColor=kme$newColor
fwrite(kme, paste0(netDir, "kMEReassignmentTable.csv"))

kme <- kme[,-grep("newModule", colnames(kme))];
kme <- kme[,-grep("newColor", colnames(kme))]
fwrite(kme, paste0(netDir, "kME.csv"))

# ---- update module data
modules <- read.csv(paste0(netDir, "modulesPreReassignment.csv"))
for (color in modules$Color) {
  modules[modules$Color == color,"nGenes"] <- sum(kme$moduleColor == color)
}
fwrite(modules, paste0(netDir, "modules.csv"))


# ---- save dendorgram post-reassignment ----
# Convert labels to colors for plotting
# dendrogram of modules and their connections to one another
allColors <- modules$Color[c(2:nrow(modules), 1)]
mergedColors <- labels2colors(netwk$colors, colorSeq = allColors)
reasigned <- read.csv(paste0(netDir, "kMEReassignmentTable.csv"))
newColors <- reasigned$newColor[match(reasigned$gene_id, names(netwk$colors))]
colorAnnotation <- data.frame(old = mergedColors,
                              new = newColors)

pdf(file =  paste0(netDir, "dendrogramAllModules.pdf"), height = 5, width = 9)
par(mfrow = c(1,1))
plotDendroAndColors(netwk$dendrograms[[1]],
                    colorAnnotation,
                    c("pre-reassignment", "post-reassignment"),
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05,
                    main = "Cluster dendrogram of WGCNA modules")
dev.off()

# ---- merge KME data with DESeq separate timepoint results ----
kme <- read.csv(paste0(netDir, "kME.csv"))
tpDE <- read.csv(paste0(resultDir, "DE/timepoints/mergedDESepTimepointsWTvsKO.csv")) %>%
  filter(geneIDs %in% kme$gene_id)
modules <- read.csv(paste0(netDir, "modules.csv"))


rownames(tpDE) <- tpDE$geneIDs
rownames(kme) <- kme$gene_id
# for all the non-DE genes, we will use pval = 1 and fc = 0 as a placeholder
nUnsig <- nrow(kme) - nrow(tpDE)
unsig <- data.frame(matrix(nrow = nUnsig, ncol = ncol(tpDE)))
rownames(unsig) <- rownames(kme)[!rownames(kme) %in% rownames(tpDE)]
colnames(unsig) <- colnames(tpDE)
unsig[,"geneIDs"] <- rownames(unsig)
unsig[,"geneNames"] <- IDToName[rownames(unsig)]
unsig[, grep("padj", colnames(unsig))] <- rep(1)
unsig[, grep("log2FC", colnames(unsig))] <- rep(0)
filledTpDE <- rbind(tpDE, unsig)
mergedKME <- cbind(kme[match(filledTpDE$geneIDs, kme$gene_id),], filledTpDE)

# check there is no mismatch
sum(mergedKME$gene_id != mergedKME$geneIDs)
# some mismatch with names, but not IDs
mergedKME[mergedKME$gene_name != mergedKME$geneNames, c("gene_name", "geneNames")]
mergedKME <- mergedKME %>% dplyr::select(-c("geneIDs", "geneNames"))
fwrite(mergedKME, paste0(netDir, "kMEMergedDERes.csv"))

allMod <- gsub("E", "", modules$Label)
# check module names line up
allMod %in% unique(mergedKME$moduleLabel)
# ncol = 9 for (odds ratio, pval, padj) for each of the 3 timepoints
DEKMECor <- data.frame(matrix(nrow = length(allMod), ncol = 9))
rownames(DEKMECor) <- allMod
colnames(DEKMECor) <- paste(rep(c("D20", "D60", "D100"), each = 3),
                            rep(c("oddsRatio", "p-value", "p-adjusted"), 3), sep = ".")
for (m in allMod) {
  for (day in c("D20", "D60", "D100")) {
    padj <- paste0(day, ".padj") # name of kme column with DE p adjusted results
    cor <- fisher.test(mergedKME$moduleLabel == m, mergedKME[[padj]] < 0.05)
    DEKMECor[m, paste0(day, ".oddsRatio")] <- cor$estimate
    DEKMECor[m, paste0(day, ".p-value")] <- cor$p.value
  }
}

# apply BH pval adjustment
for (pvalCol in grep("p-value", colnames(DEKMECor))) {
  # check the column selected is the padjusted one
  stopifnot(grepl("p-adjusted", colnames(DEKMECor[, pvalCol + 1])))
  DEKMECor[, pvalCol + 1] <- p.adjust(DEKMECor[, pvalCol], method = "BH")
}

fwrite(DEKMECor, paste0(netDir, "DEResToModCorrelation.csv"), row.names = TRUE)

# rematching with up and down regulated genes separately
DEKMECorFull <- data.frame(matrix(nrow = length(allMod), ncol = 8*3))
rownames(DEKMECorFull) <- allMod
colnames(DEKMECorFull) <- paste(rep(c("D20", "D60", "D100"), each = 8),
                              rep(c("up", "down"), each = 4),
                              rep(c("oddsRatio", "p-value", "p-adjusted", "%in"), 3),
                            sep = ".")
nDEGenes <- list()
for (m in allMod) {
  for (day in c("D20", "D60", "D100")) {
    padjCol <- paste0(day, ".padj") # name of kme column with DE p adjusted results
    logFCCol <- paste0(day, ".log2FC") 
    
    sigUpreg <- (mergedKME[[padjCol]] < 0.05 & mergedKME[[logFCCol]] > log2(1.5))
    nDEGenes[[day]]$up <- sum(sigUpreg)
    
    corUp <- fisher.test(mergedKME$moduleLabel == m, sigUpreg)
    DEKMECorFull[m, paste0(day, ".up.oddsRatio")] <- corUp$estimate
    DEKMECorFull[m, paste0(day, ".up.p-value")] <- corUp$p.value
    
    # e.g. M1, 38/2875
    DEKMECorFull[m, paste0(day, ".up.%in")] <- nrow(mergedKME[mergedKME$moduleLabel == m & sigUpreg,])*100/nrow(mergedKME[mergedKME$moduleLabel == m,])
    
    sigDownreg <- (mergedKME[[padjCol]] < 0.05 & mergedKME[[logFCCol]] < -log2(1.5))
    nDEGenes[[day]]$down <- sum(sigDownreg)
    
    corDown <- fisher.test(mergedKME$moduleLabel == m, sigDownreg)
    DEKMECorFull[m, paste0(day, ".down.oddsRatio")] <- corDown$estimate
    DEKMECorFull[m, paste0(day, ".down.p-value")] <- corDown$p.value
    
    # e.g. M1, 303/2875
    DEKMECorFull[m, paste0(day, ".down.%in")] <- nrow(mergedKME[mergedKME$moduleLabel == m & sigDownreg,])*100/nrow(mergedKME[mergedKME$moduleLabel == m,])
    
    print(corDown)
  }
}

# apply BH pval adjustment
for (pvalCol in grep("p-value", colnames(DEKMECorFull))) {
  # check the column selected is the padjusted one
  stopifnot(grepl("p-adjusted", colnames(DEKMECorFull[, pvalCol + 1])))
  DEKMECorFull[, pvalCol + 1] <- p.adjust(DEKMECorFull[, pvalCol], method = "BH")
}
fwrite(DEKMECorFull, paste0(netDir, "DEResToModCorrelationUpDownRegMinFC.csv"), row.names = TRUE)


# check the samples in meta data and MEs line up
metaNet$Sample == rownames(me)
# sort rows (samples) by day and cols (MEs) by number of genes (except grey at end)
sortedMEs <- me[order(metaNet$numericDay),
                modules$Label[c(2:nrow(modules), 1)]]
sortedMEs <- sortedMEs[c(grep("WT", rownames(sortedMEs)), grep("KO", rownames(sortedMEs))), ]

# create list of names containing day and batch
names <- gsub("D2._","D20_", rownames(sortedMEs)) %>%
  str_extract(regex("B\\d_D\\d*")) %>%
  str_replace_all("_", "\n")

colorToME <- modules$Label
names(colorToME) <- modules$Color

METoColor <- names(colorToME)
names(METoColor) <- colorToME

# ---- make heatmap of percentages ---- 
data <- read.csv(paste0(netDir, "DEResToModCorrelationUpDownRegMinFC.csv"))
mepvals <- read.csv(paste0(netDir,"padjustedMEs.csv"))
rownames(mepvals) <- gsub("ME", "M", mepvals[,1])
mepvals <- mepvals[,-1]

rownames(data) <- data[,1]
data <- data[,-1]
# no need to show grey mod
data <- data[2:18, c(grep("down", colnames(data)), grep("up", colnames(data)))]

pvals <- data[, grep("value", colnames(data))]
pvals <- pvals*18
or <- data[, grep("oddsRatio", colnames(data))]
perc <- data[, grep("..in", colnames(data))]

plotdata <- perc
mepvals <- mepvals[match(rownames(plotdata), rownames(mepvals)), ]
mepvals <- mepvals[, grep("Genotype", colnames(mepvals))]
for(j in c(1:ncol(mepvals)))
  mepvals[,j] <- ifelse(mepvals[,j] <0.05, 1, 0)
labeldata <- plotdata; labeldata[,]=" "

for (j in c(1:ncol(labeldata)))
{ 
  sig <- which((or[,j]>1)&(pvals[,j]<0.05))
  labeldata[sig,j] <- "*"
}  
pdf(paste0(netDir, "uneditedDEToWGCNAFig.pdf"))
pheatmap(plotdata,angle_col = 0,
         annotation_row=mepvals, display_numbers=labeldata)
dev.off()

# ---- 3. Run linear anova model with ~ Batch + Genotype + Day ----
MEs <- read.csv(paste0(netDir, "ME.csv"))
a <- list()
MEPval <- data.frame(matrix(nrow = 0, ncol = 4))
colnames(MEPval) <- c("Batch", # batch not needed for presentation
                      "Day", "Genotype", "GenotypeDayInteraction")
fullMECorData <- list()
for (curME in colnames(MEs[-1])) {
  lMod <- lm(MEs[, curME] ~  metaNet$Batch + metaNet$numericDay + metaNet$Genotype 
             + metaNet$numericDay:metaNet$Genotype)
  fullMECorData[[curME]] <- lMod
  a[[curME]] <- anova(lMod)
  curMEPval <- data.frame(Batch = a[[curME]]['metaNet$Batch', 'Pr(>F)'], # batch not needed for presentation
                      Day = a[[curME]]['metaNet$numericDay', 'Pr(>F)'],
                      Genotype = a[[curME]]['metaNet$Genotype', 'Pr(>F)'],
                      GenotypeDayInteraction = a[[curME]]['metaNet$numericDay:metaNet$Genotype', 'Pr(>F)'])
  rownames(curMEPval) <- curME
  MEPval <- rbind(MEPval, curMEPval)
}

saveRDS(fullMECorData, paste0(netDir, "fullMELMCorData.rds"))

MEPadj <- MEPval
for (var in colnames(MEPval)) {
  MEPadj[[var]] <- p.adjust(MEPval[[var]], method = "BH")
}
fwrite(MEPval, paste0(netDir, "pvaluesMEs.csv"), row.names = TRUE)
fwrite(MEPadj, paste0(netDir, "padjustedMEs.csv"), row.names = TRUE)

# ---- plot eigengenes ---- 
# ggplot version:
# Find cell function from http://127.0.0.1:15985/library/gridExtra/doc/tableGrob.html
# Allows style editing of individual cells of a table grob
find_cell <- function(table, row, col, name="core-fg"){
  l <- table$layout
  which(l$t==row & l$l==col & l$name==name)
}

# custom theme for the pval table
themeSmallText <- gridExtra::ttheme_default(
  core = list(fg_params = list(cex = 0.8)),
  colhead = list(fg_params = list(cex = 0.8)),
  rowhead = list(fg_params = list(cex = 0.8)))

kme <- read.csv(paste0(netDir, "kME.csv"))
modules <- modules <- read.csv(paste0(netDir, "modules.csv"))
allMod <- gsub("E", "", modules$Label)


# creating a custom legend:
leg <- data.frame(WT = c("darkolivegreen4", "darkolivegreen3", "darkolivegreen"),
           KO = c("indianred","indianred2", "indianred4"))
rownames(leg) <- c("D20", "D60", "D100")
legEmpty <- leg
legEmpty[,] <- ""
grobLeg <- tableGrob(legEmpty, theme = ttheme_default())
for (row in c(2:4)) {
  for (col in c(2:3)) {
    boxInd <- find_cell(grobLeg, row, col, "core-bg")
    textInd <- find_cell(grobLeg, row, col, "core-fg")
    grobLeg$grobs[[boxInd]]$gp$fill <- leg[row-1, col-1]
  }
}

pdf(file =  paste0(netDir, "sortedEigengenesGgplotPadjDayColour3.pdf"), height=3, width=12)
EGPlots3 <- list()
dayColors <- c(D20 = "darkblue", D60 = "darkorchid4", D100 = "deeppink4")

for (curME in colnames(sortedMEs[,1:18])) {
  n <- modules$nGenes[match(curME, modules$Label)]
  MEData <- sortedMEs %>%
    dplyr::select(all_of(curME)) %>%
    mutate(sampleLabel = names,
           Genotype = metaNet$Genotype[match(row.names(sortedMEs), metaNet$Sample)],
           Day = metaNet$Day[match(row.names(sortedMEs), metaNet$Sample)],
           GenotypeDay = paste(Genotype, Day))
  
  # main eigengene plot
  # factoring of row names keeps their original order
  p <- ggplot(MEData, aes(x = factor(rownames(MEData), levels = rownames(MEData)), y = MEData[,curME], fill = GenotypeDay)) + 
    scale_x_discrete(labels = MEData$sampleLabel) +
    #scale_fill_manual(breaks = c("WT", "KO"), values = c("darkgreen", "indianred")) +
    scale_fill_manual(breaks = c("WT D20", "WT D60", "WT D100",
                                 "KO D20", "KO D60", "KO D100"), 
                      values = c(leg$WT, leg$KO), 
                      labels = c("WT D20", "       D60", "       D100",
                                 "KO D20", "      D60", "      D100")) +
    labs(x = "Sample", y = "Eigengene value", 
         title = paste0("Module ", gsub("ME", "", curME), "; NGenes: ", n)) +
    theme(axis.text.x = element_text(size = 10
                                     #, color = dayColors[fct_inorder(MEData$Day)]
                                     ), 
          plot.title = element_text(size = 18),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12))
  
  
  p <- p + geom_bar(stat = "identity")
  cutoffs <- c(4.5, 7.5, 9.5)
  for (linePos in c(cutoffs, cutoffs + 9)[-6]) {
    p <- p + geom_vline(xintercept = linePos, linetype = "dashed", alpha = 0.4)
  }
  p
  # set up the pvalue table on the side of the eigengene plot
  # Batch is removed at this step
  curPadj <- signifNumeric(MEPadj[curME, colnames(MEPadj) != "Batch"], 3) %>%
    mutate("Genotype-Day\nInteraction" = GenotypeDayInteraction) %>%
    dplyr::select(-GenotypeDayInteraction) %>% t()
  colnames(curPadj) <- "P-adjusted"
  sigVarNums <- which(curPadj[,"P-adjusted"] < 0.05)
  plainMEPlot <- p + theme(legend.position = "none")
  pvalGrob <- tableGrob(curPadj, theme = ttheme_default())
  # highlight/bold the significant variables
  for (sigVar in sigVarNums) {
    # the + 1 accounts for the first row being the title
    ind <- find_cell(pvalGrob, sigVar + 1, 2, "core-bg")
    rowNameInd <- find_cell(pvalGrob, sigVar + 1, 1, "rowhead-fg")
    pvalGrob$grobs[ind][[1]][["gp"]] <- gpar(fill="khaki", col = "white")
    pvalGrob$grobs[rowNameInd][[1]][["gp"]] <- gpar(fontface="bold")
    
  }
  MEInfo <- ggarrange(grobLeg, pvalGrob, nrow = 2)
  full <- ggarrange(plainMEPlot, MEInfo, widths = c(3, 1))
  print(full)
  EGPlots3[[curME]]$plot <- plainMEPlot
  EGPlots3[[curME]]$info <- MEInfo
  EGPlots3[[curME]]$legend <- get_legend(p)
  EGPlots3[[curME]]$full <- full
}
dev.off()

ggarrange(EGPlots3[["ME1"]]$full, EGPlots3[["ME2"]]$full, nrow = 2,
          labels = "AUTO", font.label = list(size = 18))
ggsave(filename = paste0(netDir, "thesisFig/M1-2Eigenplot.png"), 
       width = 12, height = 6, device='png', dpi=300, bg = "white")
dev.off()

ggarrange(EGPlots3[["ME3"]]$full,
          EGPlots3[["ME4"]]$full,
          EGPlots3[["ME13"]]$full,
          nrow = 3, labels = "AUTO", font.label = list(size = 20))
ggsave(filename = paste0(netDir, "thesisFig/M3-4-13Eigenplot.png"), 
       width = 12, height = 9, device='png', dpi=300, bg = "white")
dev.off()


ggarrange(EGPlots3[["ME8"]]$full,
          EGPlots3[["ME10"]]$full,
          nrow = 2, labels = "AUTO", font.label = list(size = 20))
ggsave(filename = paste0(netDir, "thesisFig/M8-10Eigenplot.png"), 
       width = 12, height = 6, device='png', dpi=300, bg = "white")
dev.off()


pdf(file =  paste0(netDir, "sortedEigengenesGgplotPadjDayColour4.pdf"), height=3, width=12)
EGPlots4 <- list()
dayColors <- c(D20 = "darkblue", D60 = "darkorchid4", D100 = "deeppink4")

for (curME in colnames(sortedMEs[,1:18])) {
  n <- modules$nGenes[match(curME, modules$Label)]
  MEData <- sortedMEs %>%
    dplyr::select(all_of(curME)) %>%
    mutate(sampleLabel = names,
           Genotype = metaNet$Genotype[match(row.names(sortedMEs), metaNet$Sample)],
           Day = metaNet$Day[match(row.names(sortedMEs), metaNet$Sample)],
           GenotypeDay = paste(Genotype, Day))
  
  # main eigengene plot
  # factoring of row names keeps their original order
  p <- ggplot(MEData, aes(x = factor(rownames(MEData), levels = rownames(MEData)), y = MEData[,curME], fill = Genotype)) + 
    scale_x_discrete(labels = MEData$sampleLabel) +
    scale_fill_manual(breaks = c("WT D20", "WT D60", "WT D100",
                                 "KO D20", "KO D60", "KO D100"), 
                      values = c("darkolivegreen", "darkolivegreen2", "darkolivegreen4",
                                 "indianred4","indianred2", "indianred4"), 
                      labels = c("WT D20", "       D60", "       D100",
                                 "KO D20", "      D60", "      D100")) +
    labs(x = "Sample", y = "Eigengene value", 
         title = paste0("Module ", gsub("ME", "", curME), "; NGenes: ", n)) +
    theme(axis.text.x = element_text(size = 10), 
    plot.title = element_text(size = 18),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12))
  
  virdisCols <- viridis(3, alpha = 0.1, begin = 0, end = 1, direction = 1, option = "C")
  for (annValues in list(c(0.5, 4.5, virdisCols[1]), c(4.5, 7.5, virdisCols[2]), c(7.5, 9.5, virdisCols[3]))) {
    p <- p + annotate("rect", xmin = as.numeric(annValues[1]), 
                      xmax = as.numeric(annValues[2]), ymin = -Inf, ymax = Inf,
                      alpha = 0.1, fill = annValues[3]) +
      annotate("rect", xmin = as.numeric(annValues[1]) + 9, 
               xmax = as.numeric(annValues[2])+ 9, ymin = -Inf, ymax = Inf,
               alpha = 0.1, fill = annValues[3]) 
  }
  p <- p + geom_bar(stat = "identity")
  p
  # set up the pvalue table on the side of the eigengene plot
  # Batch is removed at this step
  curPadj <- signifNumeric(MEPadj[curME, colnames(MEPadj) != "Batch"], 3) %>%
    mutate("Genotype-Day\nInteraction" = GenotypeDayInteraction) %>%
    dplyr::select(-GenotypeDayInteraction) %>% t()
  colnames(curPadj) <- "P-adjusted"
  sigVarNums <- which(curPadj[,"P-adjusted"] < 0.05)
  plainMEPlot <- p + theme(legend.position = "none")
  pvalGrob <- tableGrob(curPadj, theme = ttheme_default())
  # highlight/bold the significant variables
  for (sigVar in sigVarNums) {
    # the + 1 accounts for the first row being the title
    ind <- find_cell(pvalGrob, sigVar + 1, 2, "core-bg")
    rowNameInd <- find_cell(pvalGrob, sigVar + 1, 1, "rowhead-fg")
    pvalGrob$grobs[ind][[1]][["gp"]] <- gpar(fill="khaki", col = "white")
    pvalGrob$grobs[rowNameInd][[1]][["gp"]] <- gpar(fontface="bold")
    
  }
  MEInfo <- ggarrange(get_legend(p), pvalGrob, nrow = 2)
  print(ggarrange(plainMEPlot, MEInfo, widths = c(3, 1)))
  EGPlots4[[curME]]$plot <- plainMEPlot
  EGPlots4[[curME]]$info <- MEInfo
  EGPlots4[[curME]]$legend <- get_legend(p)
}
dev.off()



# ---- 4. Go analysis of WGCNA output ----
cleanGOBarplot <- function(GORes, showCategory = 8, title = "", ylab = "",
                           xlab = "-log10 P-adjusted", hjust = 1.5, textWidth = 25) {
  GOSubset <- GORes@result[1:showCategory,] %>% filter(p.adjust < 0.05)
  ggplot(GOSubset[order(GOSubset$p.adjust, decreasing = TRUE),],
                                      aes(x = -log10(p.adjust),
                                          y = fct_inorder(Description),
                                          fill = log10(p.adjust))) +
    geom_bar(stat = "identity") +
    scale_y_discrete(labels = function(x) str_wrap(str_replace_all(x, " " , "\n"),
                                                   width = textWidth)) +
    labs(y = ylab, x = xlab,
         title = paste0("      ", title)) + 
    theme(legend.position = "none",
          plot.title.position = "plot",
          axis.title.x = element_text(hjust = hjust)
    )
  
}
 
allColors <- modules$Color[c(2:nrow(modules), 1)]
GOResultsFull <- list()

#cex1 = 0.9;
data(geneList, package="DOSE")
for (color in allColors) {
  inCurMod <- kme %>% subset(moduleColor %in% color)
  namesCurMod <- inCurMod$gene_id
  for (ontology in c("BP", "MF", "CC")) {
    GOResultsFull[[ontology]][[color]]  <- enrichGO(namesCurMod,
                                                    universe = colnames(WGCNA_TPM),
                                                    OrgDb = "org.Hs.eg.db",
                                                    keyType = "ENSEMBL",
                                                    ont = ontology,
                                                    pvalueCutoff = 0.05, # default
                                                    qvalueCutoff = 0.05, # 
                                                    minGSSize = 10, # default
                                                    maxGSSize = 500 # default)
    GOSimple[[ontology]][[color]] <- simplify(GOResultsFull[[ontology]][[color]],
                                              cutoff = 0.7, by = "p.adjust",
                                              select_fun = min)
    print(paste(color, inCurMod[1, "moduleLabel"], ontology))
  }
}

# to conserve compatability with previous code that ran on only GO BP:
GOResults <- GOResultsFull$BP

GOSimple <- list()
for (color in allColors) {
  GOSimple[[color]] <- simplify(GOResults[[color]], cutoff=0.7, by="p.adjust", select_fun=min)
}

# save top GO terms and in their simple and complex versions 
if (!dir.exists(paste0(netDir, "GO"))) {dir.create(paste0(netDir, "GO"))}
for (ont in c("BP", "MF", "CC")) {
  for (complexity in c("Full", "Simplified")) {
    pdf(file =  paste0(netDir, "GO/topGOTerms", ont, complexity, ".pdf"))
    par(mfrow = c(1,1))
    for (color in allColors) {
      curGO <- GOResultsFull[[ont]][[color]]
      if (complexity == "Simplified") {
        curGO <- GOSimple[[ont]][[color]]
      }
      # exclude modules with no significant terms
      if (sum(curGO@result$p.adjust < 0.05) > 0) {
        print(barplot(curGO, showCategory=15,
                      title = paste(colorToME[color], color,
                                    "module genes - GO", ont,
                                    str_to_lower(complexity))))
      }
    }
    dev.off()
  }
}

saveRDS(GOResultsFull, file = paste0(netDir, "allGOResultsBPMFCC.rds"))
saveRDS(GOResults, file = paste0(netDir, "allGOResults.rds"))
saveRDS(GOSimple, file = paste0(netDir, "allGOResultsSimplified.rds"))
#GOSimple2 <- readRDS( paste0(netDir, "allGOResultsSimplified.rds"))
#cleanGOBarplot(GOSimple2$turquoise, 14)

if (!exists("GOResults")) {
  GOResults <- readRDS(file =  paste0(netDir, "allGOResults.rds")) 
}

print(cleanGOBarplot(GOResults$turquoise))

color <- "lightcyan"
pdf(file =  paste0(netDir, "lightcyanModExploratoryCnetplots.pdf"))
cnetplot(GOResults[[color]], showCategory = 2, 
         cex.params = list(gene_label = 0.5))
cnetplot(GOResults[[color]], showCategory = 8, 
         cex.params = list(gene_label = 0.1, category_label=0.3))
cnetplot(GOResults[[color]], showCategory = "dorsal/ventral pattern formation")
dev.off()

#heatmap plots for different modules
module_order <- names(sortedMEs) %>% gsub("ME","", .)
sortedMEs$sample <- row.names(sortedMEs)
orderChange <- match(sortedMEs$sample, gsub("D2.","D20", metaNet$Sample))
sortedMEs$group <- metaNet$Genotype[orderChange]
sortedMEs$batch <- metaNet$Batch[orderChange]
sortedMEs$day <- metaNet$Day[orderChange]
sortedMEs$numeric_day <- as.numeric(sortedMEs$day %>% gsub("D","", .))
mME <- sortedMEs %>%
  pivot_longer(cols=starts_with("ME")) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )

if(!dir.exists(paste0(netDir, "heatmaps/"))) {dir.create(paste0(netDir, "heatmaps/"))}
pdf(file = paste0(netDir, "heatmaps/metadataVsModuleHeatmaps.pdf"))
# make sure day order is numerically sorted, not 100, 20, 60
dayPlot <-  ggplot(mME, aes(x=factor(day, levels=c("D20", "D60", "D100")), y=name, fill=value)) +
  labs(title = "Module-day Relationships", x = "Day", y = "Modules", fill="corr")
genotypePlot <- ggplot(mME, aes(x=factor(group, levels=c("WT", "KO")), y=name, fill=value)) +
  labs(title = "Module-Experimental Group Relationships", x = "Genotype", y = "Modules", fill="corr")
batchPlot <- ggplot(mME, aes(x=batch, y=name, fill=value)) +
  labs(title = "Module-Batch Relationships", x = "Batch", y = "Modules", fill="corr")
for (plot in c(dayPlot, genotypePlot, batchPlot)) {
  print(plot + geom_tile() +
          theme_bw() + 
          scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                               midpoint = 0, limit = c(-1,1)) +
          theme(axis.text.x = element_text(angle=90)))
}
dev.off()

# ---- more useful heatmap - module expression heatmap with a sidebar containing DE membership ----
DE <- readRDS(paste0(resultDir, "DE/timepoints/timeSplitDEResults.rds"))
kme <- read.csv(paste0(netDir, "kME.csv"))
modules <- modules <- read.csv(paste0(netDir, "modules.csv"))
allMod <- gsub("E", "", modules$Label)

# heamap row scaling function from https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/
cal_z_score <- function(x) {
  (x - mean(x)) / sd(x)
}

mod <- "M17"
# Specify colors
library(viridis)
library("colorspace")

turbo(2)
rocket(2)
"#440154FF" "#FDE725FF"
[1] "#000004FF" "#FCFFA4FF"
"#00204DFF" "#FFEA46FF"
"#0D0887FF" "#F0F921FF"
[1] "#30123BFF" "#7A0403FF"
"#03051AFF" "#FAEBDDFF"
[1] "#0B0405FF" "#DEF5E5FF"

plasma(2)
dayColours <- c(up = "#F0F921FF", down = "#0D0887FF", unchanged = "lightgrey")
annColors = list(
  D20 = dayColours,
  D60 = dayColours,
  D100 = dayColours,
  Genotype = c(WT = "darkgreen", KO = "indianred"),
  DifferentiationStage = c("2" = "grey10", "3" = "grey25", "4" = "grey44", "5" = "grey65", "7" = "grey80"),
  Day = c(D20 = "skyblue", D60 = "darkseagreen3", D100 = "palegreen4")
  
)
colAnn <- data.frame(DifferentiationStage = gsub("B", "", metaNet$Batch),
                     Day = factor(metaNet$Day, levels = c("D20", "D60", "D100")),
                     Genotype = metaNet$Genotype)
rownames(colAnn) <- metaNet$Sample

heatmaps <- list()
unannHeatmaps <- list()
pdf(file = paste0(netDir, "heatmaps/allModDEHeatmaps2.pdf"), width = 8, height = 6)
for (mod in allMod[c(2:18, 1)]) {
  #mod <- "M3"
  color <- modules$Color[gsub("E", "", modules$Label) == mod]
  modGenes <- kme$gene_id[kme$moduleLabel == mod]
  modExpRaw <- TPM_filtered[TPM_filtered$gene_id %in% modGenes,]
  modExpRaw$kme <- kme[[mod]][match(modExpRaw$gene_id, kme$gene_id)]
  # sorted
  modExp <- modExpRaw[order(modExpRaw$kme, decreasing = TRUE),] %>%
    remove_rownames() %>%
    column_to_rownames("gene_id") %>%
    dplyr::select(2:19) %>%
    as.matrix()
  modExpNorm <- t(apply(modExp, 1, cal_z_score))
  
  # row annotation by DE results:
  rowAnn <- data.frame(D100 = rep("unchanged", length(modGenes)),
                       D60 = rep("unchanged", length(modGenes)),
                       D20 = rep("unchanged", length(modGenes)))
  rownames(rowAnn) <- modGenes
  for (day in c("D20", "D60", "D100")) {
    curDE <- DE[[day]]$WTvsKO %>%
      as.data.frame() %>%
      filter(rownames(.) %in% modGenes) %>%
      filter(padj < 0.05)
    curDown <- curDE %>% filter(log2FoldChange < -log2(1.5)) %>% rownames()
    curUp <- curDE %>% filter(log2FoldChange > log2(1.5)) %>% rownames()
    rowAnn[curUp, day] <- "up"
    rowAnn[curDown, day] <- "down"
    rowAnn[[day]] <- factor(rowAnn[[day]], levels = c("up", "down", "unchanged"))
  }
  # sort by day before batch
  dayOrd <- c(grep("D2.", colnames(TPM_filtered)), grep("D60", colnames(TPM_filtered)), grep("D100", colnames(TPM_filtered)))
  genOrd <- c(grep("WT", colnames(TPM_filtered)[dayOrd]), grep("KO", colnames(TPM_filtered)[dayOrd]))

  mainHeatmap <- pheatmap(modExpNorm[, colnames(TPM_filtered)[dayOrd][genOrd]],
                          annotation_row = rowAnn, annotation_col = colAnn,
                          annotation_colors = annColors,
                          main = paste("Module", gsub("M", "", mod), "gene expression"),
                          show_rownames = FALSE, cluster_rows = FALSE,
                          show_colnames = FALSE, cluster_cols = FALSE,  silent = TRUE)
  
  # get combined legend for DE results (D20, D60, D100)
  leg <- list()
  leg$annColors <- annColors[!names(annColors) %in% c("D20", "D60", "D100")] 
  leg$annColors[["DE results\n(D20, D60, D100)"]] <- annColors$D20
  leg$rowAnn <- rowAnn[!names(rowAnn) %in% c("D60", "D100")]
  rownames(leg$rowAnn) <- rownames(rowAnn)
  colnames(leg$rowAnn) <- "DE results\n(D20, D60, D100)"
  
  legHeatmap <- pheatmap(modExpNorm[, colnames(TPM_filtered)[dayOrd][genOrd]], annotation_row = leg$rowAnn, annotation_col = colAnn,
                         annotation_colors = leg$annColors,
                         main = paste0("Module ", gsub("M", "", mod)), show_rownames = FALSE, cluster_rows = FALSE,
                         show_colnames = FALSE, cluster_cols = FALSE, silent = TRUE)
  # replace legend with the DE results (D20, D60, D100) legend
  mainHeatmap$gtable$grobs[[7]] <- legHeatmap$gtable$grobs[[7]]
  # move the DE results legend down a little to avoid overlap
  mainHeatmap$gtable$grobs[[7]]$children[[11]]$vjust <- mainHeatmap$gtable$grobs[[7]]$children[[11]]$vjust + 1.2
  mainHeatmap$gtable$grobs[[7]]$children[[12]]$vjust <- mainHeatmap$gtable$grobs[[7]]$children[[12]]$vjust + 2.2
  
  grid.newpage()
  print(mainHeatmap)
  
  heatmaps[[mod]] <- mainHeatmap
  
  unannHeatmap <- pheatmap(modExpNorm[, colnames(TPM_filtered)[dayOrd][genOrd]],
                           annotation_row = rowAnn, annotation_col = colAnn[colnames(colAnn) != "DifferentiationStage"],
                           annotation_colors = annColors[names(annColors) != "DifferentiationStage"],
                           main = paste("Module", gsub("M", "", mod)),
                           show_rownames = FALSE, cluster_rows = FALSE,
                           show_colnames = FALSE, cluster_cols = FALSE, silent = TRUE,
                           annotation_legend = FALSE)
  unannHeatmaps[[mod]] <- unannHeatmap
}
dev.off()

# z score is too time consuming at the moment - skip for now
zScore <- unannHeatmaps$M4$gtable$grobs[[7]]
zScore$children[[1]]$x <- zScore$children[[1]]$x + unit(0.5, "npc")
zScore$children[[2]]$x <- zScore$children[[2]]$x + unit(0.5, "npc")
zScoreGrid <- ggarrange(textGrob("Z-score", vjust = 8), zScore, nrow = 2,)


for (mod in c("M3", "M8")) {
  unannHeatmaps[[mod]]$gtable$grobs[[4]] <- grob()
  unannHeatmaps[[mod]]$gtable$grobs[[7]] <- grob()
}
unannHeatmaps$M10$gtable$grobs[[7]] <- grob()

legend <- mainHeatmap$gtable$grobs[[7]]$children
# line up legend along y axis
for (i in c(1, 4, 7, 10)) {
  legend[[i + 1]]$y <- legend[[i + 1]]$y - legend[[i]]$y[1] + unit(0.9, "npc")
  legend[[i + 2]]$y <- legend[[i + 2]]$y - legend[[i]]$y[1] + unit(0.9, "npc")
  legend[[i]]$y <- unit(0.9, "npc")
}

# the newline in DE results\n(D20, D60, D100) is not needed in this format
legend[[10]]$label <- "DE results (D20, D60, D100)"
# undo the downward shift needed to avoid overlap with a newline
legend[[11]]$vjust <- legend[[11]]$vjust - 1.2
legend[[12]]$vjust <- legend[[12]]$vjust - 2.2

col1 <- ggarrange(unannHeatmaps$M3$gtable, unannHeatmaps$M4$gtable, grob(),
                  ncol = 3, widths = c(10, 10, 1), labels = c("A", "B"),
                  font.label = list(size = 18))
col2 <- ggarrange(unannHeatmaps$M8$gtable, unannHeatmaps$M10$gtable, grob(),
                  ncol = 3, widths = c(10, 10, 1), labels = c("C", "D"),
                  font.label = list(size = 18))
col3 <- ggarrange(grob(), legend[c(1:3)], legend[c(4:6)], legend[c(10:12)], 
          ncol = 4, nrow = 1, widths = c(1, 4, 4, 5))
ggarrange(col1, 
          col2,
          col3,
          nrow = 3, heights = c(3, 3, 1))

ggsave(paste0(netDir, "thesisFig/Mods3-4-8-10Heatmaps.png"), device = "png",
       height = 7, width = 6, bg = "white")
dev.off()



thesisMod <- c("M1", "M2", "M3", "M4", "M8", "M10", "M13")
for (mod in thesisMod) {
  grid.newpage()
  show(heatmaps[[mod]])
  ggsave(paste0(netDir, "thesisFig/", mod ,"heatmap.png"), device = "png",
         plot = heatmaps[[mod]], height = 5, width = 6, bg = "white")
  dev.off()
}
print(mainHeatmap)

thesisHeatmaps <- heatmaps[thesisMod]
ggarrange(thesisHeatmaps[[1]])
b <- gList(thesisHeatmaps[[1]]$gtable$respect)
ggarrange(b)


# ---- 10. KEGG and WikiPathways enrichment - focus on only WT vs KO ----
#data(geneList, package = "DOSE")
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
hgncEntrez <- getBM(filters = "ensembl_gene_id", 
                    attributes = c("ensembl_gene_id","entrezgene_id"),
                    values = names(IDToName), mart = mart)
IDToEntrez <- hgncEntrez$entrezgene_id
names(IDToEntrez) <- hgncEntrez$ensembl_gene_id
entrezToID <- hgncEntrez$ensembl_gene_id
names(entrezToID) <- hgncEntrez$entrezgene_id
entrezToName <- IDToName[hgncEntrez$ensembl_gene_id]
names(entrezToName) <- hgncEntrez$entrezgene_id

pathwayRes <- list()
background <- IDToEntrez[colnames(WGCNA_TPM)]
print(paste("Lost background genes:", sum(is.na(background))))
# 6413 genes lost due to not having an entrez ID
background <- background[!is.na(background)]
for (color in allColors) {
  inCurMod <- kme %>%
    subset(moduleColor %in% color)
  idCurMod <- inCurMod$gene_id
  netEntrez <- IDToEntrez[idCurMod]
  print(paste("Lost genes:", sum(is.na(netEntrez))))
  netEntrz <- netEntrez[!is.na(netEntrez)]
  pathwayRes$KEGG[[color]] <- enrichKEGG(netEntrez,
                                         universe = background,
                                         organism = "hsa",
                                         keyType = "ncbi-geneid",
                                         pvalueCutoff = 0.05,
                                         qvalueCutoff = 0.05,
                                         minGSSize = 10, maxGSSize = 500)
  
  pathwayRes$WP[[color]] <- enrichWP(gene = netEntrez,
                                     organism = "Homo Sapiens",
                                     universe = background)
}

entrezToNameFun <- function(IDs) {
  unname(entrezToName[IDs])
}
entrezToIDFun <- function(IDs) {
  unname(entrezToID[IDs])
}
IDToNameFun <- function(IDs) {
  unname(IDToName[IDs])
}

fullPathway <- list()
pathwayMethods <- c("KEGG", "WP")
if (!dir.exists(paste0(netDir, "pathways"))) {
  dir.create(paste0(netDir, "pathways"))
}

if (!exists("pathwayRes")) {
  pathwayRes <- list()
  pathwayRes$KEGG <- readRDS(paste0(netDir, "pathways/fullKEGGWGCNA.rds"))
}
colorToME <- modules$Label
names(colorToME) <- modules$Color
METoColor <- names(colorToME)
names(METoColor) <- colorToME

for (method in pathwayMethods) {
  pdf(file = paste0(netDir, "pathways/top", method, "WGCNA.pdf"))
  par(mfrow = c(1,1))
  for (color in allColors) {
    print(color)
    # some pathways have no significant terms - causes error if not filtered
    if (pathwayRes[[method]][[color]]@result$p.adjust[1] < 0.05) {
      print(barplot(pathwayRes[[method]][[color]], showCategory=15, 
                    title = paste0(method, " - WGCNA module ", color, " (", colorToME[color], ")")))
    }
    
  }
  dev.off()
  
  saveRDS(pathwayRes[[method]], file = paste0(netDir, "pathways/full",
                                              method, "WGCNA.rds"))
  
  # save as csv output
  if (!dir.exists(paste0(netDir, "pathways/", method ,"CSV"))) {
    dir.create(paste0(netDir, "pathways/", method ,"CSV"))
  }
  for (color in allColors) {
    # add gene IDs and names based on "/"-separated ensemblIDs 
    curRes <- pathwayRes[[method]][[color]]@result %>%
      mutate(geneNames = lapply(lapply(str_split(geneID, "/"), entrezToNameFun),toString),
             ensembIDs = lapply(lapply(str_split(geneID, "/"), entrezToIDFun), toString))
    fwrite(curRes, paste0(netDir, "pathways/", method ,"CSV/",
                          color, method, "Res.csv"))
    fullPathway[[method]][[color]] <- curRes
  }
}

## repeat the same for GO
# save as csv output
if (!dir.exists(paste0(netDir, "GO/fullCSV/"))) {
  dir.create(paste0(netDir, "GO/fullCSV/"))
}

if (!exists("GOResultsFull")) {
  GOResultsFull <- readRDS(file = paste0(netDir, "allGOResultsBPMFCC.rds"))
}

fullGOAll <- list()
modules <- read.csv(paste0(netDir, "modules.csv"))
allColors <- modules$Color[c(2:nrow(modules), 1)]
for (color in allColors) {
  # add gene IDs and names based on "/"-separated ensemblIDs 
  for (ont in c("BP", "MF", "CC")) {
    curRes <- GOResultsFull[[ont]][[color]]@result %>%
      mutate(geneNames = lapply(lapply(str_split(geneID, "/"), IDToNameFun),toString))
    fwrite(curRes, paste0(netDir, "GO/fullCSV/", color, ont, "Res.csv"))
    fullGOAll[[ont]][[color]] <- curRes
  }
}

# ---- save for shiny ----
saveRDS(fullGOAll, file = paste0(shinyDir, "GOWGCNA.rds"))
saveRDS(fullPathway, file = paste0(shinyDir, "pathwayWGCNA.rds"))
# maintain compatability with old fullGO functionality
fullGO <- fullGOAll$BP


barplot(timeSplitWP[[day]][[reg]], showCategory=15)

# provides lookup across all timepoints and regulation levels and both KEGG
# and WP. Searches for the given term in the term descriptions
pathwayMethods <- c("KEGG")
getPathwayRes <- function(term) {
  infoCols <- c("method", "color", "sig")
  inPathways <- data.table(matrix(ncol = 16, nrow = 0))
  colnames(inPathways) <- c(infoCols, "category", "subcategory", "ID",
                            "Description", "GeneRatio", "BgRatio", "pvalue",
                            "p.adjust", "qvalue", "geneID", "Count","geneNames",
                            "ensembIDs")
  for (method in pathwayMethods) {
    for (color in allColors) {
      curRes <- read.csv(paste0(netDir, "pathways/", method ,"CSV/",
                                      color, method, "Res.csv"))
      
      matches <- curRes[grepl(term, curRes$Description, ignore.case = TRUE), ] %>% 
        mutate(method = method, color = color, sig = p.adjust < 0.05) %>%
        dplyr::select(all_of(infoCols), everything())
      inPathways <- rbind(inPathways, matches, fill = TRUE)
    }
  }
  print(inPathways[inPathways$sig, c(1:2, 7, 11)])
  print(paste(sum(inPathways$sig), "out of", nrow(inPathways)))
  return(inPathways[inPathways$sig, c(1:2, 7, 11)])
}

getGORes <- function(term) {
  infoCols <- c("ont", "color", "sig")
  inGO <- data.table(matrix(ncol = 13, nrow = 0))
  colnames(inGO) <- c(infoCols, "ID",	"Description",	"GeneRatio",	"BgRatio",	"pvalue",	
                            "p.adjust",	"qvalue",	"geneID",	"Count",	"geneNames")
  for (ont in c("BP", "CC", "MF")) {
      for (color in allColors) {
        curRes <- read.csv(paste0(netDir, "GO/fullCSV/", color, ont, "Res.csv"))
          matches <- curRes[grepl(term, curRes$Description, ignore.case = TRUE), ] %>% 
            mutate(ont = ont, color = color, sig = p.adjust < 0.05) %>%
            dplyr::select(all_of(infoCols), everything())
          inGO <- rbind(inGO, matches, fill = TRUE)
      }
  }
  print(inGO[inGO$sig, c(1:2, 5, 9)])
  print(paste(sum(inGO$sig), "out of", nrow(inGO)))
  return(inGO[inGO$sig, c(1:2, 5, 9)])
}

inPathwaysWGCNA3 <- list()
test <- "Wnt"

inPathways[[test]] <- getPathwayRes(test)
inPathwaysWGCNA3[[test]] <- getPathwayRes(test)
inPathways$DNA <- getPathwayRes("DNA")
inPathways$ASD <- getPathwayRes("ASD")

inGO <- list()
test <- "BMP"
test <- "Wnt"
test <- "hedgehog"
test <- "SHH" # none
test <- "autism"
test <- "ptch"

test <- "GABA" # green, purple, greenyellow
test <- "cell adhesion"

test <- "dorsal" # greenyellow lightcyan
test <- "synapse" # greenyellow lightcyan

test <- "frizzled" # 1 in lightcyan
test <- "smoothen" # none
test <- "patterning"
test <- "ventral"
test <- "anterior" # none
test <- "posterior" # none
test <- "motor" # none

inGO[[test]] <- getGORes(test)
inPathwaysWGCNA3[[test]] <- getPathwayRes(test)

fwrite(inGO[[test]], paste0(netDir, "wnt/inGOWGCNA2.csv"))
fwrite(inPathways[[test]], paste0(netDir, "wnt/inKEGGWGCNA2.csv"))

thesisTesting <- list()
thesisTesting$GO <- inGO
thesisTesting$Pathway <- inPathwaysWGCNA3

saveRDS(thesisTesting, paste0(netDir, "WGCNAThesisEnrichmentTest.rds"))
thesisTesting <- readRDS(paste0(netDir, "WGCNAThesisEnrichmentTest.rds"))
#saveRDS(fullPathway, file = paste0(shinyDir, "pathwayWGCNA.rds"))

test <- "migration"
test <- "gamma-aminobutyric acid"
test <- "ontogenesis"
test <- "ossification"
test <- "acetyltransferase"
test <- "histone"
test <- "transcription"

thesisTesting$GO[[test]] <- getGORes(test)
thesisTesting$Pathway[[test]] <- getPathwayRes(test)



fullWP <- readRDS(paste0(netDir, "pathways/fullWPWGCNA.rds"))
fullKEGG <- readRDS(paste0(netDir, "pathways/fullKEGGWGCNA.rds"))
GOResults <- readRDS(file =  paste0(netDir, "allGOResults.rds"))

pdf(paste0(netDir, "presentationFig/enrichmentPlots2.pdf"), height = 6, width = 7)
# mod 3, 4, 13
left <- ggarrange(cleanGOBarplot(GOResults[[allColors[3]]], title = "Module 3", hjust = 0),
                  labels = "AUTO", font.label = list(size = 18))
right <- ggarrange(cleanGOBarplot(GOResults[[allColors[4]]], title = "Module 4", hjust = 0),
                   cleanGOBarplot(GOResults[[allColors[13]]], title = "Module 13", hjust = 0),
                   nrow = 2, ncol = 1, labels = c("B", "C"), font.label = list(size = 18))
p <- ggarrange(left, right, nrow = 1, ncol = 2)
annotate_figure(p, top = text_grob(paste0("GO BP enrichment"), size = 14))

ggsave(paste0(netDir, "thesisFig/enrichmentPlotsM3-4-13.png"), height = 6, width = 7,
       dpi = 300, device = "png", bg = "white")
dev.off()

# mod 1, 2
p <- ggarrange(cleanGOBarplot(GOResults[[allColors[1]]],showCategory = 10, title = "Module 1", hjust = 1),
               cleanGOBarplot(GOResults[[allColors[2]]],showCategory = 10, title = "Module 2", hjust = 1),
               nrow = 1, ncol = 2, labels = "AUTO", font.label = list(size = 14))
annotate_figure(p, top = text_grob(paste0("GO BP enrichment"), size = 16))
ggsave(paste0(netDir, "thesisFig/enrichmentPlotsM1-2.png"), height = 6, width = 7,
       dpi = 300, device = "png", bg = "white")
dev.off()


# mod 8, 10
p <- ggarrange(cleanGOBarplot(GOResults[[allColors[8]]],showCategory = 10, title = "Module 8", hjust = 1),
               cleanGOBarplot(GOResults[[allColors[10]]],showCategory = 10, title = "Module 10", hjust = 1),
               nrow = 1, ncol = 2, labels = "AUTO", font.label = list(size = 14))
annotate_figure(p, top = text_grob(paste0("GO BP enrichment"), size = 16))
ggsave(paste0(netDir, "thesisFig/enrichmentPlotsM8-10.png"), height = 6, width = 7,
       dpi = 300, device = "png", bg = "white")
dev.off()


# height 4-9 or 3-9
pdf(paste0(netDir, "presentationFig/enrichmentPlots2.pdf"), height = 5, width = 9)

color <- "turquoise"
p <- ggarrange(cleanGOBarplot(GOResults[[color]], 10 , title = "GO"),
               cleanGOBarplot(fullKEGG[[color]], title = "KEGG"),
               cleanGOBarplot(fullWP[[color]], title = "WikiPathways"),
               nrow = 1, ncol = 3, labels = "AUTO")
annotate_figure(p, top = text_grob(paste0("M1 - ", color)))


color <- "purple"
p <- ggarrange(cleanGOBarplot(GOResults[[color]], 10 , title = "GO"),
               cleanGOBarplot(fullKEGG[[color]], title = "KEGG"),
               cleanGOBarplot(fullWP[[color]], title = "WikiPathways"),
               nrow = 1, ncol = 3, labels = "AUTO")
annotate_figure(p, top = text_grob(paste0("M2 - ", color)))

color <- "greenyellow"
p <- ggarrange(cleanGOBarplot(GOResults[[color]], 10 , title = "GO"),
               nrow = 1, ncol = 3, labels = "AUTO")
annotate_figure(p, top = text_grob(paste0("M3 - ", color)))

color <- "tan"
p <- ggarrange(cleanGOBarplot(GOResults[[color]], 10 , title = "GO", textWidth = 35),
               nrow = 1, ncol = 3, labels = "AUTO")
annotate_figure(p, top = text_grob(paste0("M4 - ", color)))



color <- "lightcyan"
p <- ggarrange(cleanGOBarplot(GOResults[[color]], 10, title = "GO"),
               cleanGOBarplot(fullKEGG[[color]], title = "KEGG"),
               cleanGOBarplot(fullWP[[color]], title = "WikiPathways"),
               nrow = 1, labels = "AUTO")
annotate_figure(p, top = text_grob(paste0("M8 - ", color)))

color <- "blue"
p <- ggarrange(cleanGOBarplot(GOResults[[color]], 15, title = "GO"),
               nrow = 1, ncol = 3, labels = "AUTO")
annotate_figure(p, top = text_grob(paste0("M10 - ", color)))

color <- "green"
p <- ggarrange(cleanGOBarplot(GOResults[[color]], title = "GO"),
               cleanGOBarplot(fullKEGG[[color]], title = "KEGG"),
               cleanGOBarplot(fullWP[[color]], title = "WikiPathways"),
               nrow = 1, ncol = 3, labels = "AUTO")
annotate_figure(p, top = text_grob(paste0("M13 - ", color)))
dev.off()
print(cleanGOBarplot(fullWP$turquoise))


pdf(file = paste0(netDir, "blackModExploratoryCnetplots.pdf"))
cnetplot(pathwayRes$KEGG[[color]], showCategory = 2, 
         cex.params = list(gene_label = 0.5))
#cnetplot(GOResultsNames[[color]], showCategory = "dorsal/ventral pattern formation")
dev.off()


print(cleanGOBarplot(GOResults$turquoise))

GO

# ---- network plotting ----
library(RCy3)
library(igraph)
library(tcltk)

modules <- read.csv(paste0(netDir, "modules.csv"))
kme <- read.csv(paste0(netDir, "kME.csv"))
# get modules with significant genotype results
sig <- read.csv(paste0(netDir, "padjustedMEs.csv")) %>% mutate(ME = X)
rownames(sig) <- sig$ME
# only necessary if we want a marker of genotype significance
genotypeSig <- sig %>%
  filter(Genotype < 0.05 | GenotypeDayInteraction < 0.05) %>%
  .$ME %>% gsub("E", "", .)
adjMatrix <- adjacency(WGCNA_TPM, type = "signed", power = 8)
#saveRDS(adjMatrix, paste0(netDir, "adjacencyMatrix"))
top50 <- list()
expr <- list()
curAdj <- list()
curAdjThin <- list()
curAdjFill <- list()
graphs <- list()
GOGraphs <- list()
curAdjGO <- list()
curAdjGOThin <- list()
DERes <- read.csv(paste0(resultDir, "DE/timepoints/mergedDESepTimepointsWTvsKO.csv"))
nameToID <- names(IDToName)
names(nameToID) <- IDToName

# NSUN5P2 was updated in Aug 2024, and both "ENSG00000106133" and "ENSG00000290831"
# are present as IDs to switch to. We select ENSG00000290831 as the default
nameToID["NSUN5P2"] <- "ENSG00000290831"

duplicateGeneNames <- IDToName[IDToName %in% IDToName[duplicated(IDToName)]] 


for (curME in unique(kme$moduleLabel)) {
  color <- modules$Color[gsub("E", "", modules$Label) == curME]
  curKMEs <- kme[kme$moduleLabel == curME,
                 colnames(kme) %in% c("gene_id", "gene_name", "moduleColor", "moduleLabel", gsub("E", "", curME))]
  colnames(curKMEs) <- c("gene_id", "gene_name", "moduleColor", "moduleLabel", "ME")
  top50[[curME]] <- curKMEs[order(curKMEs$ME, decreasing = TRUE), ][c(1:50), ]
  
  expr[[curME]] <- WGCNA_TPM[, top50[[curME]]$gene_id]
  IDs <- colnames(expr[[curME]])
  curAdj[[curME]] <- adjMatrix[IDs, IDs]
  colnames(curAdj[[curME]]) <- IDToName[IDs]
  rownames(curAdj[[curME]]) <- IDToName[IDs]
  # remove node to same node connections
  diag(curAdj[[curME]]) <- 0
  curAdjThin[[curME]] <- curAdj[[curME]]
  stackAdj <- stack(curAdj[[curME]])
  # to get top 50 bidirectional connections - get top 100 unidirectional connections
  lowest <- stackAdj[order(stackAdj$value, decreasing = TRUE)[100], ]$value
  curAdjThin[[curME]][curAdjThin[[curME]] < lowest] <- rep(0)
  
  # smaller graph with only connected nodes
  curAdjFill[[curME]] <- curAdjThin[[curME]][rowSums(curAdjThin[[curME]]) > 0, colSums(curAdjThin[[curME]]) > 0]
  graphs[[curME]] <- graph_from_adjacency_matrix(curAdjFill[[curME]], weighted = TRUE, diag = FALSE, mode = "undirected")
  g <- graphs[[curME]]
  
  # set degree, name length, and kME attributes
  degAll <- igraph::degree(g, v = igraph::V(g), mode = "all")
  g <- igraph::set_vertex_attr(g, "degree", index = igraph::V(g), value = degAll)
  g <- igraph::set_vertex_attr(g, "nameLen", index = igraph::V(g),
                               value = str_length(igraph::vertex_attr(g, "name")))
  g <- igraph::set_vertex_attr(g, "kME", index = igraph::V(g),
                               value = top50[[curME]]$ME[match(igraph::vertex_attr(g, "name"), 
                                                               top50[[curME]]$gene_name)])
  
  # get the differential expression results at D60 for colour coding
  genes <-  igraph::vertex_attr(g, "name")
  DEAtD60 <- rep("notDE", length(genes))
  names(DEAtD60) <- genes
  sigMatch <- DERes[match(genes, DERes$geneNames),] %>% filter(D60.padj < 0.05)
  DEAtD60[sigMatch$D60.log2FC > 0] <- "upreg"
  DEAtD60[sigMatch$D60.log2FC < 0] <- "downreg"
  g <- igraph::set_vertex_attr(g, "DEAtD60", index = igraph::V(g),
                              value = DEAtD60[genes])
  
  # graph with all genes in the GO results
  curGO <- read.csv(paste0(netDir, "GO/fullCSV/", color, "BPRes.csv")) %>%
    filter(p.adjust < 0.05)
  allGOGenes <- unique(unlist(str_split(curGO$geneNames, ", ")))
  curAdjGO[[curME]] <- adjMatrix[nameToID[allGOGenes], nameToID[allGOGenes]]
  
  colnames(curAdjGO[[curME]]) <- allGOGenes
  rownames(curAdjGO[[curME]]) <- allGOGenes
  # remove node to same node connections
  diag(curAdjGO[[curME]]) <- 0
  # cut off for 5% of connections
  cutOff <- sort(curAdjGO[[curME]])[round(length(curAdjGO[[curME]])*95/100)]
  
  curAdjGOThin[[curME]] <- curAdjGO[[curME]]
  curAdjGOThin[[curME]][curAdjGO[[curME]] < cutOff] <- 0
  GOG <- graph_from_adjacency_matrix(curAdjGOThin[[curME]], weighted = TRUE, diag = FALSE, mode = "undirected")
  # set degree, name length, and kME attributes
  
  degAllGO <- igraph::degree(GOG, v = igraph::V(GOG), mode = "all")
  GOG <- igraph::set_vertex_attr(GOG, "degree", index = igraph::V(GOG), value = degAllGO)
  GOG <- igraph::set_vertex_attr(GOG, "nameLen", index = igraph::V(GOG),
                               value = str_length(igraph::vertex_attr(GOG, "name")))
  GOGraphs[[curME]] <- GOG
  
  
  # get GO results of current module and select genes matching Wnt signalling
  curGORes <- read.csv(paste0(netDir, "GO/fullCSV/", color, "BPRes.csv"))
  curWntGenes <- str_split(curGORes$geneNames["Wnt signaling pathway" == curGORes$Description], ", ")
  # str split returns a list of length 1 if a single string is given as input
  if (length(curWntGenes) == 1) {
    curWntGenes <- curWntGenes[[1]]
    wntNames <- igraph::vertex_attr(g, "name")[igraph::vertex_attr(g, "name") %in% M8WNT]
    g <- igraph::set_vertex_attr(g, "inWNT", index = igraph::V(g),
                                 value = igraph::vertex_attr(g, "name") %in% M8WNT)
  } else {
    g <- igraph::set_vertex_attr(g, "inWNT", index = igraph::V(g),
                                 value = FALSE)
  }
  
  
  
  graphs[[curME]] <- g
}

createNetworkFromIgraph(GOGraphs$M8, title = paste(curME, "GO terms"))

createNetworkFromIgraph(g, title = paste(curME, "test"))
defaults <- list(NODE_SHAPE = "ellipse",
                 NODE_HEIGHT=30,
                 EDGE_TRANSPARENCY=120,
                 NODE_LABEL_POSITION="0.00,0.00")
styleName <- "WGCNATestG"
nodeWidths <- mapVisualProperty('node width','nameLen', 'c',  c(1,25), c(35,250))
nodeLabels <- mapVisualProperty('node label','id','p')
nodeFills <- mapVisualProperty('node fill color','DEAtD60','d',
                               c("notDE", "upreg", "downreg"), 
                               c('#F5EDDD', "darkgreen", "indianred"))


nodeFills <- mapVisualProperty('node fill color','degree','c',
                               c(0, 2, 12), 
                               c('#F5EDDD', "#F59777", "#F55333"))

createVisualStyle(style.name = styleName, defaults, list(nodeWidths, nodeLabels,nodeFills))
lockNodeDimensions(style.name = styleName, FALSE)

if (!dir.exists(paste0(netDir, "Cytoscape/"))) { dir.create(paste0(netDir, "Cytoscape/"))}
for (curME in names(graphs)) {
  g <- graphs[[curME]]
  
  createNetworkFromIgraph(g, title = paste(curME, "test"))
  setVisualStyle(styleName)
  wntNames <- vertex_attr(g, "name")[vertex_attr(g, "inWNT")]
  if (length(wntNames) != 0) {
    setNodeColorBypass(wntNames, "lightblue")
  }
  exportSVG(paste0(netDir, "Cytoscape/test", curME))
  exportPNG(paste0(netDir, "Cytoscape/test", curME))
}


exportVisualStyles(paste0(netDir, "visualStyles2"))


# ---- create a heatmap of the top go terms for each module ----
modules <- read.csv(paste0(netDir, "modules.csv"))
if (!exists("GORes")) {GORes <- readRDS(paste0(netDir, "allGOResults.rds"))}
enrichment <- list()
enrMatrix <- list()

# get similarity between two terms in go or pathway enrichment results 
# similarity calculated by Jaccard method - intersect over union
enrichmentJaccard <- function(termA, termB, curEnr) {
  genesA <- curEnr$geneNames[curEnr$Description == termA] %>%
    str_split_1(", ")
  genesB <- curEnr$geneNames[curEnr$Description == termB] %>%
    str_split_1(", ")
  length(intersect(genesA, genesB))/length(union(genesA, genesB))
}

if (!dir.exists(paste0(netDir, "enrichmentMatrices/"))) { dir.create(paste0(netDir, "enrichmentMatrices/"))}
for (mod in gsub("E", "", modules$Label[c(2:18, 1)])) {
  color <- modules$Color[gsub("E", "", modules$Label) == mod]
  enrichment$GO[[mod]] <- read.csv(paste0(netDir, "GO/fullCSV/", color, "BPRes.csv")) %>%
    filter(p.adjust < 0.05)
  enrichment$KEGG[[mod]] <- read.csv(paste0(netDir, "pathways/KEGGCSV/", color, "KEGGRes.csv")) %>%
    filter(p.adjust < 0.05)

  for (method in c("GO", "KEGG")) {
    curEnr <- enrichment[[method]][[mod]]
    if (nrow(curEnr) > 1) {
      curMatrix <- data.frame(matrix(nrow = nrow(curEnr),
                                     ncol = nrow(curEnr)))
      rownames(curMatrix) <- curEnr$Description
      colnames(curMatrix) <- curEnr$Description
      
      
      for (i in c(1:length(curEnr$Description))) {
        for (j in c(i:length(curEnr$Description))) {
          termA <- curEnr$Description[i]
          termB <- curEnr$Description[j]
          curMatrix[termA, termB] <- enrichmentJaccard(termB, termA, curEnr)
          curMatrix[termB, termA] <- curMatrix[termA, termB]
        }
      }
      curMatrix <- as.matrix(curMatrix)
      enrMatrix[[method]][[mod]] <- curMatrix
      fwrite(as.data.table(curMatrix, keep.rownames = TRUE), paste0(netDir, "enrichmentMatrices/", method, mod, ".csv"))
      
    }
  }
}



# mod 15 and 17 have 0 GO terms, 12 only has 1 GO term
modToSize <- data.frame(mod =         c('M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9', 'M10', 'M11', 'M13', 'M14', 'M16', 'M0'), 
                        heatmapSize = c(0.5,   0.6, 6,     4,    7,    6,    4,    3,    8,    8,     5,     8,     5,     5,     6))
pdf(paste0(netDir, "heatmaps/GOHeatmap.pdf"))
for (mod in modToSize$mod) {
  size <- modToSize$heatmapSize[modToSize$mod == mod]
  if (!is.null(enrMatrix$GO[[mod]])) {
    print(pheatmap(enrMatrix$GO[[mod]], cellwidth = size, cellheight = size, 
                 fontsize = size, legend = FALSE, treeheight_row = 0, main = mod))
  }
}
dev.off()  

size <- 5
pdf(paste0(netDir, "heatmaps/KEGGHeatmap.pdf"))
for (mod in names(enrMatrix$KEGG)) {
  print(pheatmap(enrMatrix$KEGG[[mod]], cellwidth = size, cellheight = size, 
                   fontsize = size, treeheight_row = 0, main = mod))
}
dev.off()  

pdf(paste0(netDir, "heatmaps/WPHeatmap.pdf"))
for (mod in names(enrMatrix$WP)) {
  print(pheatmap(enrMatrix$WP[[mod]], cellwidth = size, cellheight = size, 
                 fontsize = size, treeheight_row = 0, main = mod))
}
dev.off()  

# annotate graphs
# mod 8 has a wnt cluster
mod <- 'M8'
wntTerms <- c("negative regulation of Wnt signaling pathway",
              "regulation of canonical Wnt signaling pathway",
              "regulation of Wnt signaling pathway",
              "canonical Wnt signaling pathway",
              "Wnt signaling pathway",
              "cell−cell signaling by wnt")
curGO <- enrichment$GO[[mod]]
genes <- curGO$geneNames[curGO$Description %in% wntTerms] %>%
  str_split(", ") %>%
  unlist() %>%
  unique()
n <- createNetworkFromIgraph(GOGraphs[[mod]], title = paste(mod, "GO terms"))
defaults <- list(NODE_SHAPE = "ellipse",
                 NODE_HEIGHT=30,
                 EDGE_TRANSPARENCY=120,
                 NODE_LABEL_POSITION="0.00,0.00",
                 NODE_FILL_COLOR = "#F5EDDD")
styleName <- "WGCNATestGO4"
nodeWidths <- mapVisualProperty('node width','nameLen', 'c',  c(1,25), c(35,250))
nodeLabels <- mapVisualProperty('node label','id','p')
createVisualStyle(style.name = styleName, defaults, list(nodeWidths, nodeLabels))
lockNodeDimensions(style.name = styleName, FALSE)

setVisualStyle(styleName)
setNodeColorBypass(genes, "lightgreen", network = n)


mod <- 'M1'
# mod 1 has a replication cluster
repTerms <- c("chromosome segregation",
              "DNA replication",
              "DNA-templated DNA replication",
              "sister chromatid segregation",
              "nuclear chromosome segregation",
              "mitotic nuclear division",
              "mitotic sister chromatid segregation",
              "protein-DNA complex assembly",
              "nuclear division")
curGO <- enrichment$GO[[mod]]
genes <- curGO$geneNames[curGO$Description %in% repTerms] %>%
    str_split(", ") %>%
    unlist() %>%
    unique()
M1Net <- createNetworkFromIgraph(GOGraphs[[mod]], title = paste(mod, "GO terms"))
setVisualStyle(styleName)
setNodeColorBypass(genes, "lightblue", network = M1Net)

mod <- 'M3'
# mod 3 has a cluster of differentiation terms, 
# especially GABA and forebrain neurons
foreTerms <- c("forebrain generation of neurons",
               "forebrain neuron differentiation",
               "cerebral cortex neuron differentiation",
               "cerebral cortex GABAergic interneuron differentiation",
               "GABAergic neuron differentiation")
curGO <- enrichment$GO[[mod]]
genes <- curGO$geneNames[curGO$Description %in% foreTerms] %>%
    str_split(", ") %>%
    unlist() %>%
    unique()
M3Net <- createNetworkFromIgraph(GOGraphs[[mod]], title = paste(mod, "GO terms"))
setVisualStyle(styleName)
setNodeColorBypass(genes, "khaki", network = M3Net)

mod <- 'M13'
# mod 13 only has two GO terms - focus on GABA term
gabaTerms <- c("gamma-aminobutyric acid signaling pathway")
curGO <- enrichment$GO[[mod]]
genes <- curGO$geneNames[curGO$Description %in% gabaTerms] %>%
    str_split(", ") %>%
    unlist() %>%
    unique()
M13Net <- createNetworkFromIgraph(GOGraphs[[mod]], title = paste(mod, "GO terms"))
setVisualStyle(styleName)
setNodeColorBypass(genes, "pink", network = M13Net)



# plotting Go graphs
if (!exists("GORes")) {GORes <- readRDS(paste0(netDir, "allGOResults.rds"))}
GO <- GORes$lightcyan %>% filter(p.adjust < 0.05)
plotGOgraph(GO, firstSigNodes = 3)


# ---- fisher test of direct ZMYND8 binding targets to module membership ----

ChIP <- read.csv("Z:/mnt/Data1/PROJECTS/ZMYND8/IV/data.sig.chip_Exp.csv")
bindingGenes <- ChIP$EnsID

allMod <- gsub("E", "", modules$Label)
bindingCor <- data.frame(matrix(nrow = length(allMod), ncol = 3))
rownames(bindingCor) <- allMod
colnames(bindingCor) <- c("oddsRatio", "p-value", "p-adjusted")

for (m in allMod) {
  cor <- fisher.test(kme$moduleLabel == m, kme$gene_id %in% bindingGenes)
  bindingCor[m, "oddsRatio"] <- cor$estimate
  bindingCor[m, "p-value"] <- cor$p.value
}

# apply BH pval adjustment
for (pvalCol in grep("p-value", colnames(bindingCor))) {
  # check the column selected is the padjusted one
  stopifnot(grepl("p-adjusted", colnames(bindingCor[, pvalCol + 1])))
  bindingCor[, pvalCol + 1] <- p.adjust(bindingCor[, pvalCol], method = "BH")
}

fwrite(bindingCor, paste0(netDir, "correlationToZMYND8DirectTargets.csv"), row.names = TRUE)

# sanity check - filter KME for only the binding genes, and check that lots of them 
# seem to be in mod 7

bindingKMEs <- kme[kme$gene_id %in% bindingGenes,]
table(bindingKMEs$moduleLabel)
# yes - visually seems enriched in M7
bindingKMEs$gene_name[bindingKMEs$moduleLabel == "M7"]


#' create supplementary tables:
#' For each of the 7 modules discussed: 1,2,3,4,8,10,13
#'  Output the significant results for all three ontologies 
#'  Merge with pathway enrichment

colorToME <- modules$Label %>% gsub("E", "", .)
names(colorToME) <- modules$Color
METoColor <- names(colorToME)
names(METoColor) <- colorToME

IDToNameFun <- function(IDs) {
  unname(IDToName[IDs])
}

#WGCNAMerged <- list()
for (mod in c("M1","M2","M3","M4","M8","M10","M13")) {
    curData <- data.frame(matrix(nrow = 0, ncol = 12))
    colnames(curData) <- c("enrichmentType", "enrichmentSuptype", "ID",	"Description",	"GeneRatio",	"BgRatio",
                           "pvalue",	"p.adjust",	"qvalue",	"geneID",	"Count", "geneNames")
    
    for (ont in c("BP", "MF", "CC")) {
      # add geneName data
      curCSV <- read.csv(paste0(netDir, "GO/fullCSV/", METoColor[mod], ont, "Res.csv")) %>%
        filter(p.adjust < 0.05) 
      if (nrow(curCSV) > 0) {
        curCSV$enrichmentType <- "GO"
        curCSV$enrichmentSuptype <- ont
        curCSV <- curCSV %>% dplyr::select(enrichmentType, enrichmentSuptype, everything())
        curData <- rbind(curData, curCSV)
      }
    }
    curKEGG <- read.csv(paste0(netDir, "pathways/KEGGCSV/", METoColor[mod], "KEGG", "Res.csv")) %>% filter(p.adjust < 0.05)
    if (nrow(curKEGG) > 0) {
      curKEGG$enrichmentType <- "pathways"
      curKEGG$enrichmentSuptype <- "KEGG"
      curKEGG <- curKEGG %>%
        dplyr::select(!c(geneID, subcategory, category)) %>%
        mutate(geneID = ensembIDs) %>%
        dplyr::select(enrichmentType, enrichmentSuptype, !ensembIDs) 
      curData <- rbind(curData, curKEGG)
      curData <- curData %>% .[order(.$p.adjust),]
    }
    fwrite(curData, paste0(resultDir, "thesisSup/", mod, "MergedRes.csv"))
}


# Excel workbook code based off:
# https://www.r-bloggers.com/2019/08/creating-excel-workbooks-with-multiple-sheets-in-r/

# Get the file name read in as a column
read_filename <- function(fname) {
  read_csv(fname, col_names = TRUE) %>%
    mutate(filename = fname)
}
tbl <-
  list.files(path = paste0(resultDir, "thesisSup/"),
             pattern ="M*MergedRes.csv",
             full.names = TRUE) %>%
  map_df(~read_filename(.))
tbl$filename <- tbl$filename %>% gsub(".*/", "", .) %>% gsub(".csv", "", .)

mylist <- tbl %>% split(.$filename)
names(mylist)

wb <- createWorkbook()
lapply(seq_along(mylist), function(i){
  addWorksheet(wb=wb, sheetName = names(mylist[i]))
  writeData(wb, sheet = i, mylist[[i]][-length(mylist[[i]])])
})
# add merged DE results

saveWorkbook(wb, paste0(resultDir, "thesisSup/WGCNAResults.xlsx"), overwrite = TRUE)

