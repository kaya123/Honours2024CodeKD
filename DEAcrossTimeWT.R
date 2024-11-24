# Overview from https://docs.google.com/document/d/1AJ1IW-3AZ9_yB2gNR8Ri7irl-GU-JHTVbN1nZKAPoV8/edit
# Differential expression analysis
# Effect of timepoint within WT samples, correcting for batch
# Subset for WT samples, then:
#  ~ batch +Day


# Proposed program structure:
# 1. Take in preprocessed data from "Z:/mnt/Data1/PROJECTS/ZMYND8/Organoid_BulkRNA-seq/Analysis_Results_Kaya/preprocessing/preprocessedData.rds")
# 2. Run DESeq2 ~batch+Day
# 3. Contrast days

# ---- Step 0: Set up and package management ---
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2", suppressUpdates = TRUE) 
BiocManager::install("clusterProfiler") # GO enrichment
BiocManager::install("org.Hs.eg.db") # GO annotation

library(ggplot2)
library(dplyr)
library(tidyverse)
library(DESeq2)
library(patchwork)
library(data.table)
library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(stringr)
library(gridExtra)
library(grid)
library(cowplot)

library(ggpubr)

# constants:
resultDir <- "Z:/mnt/Data1/PROJECTS/ZMYND8/Organoid_BulkRNA-seq/Analysis_Results_Kaya/"
dayPairs <- c("D20vsD60", "D20vsD100", "D60vsD100")

# ---- 1. Take in preprocessed data -----
preprocessedData <- readRDS(paste0(resultDir, "preprocessing/preprocessedData.rds"))
counts <- preprocessedData$filteredCounts
meta <- preprocessedData$metadata

# ---- 2. Subset for WT -----
metaWT <- as_tibble(meta[meta$Genotype == "WT", ])
countsWT <- round(counts[, meta$Genotype == "WT"])

# ---- 3. Run DESeq2 ~batch+Day ----
rownames(metaWT) <- metaWT$Sample
model.matrix(~ Batch + Day, data = metaWT)
WT_DE <- DESeqDataSetFromMatrix(countData = countsWT, 
                                  colData = metaWT,
                                  design = ~ Batch + Day)

WT_DE$Day <- relevel(WT_DE$Day, ref = "D20")
WT_DEFull <- DESeq(WT_DE)

WT_DERes <- list()
WT_DERes$D20vsD60 <- results(WT_DEFull, contrast = c("Day", "D60", "D20"))
WT_DERes$D20vsD100 <- results(WT_DEFull, contrast = c("Day", "D100", "D20"))
WT_DERes$D60vsD100 <- results(WT_DEFull, contrast = c("Day", "D100", "D60"))

# ---- 4. Save DE genes ----
saveRDS(WT_DERes, file = paste0(resultDir, "DE/WT/WT_DEResults.rds"))
WTSigRes <- list()
for (dayPair in dayPairs) {
  WTSigRes[[dayPair]] <- WT_DERes[[dayPair]] %>%
    as.data.frame() %>%
    .[order(.$padj), ] %>%
    filter(.$padj < 0.05 & abs(.$log2FoldChange) > log2(1.5)) %>%
    rownames_to_column(var = "geneID") %>%
    mutate(geneName = IDToName[geneID]) %>%
    select(geneID, geneName, everything())
  fwrite(WTSigRes[[dayPair]], file = paste0(resultDir, "DE/WT/WT", dayPair, "DEGenes.csv"))
}


# ---- 5. volcano plots ----
pdf(file = paste0(resultDir, "DE/WT/WTVolcano.pdf"))
vcPlots = list()
i <- 1
for (dayPair in dayPairs) {
    fc <- WT_DERes[[dayPair]]$log2FoldChange
    padj <- WT_DERes[[dayPair]]$padj
    fcPvalData <- data.frame(fc = fc, padj = padj, sig = padj < 0.05) %>%
      filter(!is.na(padj)) %>%
      filter(abs(fc) > log2(1.5))
    vcPlots[[i]] = ggplot(fcPvalData, aes(fc, -log2(padj), color = sig)) +
      geom_point(size = 0.5) +
      scale_color_manual(values = c("darkred", "black"), 
                         breaks = TRUE) +
      labs(x = "log2 fold change", y = "-log2 p-adjusted",
           color = "P-adjusted < 0.05", title = paste(dayPair)) +
      guides(color = guide_legend(override.aes = list(size = 3)))
    i <- i+1
}

# combine the three timpoints with a single legend
combined <- vcPlots[[1]] + vcPlots[[2]] + vcPlots[[3]] &
  theme(legend.position = "bottom")
volcanoWT <- combined + plot_layout(guides = "collect")
print(volcanoWT)
dev.off()

# ---- 8. Export differentially expressed genes, WT vs KO only, meeting p adjusted < 0.05, abs(log2FC)>log2(1.5) in any timepoint ----
sigWTDE <- list()
for (dayPair in dayPairs) {
  sigWTDE[[dayPair]] <- WT_DERes[[dayPair]] %>%
    as.data.frame() %>%
    .[order(.$padj), ] %>%
    filter(.$padj < 0.05 & abs(.$log2FoldChange) > log2(1.5)) %>%
    rownames_to_column(var = "geneID") %>%
    .$geneID
}
allSigWTDE <- unique(c(sigWTDE$D20vsD60, sigWTDE$D20vsD100, sigWTDE$D60vsD100))
mergedWTData <- data.frame(row.names = allSigWTDE, geneIDs = allSigWTDE,
                           geneNames = IDToName[allSigWTDE])


D20vsD60 <- WT_DERes$D20vsD60[allSigWTDE,]
D20vsD100 <- WT_DERes$D20vsD100[allSigWTDE,]
D60vsD100 <- WT_DERes$D60vsD100[allSigWTDE,]

mergedWTData <- mergedWTData %>%
  mutate(D20vsD60.log2FC = D20vsD60$log2FoldChange,
         D20vsD60.padj = D20vsD60$padj,
         D20vsD100.log2FC = D20vsD100$log2FoldChange,
         D20vsD100.padj = D20vsD100$padj,
         D60vsD100.log2FC = D60vsD100$log2FoldChange,
         D60vsD100.padj = D60vsD100$padj) %>%
  .[order(.$D20vsD60.padj),] %>%
  replace_na(list(D20vsD60.padj = 1, D20vsD100.padj = 1, D60vsD100.padj = 1))
fwrite(mergedWTData, paste0(resultDir, "DE/WT/mergedDEAcrossTimepointsWT.csv"))
saveRDS(WT_DERes, file = paste0(resultDir, "DE/WT/WT_DEResults.rds"))

# ---- 5. GO enrichment - split into up and down regulated ----
data(geneList, package = "DOSE")
WT_DEGeneNames <- list()
WT_GORes <- list()
WTSigRes <- list()
for (dayPair in dayPairs) {
  WTSig <- read.csv(paste0(resultDir, "DE/WT/WT", dayPair, "DEGenes.csv"))
  WT_DEGeneNames[[dayPair]] <- list()
  WT_GORes[[dayPair]] <- list()
  
  # upregulated 
  WT_DEGeneNames[[dayPair]]$up <- WTSig$geneID[WTSig$log2FoldChange > 0]
  WT_GORes[[dayPair]]$up <- enrichGO(WT_DEGeneNames[[dayPair]]$up, 
                                     universe = rownames(countsWT),
                             OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP",
                             pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                             minGSSize = 10, maxGSSize = 500)
  # downregulated 
  WT_DEGeneNames[[dayPair]]$down <- WTSig$geneID[WTSig$log2FoldChange < 0]
  WT_GORes[[dayPair]]$down <- enrichGO(WT_DEGeneNames[[dayPair]]$down, 
                                       universe = rownames(countsWT),
                                     OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP",
                                     pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                                     minGSSize = 10, maxGSSize = 500)
}
# Save GO top 15 enriched BP term graphs
pdf(file = paste0(resultDir, "GO/topBPTermsWTTimepointsSepUpDownReg.pdf"))
par(mfrow = c(1,1))
for (dayPair in dayPairs) {
  print(barplot(WT_GORes[[dayPair]]$up, showCategory = 15,
                title = paste("GO BP - WT upregulated ", dayPair)))
  print(barplot(WT_GORes[[dayPair]]$down, showCategory = 15,
                title = paste("GO BP - WT downregulated", dayPair)))
}
dev.off()

# Save full GO output data 
saveRDS(WT_GORes, file = paste0(resultDir, "GO/fullGOAcrossWT.rds"))
# if the user has run this code before and is rerunning some analysis,
# it is faster to load GO output from rds
if (!exists("WT_GORes")) {
  WT_GORes <- readRDS(paste0(resultDir, "GO/fullGOAcrossWT.rds"))
}

# Save go output as CSV
if (!dir.exists(paste0(resultDir, "GO/WTFullCSV/"))) {
  dir.create(paste0(resultDir, "GO/WTFullCSV/"))
}
for (dayPair in dayPairs) {
  for (dir in c("Up", "Down")) {
    fwrite(WT_GORes[[dayPair]][[str_to_lower(dir)]]@result,
           paste0(resultDir, "GO/WTFullCSV/", dayPair, dir, "GOBPRes.csv"))
  }
}

# ---- 6. Create combined figure 1 for practice thesis ----
plotsVC <- list()
plotsVCNoL <- list()
i <- 1
for (dayPair in dayPairs) {
  fc <- WT_DERes[[dayPair]]$log2FoldChange
  padj <- WT_DERes[[dayPair]]$padj
  fcPvalData <- data.frame(fc = fc, padj = padj, sig = padj < 0.05) %>%
    filter(!is.na(padj)) %>%
    filter(abs(fc) > log2(1.5))
  nUpreg <- sum(fcPvalData$fc > 0 & fcPvalData$sig)
  nDownreg <- sum(fcPvalData$fc < 0 & fcPvalData$sig)
  plotsVC[[i]] <- ggplot(fcPvalData, aes(fc, -log2(padj), color=sig)) +
    geom_point(size = 0.5) +
    scale_color_manual(values = c("darkred", "black"), 
                       breaks = TRUE) +
    labs(x = "log2 fold change", y = "-log2 P-adjusted", color = "P-adjusted < 0.05",
         title = paste(dayPair), 
         subtitle = paste0("Significantly up: ",
                           nUpreg, "\n                down: ", nDownreg)) +
    theme(legend.position = "bottom") +
    guides(color = guide_legend(override.aes = list(size = 3))) 
  plotsVCNoL[[i]] <- plotsVC[[i]] + theme(legend.position = "none")
  
  i <- i+1
}

GOplots <- list()
for (dayPair in dayPairs) {
  GOplots[[dayPair]] <- list()
  for (dir in c("up", "down")) {
    GOSubset <- WT_GORes[[dayPair]][[dir]]@result[1:10,] %>% filter(p.adjust < 0.05)
    breaks <- c(max(GOSubset$p.adjust), min(GOSubset$p.adjust))
    labels <- c(max(signif(GOSubset$p.adjust, 2)), min(signif(GOSubset$p.adjust, 2)))
    if (dayPair == "D60vsD100" & dir == "up") {
      breaks <- c(max(GOSubset$p.adjust))
      labels <- c(max(signif(GOSubset$p.adjust, 3)))
    } 
    GOplots[[dayPair]][[dir]] <- ggplot(GOSubset[order(GOSubset$p.adjust, decreasing = TRUE),],
                                        aes(x = -log10(p.adjust),
                                            y = fct_inorder(Description),
                                            fill = log10(p.adjust))) +
      geom_bar(stat = "identity") +
      scale_y_discrete(labels = function(x) str_wrap(str_replace_all(x, " " , "\n"),
                                                     width = 25)) +
      labs(y = "GO description", x = "-log10 P-adjusted         ",
           title = paste0("   ", dayPair, " ", dir, "reg")) + 
      theme(legend.position = "none",
           plot.title.position = "plot",
           axis.title.x = element_text(hjust = 1.5)
      )
    
  } 
}


top <- ggarrange(plotsVCNoL[[1]], plotsVCNoL[[2]], plotsVCNoL[[3]],
                 labels = c("A", "B", "C"), ncol = 3, nrow = 1,
                 common.legend = TRUE, legend = "bottom")

bottom <- ggarrange(GOplots$D20vsD60$up, GOplots$D20vsD60$down,
                    labels = c("D", "E"), ncol = 2)
pdf(file = paste0(resultDir, "DE/WT/WTVolcanoAndGO.pdf"), height = 10, width = 8)  
ggarrange(top, bottom, ncol = 1, heights = c(1,1))
dev.off()

ggarrange(top, bottom, ncol = 1, heights = c(1,1))
ggsave(paste0(resultDir, "DE/WT/WTVolcanoAndGO.png"), height = 10, width = 8,
       device = "png", dpi = 300, background = "white")


