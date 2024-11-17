# Overview from https://docs.google.com/document/d/1AJ1IW-3AZ9_yB2gNR8Ri7irl-GU-JHTVbN1nZKAPoV8/edit
# Differential expression analysis
# KO vs WT and Het vs WT at each timepoint controlling for batch
# Subset by Day
# ~ Batch + Genotype
# > contrasts to extract KO vs WT and Het vs WT
# Compare with scRNA-seq results at d20, d60, d100.


# Proposed program structure:
# 1. Take in preprocessed data from "Z:/mnt/Data1/PROJECTS/ZMYND8/Organoid_BulkRNA-seq/Analysis_Results_Kaya/preprocessing/preprocessedData.rds")
# 2. Pre-filtering - should this be performed in the preprocessing? MOVED TO PREPROCESSING
# 3. Subset data into the three timepoint categories
# 4. Run DESeq2 with ~ Batch + Genotype
# 5. Use contrast to find KO vs WT and Het vs WT (described here https://www.youtube.com/watch?v=X6p3E-QTcUc)
# 6. Export histograms of expression and untransformed and rlog transformed data 

# ---- Step 0: Set up and package management ---
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2", suppressUpdates=TRUE) 

# Plotting
library(ggplot2)
library(gplots) # for heatmap.2
library(RColorBrewer)
library(patchwork)
library(ggpubr)
# Data manipulation
library(dplyr)
library(tidyverse)
library(data.table)
library(stringr)
# DE and GO
library(DESeq2)
library(clusterProfiler)
library(enrichplot)
library(biomaRt)
# Save supplementary data
library(openxlsx)

resultDir <- "Z:/mnt/Data1/PROJECTS/ZMYND8/Organoid_BulkRNA-seq/Analysis_Results_Kaya/"
shinyDir <- "Z:/mnt/Data1/PROJECTS/ZMYND8/Organoid_BulkRNA-seq/Scripts_Kaya/shiny/pathwayEnrichment/"

timepointDir <- paste0(resultDir, "DE/timepoints/")
if (!dir.exists(timepointDir)) { dir.create(timepointDir) }

# ---- 1. Take in preprocessed data -----
preprocessedData <- readRDS(paste0(resultDir, "preprocessing/preprocessedData.rds"))
counts <- preprocessedData$filteredCounts
meta <- preprocessedData$metadata
IDToName <- preprocessedData$IDToName

# ---- 2. Subset data into the three timepoint categories ----
days <- c("D20", "D60", "D100")
timeSplit <- list()
timeSplitData <- list()

for (day in days) {
  timeSplit[[day]] <- meta[meta$Day == day, ]
  timeSplit[[day]]$Genotype <- factor(timeSplit[[day]]$Genotype)
  timeSplit[[day]]$Batch <- factor(timeSplit[[day]]$Batch)
  
  timeSplitData[[day]] <- round(counts[colnames(counts) %in% timeSplit[[day]]$Sample])
}

# ---- 3. Run DESeq2 with ~Batch + Genotype, using contrast to find KO vs WT and Het vs WT ----
timeSplitDE <- list()
timeSplitDEFull <- list()
timeSplitDERes <- list() 
timeSplitDErlog <- list()

for (day in days) {
  timeSplitDE[[day]] <- DESeqDataSetFromMatrix(countData = timeSplitData[[day]],
                                             timeSplit[[day]],
                                             design = ~ Batch + Genotype)
  timeSplitDE[[day]]$Genotype <- relevel(timeSplitDE[[day]]$Genotype, ref = "WT")
  timeSplitDEFull[[day]] <- DESeq(timeSplitDE[[day]])
  
  timeSplitDERes[[day]] <- list()
  timeSplitDERes[[day]]$WTvsKO <- results(timeSplitDEFull[[day]],
                                          contrast = c("Genotype", "KO", "WT"))
  timeSplitDERes[[day]]$WTvsHet <- results(timeSplitDEFull[[day]],
                                           contrast = c("Genotype", "Het", "WT"))
  
  timeSplitDErlog[[day]] <-rlog(timeSplitDEFull[[day]])
}


# ---- 4. Export MA plots of all timepoints WTvsHet and WTvsKO ----
pdf(file = paste0(timepointDir, "MAPlots.pdf"))
par(mfrow = c(3,2))
for (day in days) {
  plotMA(timeSplitDERes[[day]]$WTvsHet, ylim = c(-2,2), 
         main = paste("WT vs Het", day))
  plotMA(timeSplitDERes[[day]]$WTvsKO, ylim = c(-2,2),
         main = paste("WT vs KO", day))
  
}  
dev.off()

# ---- 5. Export histograms of expression of rlog transformed data ----
pdf(file = paste0(timepointDir, "normalisedExpressionDist.pdf"))
for (day in days) {
  longer <- assay(timeSplitDErlog[[day]]) %>% 
    as_tibble(rownames = "GeneID") %>%
    pivot_longer(-1, names_to = "Sample", values_to = "rlogExpression")
  
   print(ggplot(longer, aes(x = Sample, y = rlogExpression)) +
    geom_violin() +
    geom_point() +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(title = paste("Normalised expression", day)))
}
dev.off()

# ---- 6. Export differentially expressed genes in all 6 categories ----
saveRDS(timeSplitDERes, file = paste0(timepointDir, "timeSplitDEResults.rds"))
# if only one section of the code is being run, load previously saved DE data
if (!exists("timeSplitDERes")) {
  timeSplitDERes <- readRDS(paste0(resultDir, "DE/timepoints/timeSplitDEResults.rds"))
}
for (day in days) {
  for (group in c("WTvsHet", "WTvsKO")) {
    sigRes <- timeSplitDERes[[day]][[group]] %>%
      as.data.frame() %>%
      .[order(.$padj), ] %>%
      filter(.$padj < 0.05 & abs(.$log2FoldChange) > log2(1.5)) %>%
      rownames_to_column(var = "geneID") %>%
      mutate(geneName = IDToName[geneID]) %>%
      dplyr::select(geneID, geneName, everything())
    fwrite(sigRes, paste0(timepointDir, "DEGenes", day, group, ".csv"))
  }
}

# ---- 8. Export all differentially expressed genes, WT vs KO only, meeting p adjusted < 0.05, abs(log2FC)>log2(1.5) in any timepoint ----
sigDEGenes <- list()
for (day in days) {
  sigDEGenes[[day]] <- timeSplitDERes[[day]]$WTvsKO %>%
    as.data.frame() %>%
    .[order(.$padj), ] %>%
    filter(.$padj < 0.05 & abs(.$log2FoldChange) > log2(1.5)) %>%
    rownames_to_column(var = "geneID") %>%
    .$geneID
}
allSigDE <- unique(c(sigDEGenes$D20, sigDEGenes$D60, sigDEGenes$D100))
mergedDE <- data.frame(row.names = allSigDE, geneIDs = allSigDE,
                       geneNames = IDToName[allSigDE])

D20Sig <- timeSplitDERes$D20$WTvsKO[allSigDE, ]
D60Sig <- timeSplitDERes$D60$WTvsKO[allSigDE, ]
D100Sig <- timeSplitDERes$D100$WTvsKO[allSigDE, ]
mergedDE <- mergedDE %>% 
  mutate(D20.log2FC = D20Sig$log2FoldChange,
         D20.padj = D20Sig$padj,
         D60.log2FC = D60Sig$log2FoldChange,
         D60.padj = D60Sig$padj,
         D100.log2FC = D100Sig$log2FoldChange,
         D100.padj = D100Sig$padj) %>%
  .[order(.$D20.padj),] %>%
  replace_na(list(D20.padj = 1, D60.padj = 1, D100.padj = 1))
fwrite(mergedDE, paste0(timepointDir, "mergedDESepTimepointsWTvsKO.csv"))

# ---- 7. Export -log padj against fold change ----
pdf(file = paste0(timepointDir, "fcVsPadjSepTimepoints.pdf"))
plots <- list()
i <- 1
for (day in days) {
  for (group in c("WTvsHet", "WTvsKO")) {
    fc <- timeSplitDERes[[day]][[group]]$log2FoldChange
    padj <- timeSplitDERes[[day]][[group]]$padj
    fcPvalData <- data.frame(fc = fc, padj = padj, sig = padj < 0.05) %>%
      filter(!is.na(padj))
    plots[[i]] = ggplot(fcPvalData, aes(fc, -log2(padj), color = sig)) +
      geom_point() +
      scale_color_manual(values=c("black", "darkred")) +
      labs(x = "log2 fold change", y = "-log2 P-adj", color = "P-adj <0.05",
           title = paste(group, day))
    i <- i+1
  }
}
wrap_plots(plots, ncol = 2)
dev.off()

# ---- 8. Export heatmaps of most DE genes and their expression across WT & KO ----
pdf(file = paste0(timepointDir, "heatmapTop35Genes.pdf"))
for (day in days) {
  WTKOSamples <- meta$Sample[meta$Genotype != "Het" & meta$Day == day]
  topGenes <- timeSplitDERes[[day]]$WTvsKO %>%
    .[order(.$padj), ] %>%
    head(35) %>%
    rownames()
  topGenesData <- assay(timeSplitDErlog[[day]])[topGenes, WTKOSamples]
  colnames(topGenesData) <- gsub("_GOK.*","",colnames(topGenesData))
  rownames(topGenesData) <- IDToName[rownames(topGenesData)]
  
  heatmap.2(topGenesData, scale = "row",
            trace = "none", dendrogram = "column",
            col = colorRampPalette(rev(brewer.pal(8, "RdBu")))(255),
            margins = c(8,8), cexCol = 1, cexRow = 0.8,
            main = paste("           ", day, "35 most significantly DE genes"))
}
dev.off()

# ---- 9. GO - focus on only WT vs KO ----
data(geneList, package = "DOSE")
timeSplitGeneIDs <- list()
timeSplitGORes <- list()
group <- "WTvsKO"
days <- c("D20", "D60", "D100")
for (day in days) {
  upDownReg <- timeSplitDERes[[day]][[group]] %>%
    as.data.frame() %>% 
    filter(!is.na(padj) & padj < 0.05)
  timeSplitGeneIDs[[day]]$up <- upDownReg %>%
    filter(log2FoldChange > log2(1.5)) %>%
    rownames()
  timeSplitGeneIDs[[day]]$down <- upDownReg %>%
    filter(log2FoldChange < -log2(1.5)) %>%
    rownames()
}

for (day in days) {
  for (reg in c("up", "down")) {
    timeSplitGORes[[day]][[reg]] <- enrichGO(timeSplitGeneIDs[[day]][[reg]],
                                             universe = rownames(counts),
                                               OrgDb = "org.Hs.eg.db", 
                                               keyType = "ENSEMBL", ont = "BP",
                                               pvalueCutoff = 0.05, 
                                               qvalueCutoff = 0.05,
                                               minGSSize = 10, maxGSSize = 500)
  }
}
otherGO <- list()
for (day in days) {
    for (reg in c("up", "down")) {
        for (ont in c("MF", "CC", "BP")) {
            otherGO[[day]][[reg]][[ont]] <- enrichGO(timeSplitGeneIDs[[day]][[reg]],
                                                 universe = rownames(counts),
                                                 OrgDb = "org.Hs.eg.db", 
                                                 keyType = "ENSEMBL", ont = ont,
                                                 pvalueCutoff = 0.05, 
                                                 qvalueCutoff = 0.05,
                                                 minGSSize = 10, maxGSSize = 500)
        }
    }
}

saveRDS(timeSplitGORes, file = paste0(resultDir, "GO/fullGOSplitTimepointsWTKO.rds"))
saveRDS(otherGO, file = paste0(resultDir, "GO/fullGOMFCCBPSplitTimepointsWTKO.rds"))

# if the user has run this code before and is rerunning some analysis after this point, it is faster to load GO output from the file
if (!exists("timeSplitGORes")) {
    timeSplitGORes <- readRDS(file = paste0(resultDir, "GO/fullGOSplitTimepointsWTKO.rds"))
}

# Save GO output as CSV:
if (!dir.exists(paste0(resultDir, "GO/sepTimepointsFullCSV"))) {
  dir.create(paste0(resultDir, "GO/sepTimepointsFullCSV"))
}
for (day in days) {
    for (reg in c("up", "down")) {
        fwrite(timeSplitGORes[[day]][[reg]]@result, 
               paste0(resultDir, "GO/sepTimepointsFullCSV/", day, group, 
                      str_to_title(reg), "GOBPRes.csv"))
      
        for (ont in c("MF", "CC")) {
            fwrite(otherGO[[day]][[reg]][[ont]]@result, 
                paste0(resultDir, "GO/sepTimepointsFullCSV/", day, group, 
                    str_to_title(reg), "GO", ont, "Res.csv"))
        }
    }
}

# save plots of top 15 GO terms
pdf(file = paste0(resultDir, "GO/topBPTermsEachTimepointWTKO.pdf"))
par(mfrow = c(1,1))
for (day in days) {
    for (reg in c("up", "down")) {
      print(barplot(timeSplitGORes[[day]][[reg]], showCategory=15, 
          title = paste0("GO BP - WT vs KO ", reg, "reg ", day)))
    }
}
dev.off()

pdf(file = paste0(resultDir, "GO/topMFCCTermsEachTimepointWTKO.pdf"))
par(mfrow = c(1,1))
for (day in days) {
    for (reg in c("up", "down")) {
        print(barplot(otherGO[[day]][[reg]]$MF, showCategory=15, 
                      title = paste0("GO MF - WT vs KO ", reg, "reg ", day)))
        print(barplot(otherGO[[day]][[reg]]$CC, showCategory=15, 
                      title = paste0("GO CC - WT vs KO ", reg, "reg ", day)))
    }
}
dev.off()


# Make combined GO plots 
# Up d20   up d60   up d100
# down d20 down d60 down d100
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

p <- list()
i <- 0
for (day in days) {
    for (reg in c("up", "down")) {
        i <- i + 1
        p[[i]] <- cleanGOBarplot(timeSplitGORes[[day]][[reg]],
                       title = paste0(day, " ", reg, "regulated GO BP"),
                       textWidth = 30, showCategory = 5,
                       hjust = 1)
        
    }
}
ggarrange(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], p[[6]], 
          nrow = 3, ncol = 2, labels = "AUTO", font.label = list(size = 18))
ggsave(paste0(resultDir, "GO/topTerms.png"), device = "png",
       height = 8, width = 6, bg = "white")

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

for (day in days) {
  for (reg in c("up", "down")) {
    timeSplitEntrez <- IDToEntrez[timeSplitGeneIDs[[day]][[reg]]]
    pathwayRes$KEGG[[day]][[reg]] <- enrichKEGG(timeSplitEntrez,
                                                universe = IDToEntrez[rownames(counts)],
                                             organism = "hsa", 
                                             keyType = "ncbi-geneid",
                                             pvalueCutoff = 0.05, 
                                             qvalueCutoff = 0.05,
                                             minGSSize = 10, maxGSSize = 500)
    pathwayRes$WP[[day]][[reg]] <- enrichWP(timeSplitEntrez,
                                              organism = "Homo sapiens")
  }
}

IDToNameFun <- function(IDs) {
  unname(IDToName[IDs])
}

entrezToNameFun <- function(IDs) {
  unname(entrezToName[IDs])
}
entrezToIDFun <- function(IDs) {
  unname(entrezToID[IDs])
}


# ---- 11. save shinyapp data ----
fullPathway <- list()
pathwayMethods <- c("KEGG")
for (method in pathwayMethods) {
  pdf(file = paste0(resultDir, "pathways/top", method, "EachTimepointWTKO.pdf"))
  par(mfrow = c(1,1))
  for (day in days) {
    for (reg in c("up", "down")) {
      print(barplot(pathwayRes[[method]][[day]][[reg]], showCategory=15, 
                    title = paste0(method, " - WT vs KO ", reg, "reg ", day)))
    }
  }
  dev.off()
  
  saveRDS(pathwayRes[[method]], file = paste0(resultDir, "pathways/full",
                                              method, "SplitTimepointsWTKO.rds"))
  
  # save as csv output
  if (!dir.exists(paste0(resultDir, "pathways/", method ,"SepTimepointsCSV"))) {
    dir.create(paste0(resultDir, "pathways/", method ,"SepTimepointsCSV"))
  }
  for (day in days) {
    for (reg in c("up", "down")) {
      # add gene IDs and names based on "/"-separated ensemblIDs 
      curRes <- pathwayRes[[method]][[day]][[reg]]@result %>%
        mutate(geneNames = lapply(lapply(str_split(geneID, "/"), entrezToNameFun),toString),
               ensembIDs = lapply(lapply(str_split(geneID, "/"), entrezToIDFun), toString))
      fwrite(curRes, paste0(resultDir, "pathways/", method ,"SepTimepointsCSV/",
                             day, group, str_to_title(reg), method, "Res.csv"))
      fullPathway[[method]][[day]][[reg]] <- curRes
    }
  }
}
fullShinyGO <- list()
for (ont in c("MF", "BP", "CC")) {
  for (day in days) {
    for (reg in c("up", "down")) {
      # add gene IDs and names based on "/"-separated ensemblIDs 
      curRes <- otherGO[[day]][[reg]][[ont]]@result %>%
        mutate(geneNames = lapply(lapply(str_split(geneID, "/"), IDToNameFun),toString))
      fullShinyGO[[ont]][[day]][[reg]] <- curRes
    }
  }
}

saveRDS(fullShinyGO, file = paste0(shinyDir, "DEGOShiny.rds"))
saveRDS(fullPathway, file = paste0(shinyDir, "pathwayData2.rds"))
barplot(timeSplitWP[[day]][[reg]], showCategory=15)

# provides lookup across all timepoints and regulation levels and both KEGG
# and WP (currently paused WP compatability).
# Searches for the given term in the term descriptions
getPathwayRes <- function(term) {
  infoCols <- c("method", "day", "reg", "sig")
  inPathways <- data.table(matrix(ncol = 17, nrow = 0))
  colnames(inPathways) <- c(infoCols, "category", "subcategory", "ID",
                            "Description", "GeneRatio", "BgRatio", "pvalue",
                            "p.adjust", "qvalue", "geneID", "Count","geneNames",
                            "ensembIDs")
  for (method in pathwayMethods) {
    for (day in days) {
      for (reg in c("up", "down")) {
        curRes <- read.csv(paste0(resultDir, "pathways/", method, "SepTimepointsCSV/",
                                  day, group, str_to_title(reg), method, "Res.csv"))
        
        matches <- curRes[grepl(term, curRes$Description, ignore.case = TRUE), ] %>% 
          mutate(method = method, day = day, reg = reg, sig = p.adjust < 0.05) %>%
          dplyr::select(all_of(infoCols), everything())
        inPathways <- rbind(inPathways, matches, fill = TRUE)
      }
    }
  }
  print(inPathways[inPathways$sig, c(1:3, 8, 12)])
  print(paste(sum(inPathways$sig), "out of", nrow(inPathways)))
  return(inPathways[inPathways$sig, c(1:3, 8, 12)])
}

getGORes <- function(term) {
    infoCols <- c("ont", "day", "reg", "sig")
    inGO <- data.table(matrix(ncol = 13, nrow = 0))
    colnames(inGO) <- c(infoCols, "ID",
                              "Description", "GeneRatio", "BgRatio", "pvalue",
                              "p.adjust", "qvalue", "geneID", "Count")
    for (ont in c("BP", "MF", "CC")) {
        for (day in days) {
            for (reg in c("up", "down")) {
                curRes <- read.csv(paste0(resultDir, "GO/sepTimepointsfullCSV/", day, "WTvsKO",
                                          str_to_title(reg), "GO", ont, "Res.csv"))
                
                matches <- curRes[grepl(term, curRes$Description, ignore.case = TRUE), ] %>% 
                    mutate(ont = ont, day = day, reg = reg, sig = p.adjust < 0.05) %>%
                    dplyr::select(all_of(infoCols), everything())
                inGO <- rbind(inGO, matches, fill = TRUE)
            }
        }
    }
    print(inGO[inGO$sig, c(1:3, 6, 10)])
    print(paste(sum(inGO$sig), "out of", nrow(inGO)))
    return(inGO[inGO$sig, c(1:3, 6, 10)])
}

curRes <- pathwayRes[[method]][[day]][[reg]]@result
curRes %>% grepl("ASD", curRes$Description)
inPathways2 <- list()
test <- "WNT"
inPathways2[[test]] <- getPathwayRes(test)

inGOs2 <- list()
test <- "WNT"
test <- "BMP" # 5 in down D20 and D60 each
test <- "hedgehog" # 0
test <- "GABA" # 3, 22
test <- "dorsal" # 7/43 - up and down
test <- "SHH" # 
test <- "smoothen" # none
test <- "ventral"


inPathways2[[test]] <- getPathwayRes(test)
inGOs2[[test]] <- getGORes(test)

getGORes(test)

thesisTestingDE <- list()
thesisTestingDE$GO <- inGOs2
thesisTestingDE$Pathway <- inPathways2

saveRDS(thesisTestingDE, paste0(resultDir, "GO/DEThesisEnrichmentTest.rds"))


fwrite(inGOs[[test]], paste0(resultDir, "GO/wnt/inGO.csv"))
fwrite(inPathways[[test]], paste0(resultDir, "GO/wnt/inKEGG.csv"))

inPathways$DNA <- getPathwayRes("DNA")
inPathways$ASD <- getPathwayRes("ASD")



# ---- Answer question: is the switch to ventral progenitors toward an MGE or LGE type of progenitors ----
getGeneDEs <- function(geneList) {
  res <- data.frame(genes = geneList)
  res$IDs <- nameToID[res$genes]
  
  for (day in c("D20", "D60", "D100")) {
    sigDE <- read.csv(paste0(resultDir, "DE/DEGenes", day, "WTvsKO.csv"))
    setDT(sigDE) 
    res[[day]] <- sigDE[match(res$IDs, sigDE$geneID)]$log2FoldChange
  }
  return(res)
}

# Set up reverse look up to IDToName
nameToID <- names(IDToName)
names(nameToID) <- IDToName[nameToID]

# MGE-Specific Marker Genes:
MGEs <- getGeneDEs(c("NKX2-1", "LHX6", "GAD1", "GAD2", "DLX1", "DLX2",
                            "SST", "LHX8", "CRABP1"))

# LGE-Specific Marker Genes:
LGEs <- getGeneDEs(c("MEIS2", "EBF1", "PCP4", "ISL1", "ZFHX3", "ZFHX4"))

# hindbrain markers:
HBs <- getGeneDEs(c("EGR2", "IRX3", "PHOX2B", "LHX1", "HOXA2", "HOXB2",
                    "HOXB3", "PAX2", "FGF8", "BARHL1", "EN1", "OTX2", 
                    "EN2", "GBX2", "HOXB1", "HOXB4"))
                  
# ---- save supplementary ----

IDToNameFun <- function(IDs) {
  unname(IDToName[IDs])
}

DE <- list()
for (day in c("D20", "D60", "D100")) {
  for (reg in c("up", "down")) {
    curData <- data.frame(matrix(nrow = 0, ncol = 12))
    colnames(curData) <- c("enrichmentType", "enrichmentSuptype", "ID",	"Description",	"GeneRatio",	"BgRatio",
                           "pvalue",	"p.adjust",	"qvalue",	"geneID",	"Count", "geneNames")
    
    for (ont in c("BP", "MF", "CC")) {
      # add geneName data
      curCSV <- read.csv(paste0(resultDir, "GO/sepTimepointsFullCSV/", day, "WTvsKO", 
                                str_to_title(reg), "GO", ont, "Res.csv")) %>%
        filter(p.adjust < 0.05) %>%
        mutate(geneNames = lapply(lapply(str_split(geneID, "/"), IDToNameFun), toString))
      #fwrite(curCSV, paste0(resultDir, "thesisSup/", day,
      #                      str_to_title(reg), "GO", ont, "Res.csv"))
      curCSV$enrichmentType <- "GO"
      curCSV$enrichmentSuptype <- ont
      curCSV <- curCSV %>% dplyr::select(enrichmentType, enrichmentSuptype, everything())
      curData <- rbind(curData, curCSV)
    }
    curKEGG <- read.csv(paste0(resultDir, "pathways/KEGGSepTimepointsCSV/", day, "WTvsKO", 
                              str_to_title(reg), "KEGG", "Res.csv")) %>% filter(p.adjust < 0.05)
    curKEGG$enrichmentType <- "pathways"
    curKEGG$enrichmentSuptype <- "KEGG"
    curKEGG <- curKEGG %>%
      dplyr::select(!c(geneID, subcategory, category)) %>%
      mutate(geneID = ensembIDs) %>%
      dplyr::select(enrichmentType, enrichmentSuptype, !ensembIDs) 
    curData <- rbind(curData, curKEGG)
    curData <- curData %>% .[order(.$p.adjust),]
    fwrite(curData, paste0(resultDir, "thesisSup/", day,
                                                str_to_title(reg), "MergedRes.csv"))
    
  }
}
# Excel workbook code from:
# https://www.r-bloggers.com/2019/08/creating-excel-workbooks-with-multiple-sheets-in-r/

# Get the file name read in as a column
read_filename <- function(fname) {
  read_csv(fname, col_names = TRUE) %>%
    mutate(filename = fname)
}
tbl <-
  list.files(path = paste0(resultDir, "thesisSup/"),
             pattern ="*MergedRes.csv",
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

saveWorkbook(wb, paste0(resultDir, "thesisSup/DETimepointEnrichmentResults.xlsx"), overwrite = TRUE)
