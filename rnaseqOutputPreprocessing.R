# Overview from https://docs.google.com/document/d/1AJ1IW-3AZ9_yB2gNR8Ri7irl-GU-JHTVbN1nZKAPoV8/edit
# In (RNA-seq data): 
# - Gene Counts (for DEseq)
# - Gene TPM (for WGCNA)
# - Transcript Counts/TPM (for ZMYND8 isoforms)
# Perform:
# - Filtering
# - Log2 transformation
# - Normalize to RPKM
# - PCA
# - Format metadata and save 

# Proposed metadata structure:
# R data structure: {
#     Unfiltered gene counts,
#     Filtered gene counts,
#     filtered RPKM,
#     log2 transformed RPKM filtered,
#     Gene TPM of filtered genes,
#     GeneID to gene name,
#     Transcript TPM,
#     Transcript counts,
#     GeneID to gene name,
#     Metadata
#}

# Proposed program structure:
# 1. Take in gene counts, gene lengths, gene TPM, and transcript counts
# 2. Create pre-filtered DESeq2 data, counts with min 10 in min 3 samples
# 3. Calculate RPKM based on gene counts and gene lengths
# 4. Threshold RPKM at 0.5 in min 2 samples
# 5. log2 transform RPKM 
# 6. Save RPKM and gene counts to metadata structure
# 7. Search genbank and ensembl for ensembl IDs and save
# 8. Output PCA plots
# 9. Process transcript data [TODO]


# ---- Step 0: Set up and package management ---
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")


BiocManager::install("genefilter")
BiocManager::install("biomaRt")
BiocManager::install("furrr")
BiocManager::install("clipr")
BiocManager::install("tidyr")
BiocManager::install("ggpubr")


install.packages("dplyr")
install.packages("tidyverse") # data manipulation
install.packages("ggplot2")


library(ggplot2)
library(data.table)
library(dplyr)
library(tidyverse)
library(tidyr)
library(biomaRt)


library(varhandle)
library(magrittr) 

inputDir <- "Z:/mnt/Data1/PROJECTS/ZMYND8/Organoid_BulkRNA-seq/NFcore_results/star_salmon/"
outputDir <- "Z:/mnt/Data1/PROJECTS/ZMYND8/Organoid_BulkRNA-seq/Analysis_Results_Kaya/preprocessing/"

# ---- 1. Take in gene counts, gene lengths, gene TPM, and transcript counts ----
fullGeneCountsData <- readRDS(paste0(inputDir, "salmon.merged.gene_counts.rds"))
counts <- fullGeneCountsData@assays@data@listData[["counts"]]

# if this script has been run before, there is no need to recreate everything
if(file.exists(paste0(outputDir, "preprocessedData.rds"))){
  preprocessedData <- readRDS(paste0(outputDir, "preprocessedData.rds"))
} else {
  preprocessedData <- list()
}
preprocessedData$unfilteredCounts <- counts
geneLengths <- read.delim(paste0(inputDir, "salmon.merged.gene_lengths.tsv"))
geneTPM <- read.delim(paste0(inputDir, "salmon.merged.gene_tpm.tsv"))

fullTranscriptCountsData <- readRDS(paste0(inputDir, "salmon.merged.transcript_counts.rds"))
transcriptTPM <- fullTranscriptCountsData@assays@data@listData[["abundance"]]
transcriptCounts <- fullTranscriptCountsData@assays@data@listData[["counts"]]

# ---- 2. Create pre-filtered DESeq2 data, counts over 10 in min 3 samples ----
filteredCounts <- counts[rowSums(counts > 10) >= 3,]
preprocessedData$filteredCounts <- filteredCounts

# ---- 3. Calculate RPKM based on gene counts and gene lengths ----
# from Javier's RPKM script
sequencingDepth <- colSums(counts)/10^6
###Pre-normalization // standardization to counts per million
cpm <- counts
for (j in 1:27)
  cpm[,j] <- cpm[,j]/sequencingDepth[j]
##remove genes without length information
sum(rowSums(is.na(geneLengths[, 2:29]))) # none without length
RPKM <- cpm %>% filter(rowSums(is.na(geneLengths[, 2:29])) == 0)
###check % of genes with length assignation
nrow(RPKM)/nrow(cpm)*100
###Compute RPKM
RPKM <- (RPKM/geneLengths[, 3:29])*1e3

# ---- 4. Threshold RPKM at 0.5 in min 2 samples ----
overThreshold <- rowSums(RPKM > 0.5) >= 2
filteredRPKM <- RPKM[overThreshold,]
preprocessedData$filteredRPKM <- filteredRPKM

# ---- 5. log2 transform RPKM (with addition of pseudocount = 1) ----
log2RPKM <- log2(filteredRPKM + 1)
preprocessedData$log2RPKM <- log2RPKM

# ---- 6. Save metadata ----
metadata <- data.frame(Sample = names(counts)) %>%
  mutate(Type = gsub("_GOK.*","", Sample) %>% gsub("T","T_522", .) %>% 
           gsub("D2._","D20_", .)) %>%
  separate_wider_delim(Type, delim = "_", 
                       names = c("Batch", "Day", "Genotype", "Clone")) %>% 
  mutate(Genotype = factor(.$Genotype), Batch = factor(.$Batch),
         Day = factor(.$Day))
preprocessedData$metadata <- metadata

# ---- 7. Output PCA plot ---- 
data <- log2RPKM
day <- metadata$Day[metadata$Sample == colnames(data)]
group <- metadata$Group[metadata$Sample == colnames(data)]
pca <- princomp(data, cor = TRUE)$loadings

pdf(file = paste0(outputDir, "PCA.pdf"))
par(mfrow = c(1,1))
ggplot(pca[, 1-2], aes(x = pca[,1],y = pca[,2], color = day, shape = group)) +
  geom_point() +
  scale_color_discrete(limits = c('D20', 'D60', 'D100')) +
  scale_shape_manual(values = c(4, 16, 17)) +
  labs(x = "PC1", y = "PC2", title = "PCA on log2 transformed filtered RPKM", 
       color = "Day", shape = "Group")
dev.off()

# ---- 8. Threshold TPM at 0.5 in min 2 samples----
WTKO <-  metadata$Sample[metadata$Genotype != "Het"]
longer <- log10(geneTPM[-infoCols, WTKO]) %>% 
  as_tibble(rownames = "GeneID") %>%
  pivot_longer(-1, names_to = "Sample", values_to = "log10Expression")

ggplot(longer, aes(x = Sample, y = log10Expression)) +
        geom_violin() +
        geom_point() +
        theme(axis.text.x = element_text(angle = 90)) +
        labs(title = paste("Normalised expression"))

ggplot(log10(geneTPM[-infoCols, WTKO]), aes(x = Sample)) +
  geom_histogram()
ggplot(longer) +
  geom_histogram(aes(x = log10Expression)) +
  facet_wrap(~ Sample)

hist(longer$log10Expression, )
preprocessedData$geneTPM <- geneTPM
infoCols <- c(1,2)
overThresholdTPM <- rowSums(geneTPM[-infoCols] > 0.5) >= 2
sum(rowSums(geneTPM[-infoCols] > 0.5) >= 2)

sum(rowSums(geneTPM[-infoCols] > 1) >= 2)
sum(rowSums(geneTPM[-infoCols] > 5) >= 2)
sum(rowSums(geneTPM[-infoCols] > 10) >= 2)

over5TPM <- geneTPM[rowSums(geneTPM[-infoCols] > 5) >= 2, ]
over1TPM <- geneTPM[rowSums(geneTPM[-infoCols] > 1) >= 2, ]

filteredGeneTPM <- geneTPM[overThresholdTPM, ]
preprocessedData$filteredGeneTPM <- filteredGeneTPM

# ---- 9. log2 transform TPM (with addition of pseudocount = 1) ----
log2GeneTPM <- filteredGeneTPM
log2GeneTPM[-infoCols] <- log2(filteredGeneTPM[-infoCols] + 1)
preprocessedData$log2GeneTPM <- log2GeneTPM


# ---- 10. Save ensemblID to gene name conversion ----
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- row.names(preprocessedData$unfilteredCounts)
hgncSym <- getBM(filters = "ensembl_gene_id", 
                 attributes = c("ensembl_gene_id","hgnc_symbol"),
                 values = genes, mart = mart)

geneIDToName <- hgncSym$hgnc_symbol
# when the gene name isn't known, keep the ensemblID
emptyNames <- is.na(hgncSym$hgnc_symbol) | hgncSym$hgnc_symbol == ""
geneIDToName[emptyNames] <- hgncSym$ensembl_gene_id[emptyNames]
names(geneIDToName) <- hgncSym$ensembl_gene_id

# some (127) ensemblIDs are not recognised by getBM and are lost
# include these with their ID as the gene name
unknownIDs <- genes[!genes %in% hgncSym$ensembl_gene_id]
geneIDToName[unknownIDs] <- unknownIDs

# some gene names are duplicates - should these be saved as ensemblIDs too?

preprocessedData$IDToName <- geneIDToName

# ---- 11. Save transcript Counts/TPM ---- 
preprocessedData$transcriptTPM <- transcriptTPM
preprocessedData$transcriptCounts <- transcriptCounts

# ---- 12. Export data as rds for processing and csv for use in supplementary tables
saveRDS(preprocessedData, file = paste0(outputDir, "preprocessedData.rds"))



