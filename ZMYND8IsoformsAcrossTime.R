#' ZMYND8 isoform analysis
#' 1. take in ipsc and organoid isoform data
#' 2. filter to only ZMYND8 isoforms with expression > 1 in min 2 samples
#' 3. test association of expression with timepoint
#' 4. plot pheatmap summary of isoform %
#' 5. test batch correction

# ---- Step 0: Set up and package management ---
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")

# plotting
library(pheatmap)
library(ggplot2)
library(ggpubr)
# data manipulation
library(dplyr)
library(tidyr)
library(data.table)
# batch correction
library(limma)

zmynd8Dir <- "Z:/mnt/Data1/PROJECTS/ZMYND8/"
inDir <- paste0(zmynd8Dir, "ZMYND8_isoform_quant/")
ips <- read.table(paste0(inDir, "results_ipsc/star_salmon/salmon.merged.transcript_tpm.tsv"), header=TRUE)
ipsSamples <- read.csv(paste0(zmynd8Dir, "iPSC_BulkRNA-seq/SampleInfo.csv"))

for (j in c(1:nrow(ipsSamples)))
{
  k <- grep(ipsSamples$Sample.Name[j], colnames(ips))
  colnames(ips)[k] <- ipsSamples$External.ID[j]
}
org <- read.table(paste0(inDir, "results_bulkOrg/star_salmon/salmon.merged.transcript_tpm.tsv"),
               header=TRUE)
g <- "ENSG00000101040" # ZMYND8 ensembl ID
ips <- ips[grep(g, ips$gene_id),]
org <- org[grep(g, org$gene_id),]
colnames(org) <- gsub(".markdup.sorted", "", colnames(org))

merged <- merge(ips, org[,-2], by="tx")
merged$tx <- gsub("ENSG00000101040.15", "", merged$tx)
merged$tx <- gsub("|20|-", "", merged$tx)
merged$tx <- gsub("|", "", merged$tx, fixed=TRUE)
rownames(merged) <- merged$tx

# subset for WT samples
merged <- merged[,c(1,2,grep("wt", colnames(merged), ignore.case = TRUE))]
#filter for > 1TPM in at least 2 samples
merged <- merged[rowSums(merged[, -c(1:2)] > 1) >=2, ]

# update sample names
colnames(merged)[-c(1,2)] <- colnames(merged[, -c(1,2)]) %>%
  gsub("_GOK.*", "", .) %>%
  gsub("D2.", "D20", .) %>%
  gsub("wt", "iPSC_WT", .)

# ---- linear model ----
expr <- merged[, -c(1:2)]
logExpr <- log2(expr + 0.1)
timepoint <- colnames(logExpr) %>%
  gsub("_WT.*", "", .) %>%
  gsub("B._", "", .) %>%
  gsub("D[61]00?", "D60-100", .) %>%
  as.factor() 
sigRes <- data.frame(matrix(nrow = nrow(logExpr),
                            ncol = 6))
rownames(sigRes) <- rownames(logExpr)
colnames(sigRes) <- paste(rep(c("iPSCvsD20", "D20vsD60-100"), each = 3),
                          rep(c("t", "pval", "padjusted"), 2), sep = ".")
# check normal dist
hist(unlist(logExpr, use.names = F))
for (tscript in rownames(expr)) {
  print(tscript)
  tpLm <- lm(t(logExpr)[,tscript] ~ timepoint)
  lmSum <- summary(tpLm)
  D20vsiPSC <- lmSum$coefficients["timepointiPSC",]
  # when swapping order of comparison, estimate and t val sign is swapped
  iPSCvsD20 <- D20vsiPSC
  iPSCvsD20[c("Estimate", "t value")] <- D20vsiPSC[c("Estimate", "t value")] * -1
  D20vsD60100 <- lmSum$coefficients["timepointD60-100",]
  # note empty NA left for the p-adjusted val
  sigRes[tscript, ] <- c(iPSCvsD20[c("t value", "Pr(>|t|)")], NA,
                         D20vsD60100[c("t value", "Pr(>|t|)")], NA)
}

# adjust for multiple testing
sigRes$iPSCvsD20.padjusted <- p.adjust(sigRes$iPSCvsD20.pval, method = "BH")
sigRes$`D20vsD60-100.padjusted` <- p.adjust(sigRes$`D20vsD60-100.pval`, method = "BH")

fwrite(sigRes, paste0(inDir, "Kaya_thesis/lmIsoformExpVsTimepointsLog2.csv"), row.names = TRUE)

sigiPSCD20 <- rownames(sigRes[sigRes$iPSCvsD20.padjusted < 0.05,])
sigD20D60 <- rownames(sigRes[sigRes$`D20vsD60-100.padjusted` < 0.05,])

# ---- plotting ----
# transform to % of expression accounted for by each isoform
mergedTr <- merged
for (j in c(3:ncol(merged)))
  mergedTr[,j] <- (merged[,j]/sum(merged[,j]))*100

exprTr <- mergedTr[, -c(1:2)]
fwrite(exprTr, paste0(inDir, "Kaya_thesis/isoExprAsPercentage.csv"), row.names = TRUE)

pheatmap(exprTr)

# Modify ordering of the clusters using clustering callback option from pheatmap manual
callback <- function(hc, mat){
  sv <- svd(t(mat))$u[,3] # this number can be modified to change ordering
  dend <- reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

rowAnn <- data.frame(row.names = rownames(exprTr),
                     iPSCvsD20 = factor(rownames(exprTr) %in% sigiPSCD20, levels = c("TRUE", "FALSE")),
                     D20vsD60 = factor(rownames(exprTr) %in% sigD20D60, levels = c("TRUE", "FALSE")))
allTimepoints <- colnames(exprTr) %>%
  gsub("_WT.*", "", .) %>%
  gsub("B._", "", .) %>%
  factor(levels = c("iPSC", "D20", "D60", "D100")) 

colAnn <- data.frame(row.names = colnames(exprTr),
                     Timepoint = allTimepoints)

p <- pheatmap(exprTr, annotation_row = rowAnn, annotation_col = colAnn, 
         annotation_colors = list(iPSCvsD20 = c("TRUE" = "darkgreen",
                                            "FALSE" = "white"),
                                  D20vsD60 = c("TRUE" = "darkolivegreen",
                                               "FALSE" = "lightgrey")),
         annotation_names_row = FALSE,
         annotation_names_col = FALSE,
         cutree_cols = 4,
         silent = TRUE,
         clustering_callback = callback)
# replace the legend 2 and 3 names with descriptive labels
p$gtable$grobs[[8]]$children[[4]]$label <- "Signficant\nD20-D60\nswitch"
p$gtable$grobs[[8]]$children[[7]]$label <- "Signficant\niPSC-D20\nswitch"

# move legend 2 down
p$gtable$grobs[[8]]$children[[5]]$vjust <- p$gtable$grobs[[8]]$children[[5]]$vjust + 2.1
p$gtable$grobs[[8]]$children[[6]]$vjust <- p$gtable$grobs[[8]]$children[[6]]$vjust + 3.9

# move legend 3 further down
p$gtable$grobs[[8]]$children[[7]]$vjust <- p$gtable$grobs[[8]]$children[[7]]$vjust + 1
p$gtable$grobs[[8]]$children[[8]]$vjust <- p$gtable$grobs[[8]]$children[[8]]$vjust + 4.5
p$gtable$grobs[[8]]$children[[9]]$vjust <- p$gtable$grobs[[8]]$children[[9]]$vjust + 8.75
ggarrange(p$gtable)
ggsave(paste0(inDir, "Kaya_thesis/isoformPheatmap.png"), device = "png",
      height = 6, width = 7, bg = "white")


# using batch corrected data instead
merged <- merge(ips, org[,-2], by="tx")
merged$tx <- gsub("ENSG00000101040.15", "", merged$tx)
merged$tx <- gsub("|20|-", "", merged$tx)
merged$tx <- gsub("|", "", merged$tx, fixed=TRUE)
rownames(merged) <- merged$tx
colnames(merged)[-c(1,2)] <- colnames(merged[, -c(1,2)]) %>%
  gsub("_GOK.*", "", .) %>%
  gsub("D2.", "D20", .) %>%
  gsub("wt", "iPSC_WT", .)

merged[, -c(1,2)] <- log2(merged[, -c(1,2)] + 0.1)
merged <- merged[, grep("WT", colnames(merged), ignore.case = TRUE)]
batch <- c(rep(1, 3), rep(2, 9))
adjusted_merged.wt <- removeBatchEffect(merged, batch = batch)
adjusted_merged.wt <- adjusted_merged.wt[rownames(expr),]

#Supp Fig
pdf(paste0(inDir, "Kaya_thesis/supplementaryBatchEffects.pdf"), height = 5, width = 7)
pheatmap(cor(merged[rownames(expr),]), main = "cor(logExpr)\n")
pheatmap(cor(adjusted_merged.wt), main = "cor(batch adjusted logExpr)\n")
pheatmap(adjusted_merged.wt, main = "batch adjusted logExpr\n")
dev.off()

timepoint <- colnames(logExpr) %>%
  gsub("_WT.*", "", .) %>%
  gsub("B._", "", .) %>%
  gsub("D[61]00?", "D60-100", .) %>%
  as.factor() 
sigRes2 <- data.frame(matrix(nrow = nrow(adjusted_merged.wt),
                            ncol = 6))
rownames(sigRes2) <- rownames(adjusted_merged.wt)
colnames(sigRes2) <- paste(rep(c("iPSCvsD20", "D20vsD60-100"), each = 3),
                          rep(c("t", "pval", "padjusted"), 2), sep = ".")
# check normal dist
hist(unlist(adjusted_merged.wt, use.names = F))
for (tscript in rownames(adjusted_merged.wt)) {
  print(tscript)
  tpLm <- lm(t(adjusted_merged.wt)[,tscript] ~ timepoint)
  lmSum <- summary(tpLm)
  D20vsiPSC <- lmSum$coefficients["timepointiPSC",]
  # when swapping order of comparison, estimate and t val sign is swapped
  iPSCvsD20 <- D20vsiPSC
  iPSCvsD20[c("Estimate", "t value")] <- D20vsiPSC[c("Estimate", "t value")] * -1
  D20vsD60100 <- lmSum$coefficients["timepointD60-100",]
  # note empty NA left for the p-adjusted val
  sigRes2[tscript, ] <- c(iPSCvsD20[c("t value", "Pr(>|t|)")], NA,
                         D20vsD60100[c("t value", "Pr(>|t|)")], NA)
}

# adjust for multiple testing
sigRes2$iPSCvsD20.padjusted <- p.adjust(sigRes2$iPSCvsD20.pval, method = "BH")
sigRes2$`D20vsD60-100.padjusted` <- p.adjust(sigRes2$`D20vsD60-100.pval`, method = "BH")
fwrite(sigRes2, paste0(inDir, "Kaya_thesis/lmIsoformBatchCorrected.csv"), row.names = TRUE)
sigiPSCD20 <- rownames(sigRes2[sigRes2$iPSCvsD20.padjusted < 0.05,])
sigD20D60 <- rownames(sigRes2[sigRes2$`D20vsD60-100.padjusted` < 0.05,])

rowAnn <- data.frame(row.names = rownames(exprTr),
                     iPSCvsD20 = factor(rownames(exprTr) %in% sigiPSCD20, levels = c("TRUE", "FALSE")),
                     D20vsD60 = factor(rownames(exprTr) %in% sigD20D60, levels = c("TRUE", "FALSE")))


p <- pheatmap(adjusted_merged.wt, annotation_row = rowAnn, annotation_col = colAnn, 
              annotation_colors = list(iPSCvsD20 = c("TRUE" = "darkgreen",
                                                     "FALSE" = "white"),
                                       D20vsD60 = c("TRUE" = "darkolivegreen",
                                                    "FALSE" = "lightgrey")),
              annotation_names_row = FALSE,
              annotation_names_col = FALSE,
              silent = TRUE,
              cutree_cols = 4,
              clustering_callback = callback)
ggarrange(p$gtable)
ggsave(paste0(inDir, "Kaya_thesis/isoformPheatmapSupplementary.png"), device = "png",
       height = 6, width = 7, bg = "white")
