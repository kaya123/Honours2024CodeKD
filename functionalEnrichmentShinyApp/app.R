#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(dplyr)
library(tidyverse)
library(data.table)
library(stringr)
library(DOSE)
#rm(list=ls())
#setwd("//rna2/share/mnt/Data1/PROJECTS/ZMYND8/Organoid_BulkRNA-seq/Scripts_Kaya/shiny/pathwayEnrichment")
allData <- list()
allData$DE$GO <- readRDS("DEGOShiny.rds")
allData$DE$path <- readRDS("pathwayData2.rds")
allData$WGCNA$GO <- readRDS("GOWGCNA.rds")
allData$WGCNA$path <- readRDS("pathwayWGCNA.rds")

# getting the direction of regulation from WT to KO in modules
MEs <- read.csv("ME.csv")
modDir <- rep("none", ncol(MEs) - 1)
names(modDir) <- colnames(MEs)[-1] 
padjMEs <- read.csv("padjustedMEs.csv", row.names = 1)
colnames(MEs)[-1] 
sigMEs <- colnames(MEs)[-1][padjMEs[colnames(MEs)[-1],]$Genotype < 0.05]
up <- sigMEs[colSums(MEs[grep("KO", MEs$Sample), sigMEs]) > colSums(MEs[grep("WT", MEs$Sample), sigMEs])]
down <- sigMEs[!sigMEs %in% up]
modDir[up] <- "up"
modDir[down] <- "down" 
names(modDir) <- gsub("E", "", names(modDir))

# note - modules was manually copied across, check for compatability when rerunning WGCNA
modules <- read.csv("modules.csv")
colorToM <- modules$Label %>% gsub("E", "", .)
names(colorToM) <- modules$Color

infoCols <- c("method", "dayMod", "reg", "sig")
pathwayMethods <- c("KEGG")
goOnts <- c("BP", "CC", "MF")
days <- c("D20", "D60", "D100")
mods <- names(allData$WGCNA$GO$BP)
options(digits = 3, scipen = -2)

getEnrichmentRes <- function(term, pathwayOrGO = c("path", "GO"), type = c("WGCNA", "DE")) {
  fullIn <- allData[[type]][[pathwayOrGO]]
  if (pathwayOrGO == "path") {
    inEnrichment <- data.table(matrix(ncol = 17, nrow = 0))
    colnames(inEnrichment) <- c(infoCols, "category", "subcategory", "ID",
                              "Description", "GeneRatio", "BgRatio", "pvalue",
                              "p.adjust", "qvalue", "geneID", "Count","geneNames",
                              "ensembIDs")
    methods <- pathwayMethods
  } else {
    inEnrichment <- data.table(matrix(ncol = 14, nrow = 0))
    colnames(inEnrichment) <- c(infoCols, "ID",	"Description",	"GeneRatio",	"BgRatio",	"pvalue",	
                        "p.adjust",	"qvalue",	"geneID",	"Count",	"geneNames")
    methods <- goOnts
  }
  
  stopifnot(type %in% c("WGCNA", "DE"))
  
  for (method in methods) {
    if (type == "DE") {
      for (day in days) {
        for (reg in c("up", "down")) {
          curRes <- fullIn[[method]][[day]][[reg]]
          matches <- curRes[grepl(term, curRes$Description, ignore.case = TRUE), ] %>% 
            mutate(method = method, dayMod = day, reg = reg, sig = p.adjust < 0.05) %>%
            dplyr::select(all_of(infoCols), everything())
          inEnrichment <- rbind(inEnrichment, matches, fill = TRUE)
        }
      }
    } else { # WGCNA
      for (mod in mods) {
        curRes <- fullIn[[method]][[mod]]
        matches <- curRes[grepl(term, curRes$Description, ignore.case = TRUE), ] %>% 
          mutate(method = method, dayMod = colorToM[mod], reg = modDir[colorToM[mod]], sig = p.adjust < 0.05) %>%
          dplyr::select(all_of(infoCols), everything())
        inEnrichment <- rbind(inEnrichment, matches, fill = TRUE)
      }
    }
  }
  inEnrichment <- data.frame(inEnrichment)
  inEnrichment[, c(infoCols, "Description", "p.adjust", "geneNames", "GeneRatio")]
}

getAllRes <- function(term) {
  allEnrichment <- rbind(getEnrichmentRes(term, "path", "DE"),
                         getEnrichmentRes(term, "path", "WGCNA"),
                         getEnrichmentRes(term, "GO", "DE"),
                         getEnrichmentRes(term, "GO", "WGCNA"))
  
  sigRatio <- paste(sum(allEnrichment$sig), "significant results out of", nrow(allEnrichment))
  displayData <- allEnrichment[allEnrichment$sig,] %>%
    mutate(p.adjust = signif(p.adjust, digits = 3)) %>%
    select(-c(sig, GeneRatio))
  colnames(displayData)[2] <- "day/mod"
  return(list(table=displayData, text=sigRatio))
}

ui <- fluidPage(

    titlePanel("WT vs KO DESeq2 and WGCNA enrichment results: pathway (KEGG) and GO (BP, CC, MF)"),

    sidebarLayout(
        sidebarPanel(
            textInput("term",
                        "Enter search term: ",
                        value = "Wnt"),
            textOutput("sigRatio"),
            width = 2
        ),

        mainPanel(
          DT::DTOutput("sigTable", ),
          width = 10
        )
    )
)

server <- function(input, output) {
    output$sigTable <- DT::renderDT({getAllRes(input$term)$table},
                                    options = list(
                                      autoWidth = TRUE,
                                      scrollX = TRUE,
                                      columnDefs = list(list(width = "470px", targets = c(5)))
                                    ), rownames = FALSE)
    output$sigRatio <- renderText({getAllRes(input$term)$text})
}

# Run the application 
shinyApp(ui = ui, server = server)
