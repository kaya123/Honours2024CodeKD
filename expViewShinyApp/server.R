#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)
library(rsconnect)
library(vroom)
library(tibble)
library(shinyFiles)


values <- reactiveValues()

# ---- 1. Take in data from organoid paper and from WT DE -----
shinyData <- readRDS("shinyData.rds")
imagingData <- readRDS("IHCPaths.rds")

#tpm <- shinyData$TPM[c("gene_id", "gene_name", meta$Sample)]
markerTables <- shinyData$markerTables
sigDEs <- shinyData$sigDEs
markers <- names(markerTables)
genes <- shinyData$TPM$gene_name
shinyData$metaData$Clone <- gsub("522", "0", shinyData$metaData$Clone)
volumes = c(wd='.') # data0='/mnt/Data0', data1='/mnt/Data1')

function(input, output, session) {
  values$TPMIn <- read.delim("salmon.merged.gene_tpm.tsv")
  values$metaIn <- shinyData$metaData
  values$colorBy <- "Genotypes"
  
  shinyFileChoose(input, 'TPMIn', root=volumes, filetypes=c('tsv', 'csv'))
  shinyFileChoose(input, 'metaIn', root=volumes, filetypes=c('tsv', 'csv'))
  
  
  observe({
    if(!is.null(input$TPMIn) && !is.integer(input$TPMIn)) {
      file_selected<-parseFilePaths(volumes, input$TPMIn)
      path<-as.character(file_selected$datapath)
      ext <- tools::file_ext(path)
      values$TPMIn <- switch(ext,
                             csv = vroom::vroom(path, delim = ","),
                             tsv = vroom::vroom(path, delim = "\t"),
                             validate("Invalid file; Please upload a .csv or .tsv file")
      )
    }
  })
  
  observe({
    if(!is.null(input$metaIn) && !is.integer(input$metaIn)) {
      file_selected<-parseFilePaths(volumes, input$metaIn)
      path<-as.character(file_selected$datapath)
      ext <- tools::file_ext(path)
      values$metaIn <- switch(ext,
                              csv = vroom::vroom(path, delim = ","),
                              tsv = vroom::vroom(path, delim = "\t"),
                              validate("Invalid file; Please upload a .csv or .tsv file")
      )
    }
  })
  
  observe({
    # assume the gene name column is the last non numeric column
    # this is not always the case
    dataCol <- names(dplyr::select_if(values$TPMIn, is.numeric))
    infoCol <- names(values$TPMIn)[!names(values$TPMIn) %in% dataCol]
    geneNameCol <- infoCol[-1]
    #print(geneNameCol)
    updateSelectizeInput(session, 'gene', choices = values$TPMIn[[geneNameCol]], server = TRUE, selected = values$TPMIn[[geneNameCol]][1])
    updateSelectInput(session, 'colorBy', choices = colnames(values$metaIn), , selected = 'Genotype')
    updateSelectInput(session, 'xAxis', choices = colnames(values$metaIn), selected = 'Day')
  })
  
  observe({
    # update the "select groups to colour/group x-axis" box
    values$colorBy <- input$colorBy
    values$xAxisBy <- input$xAxis
  })
  
  
  output$xAxis <- renderUI({
    selectInput("xAxis", "Group x-axis by", choices = colnames(values$metaIn),
                selected = "Day")
  })
  
  
  output$colorGroups <- renderUI({
    selectInput("colorGroups", paste("Select", values$colorBy, "(colour)") , choices = unique(values$metaIn[[values$colorBy]]),
                selected = unique(values$metaIn[[values$colorBy]]), multiple = TRUE)
  })
  
  output$xAxisGroups <- renderUI({
    selectInput("xAxisGroups", paste("Select", values$xAxisBy, "(x-axis)") , choices = unique(values$metaIn[[values$xAxisBy]]),
                selected = unique(values$metaIn[[values$xAxisBy]]), multiple = TRUE)
  })
  
  
  # commented out updateSelectizeInput runs much faster but doesn't change gene
  # names based on the input file
  #updateSelectizeInput(session, 'gene', choices = genelist, server = TRUE)
  output$exprPlot <- renderPlot({
    req(input$gene)
    req(input$colorGroups)
    req(input$xAxisGroups)
    
    gene <- input$gene
    
    dataCol <- names(dplyr::select_if(values$TPMIn, is.numeric))
    infoCol <- names(values$TPMIn)[!names(values$TPMIn) %in% dataCol]
    geneNameCol <- infoCol[-1]
    genes <- values$TPMIn[[geneNameCol]]
    
    curMeta <- values$metaIn %>%
      filter(.[[values$xAxisBy]] %in% input$xAxisGroups) %>%
      filter(.[[values$colorBy]] %in% input$colorGroups)
    
    if (gene %in% genes) { # TODO: case where this is false not handled
      plotdata <- data.frame(
        colnames(values$TPMIn[, dataCol]),
        t(values$TPMIn[which(genes %in% gene), dataCol]))
      colnames(plotdata) = c("Sample", "TPM")
      plotdata <- merge(plotdata, curMeta, by = "Sample")
      
      # the uploaded metadata may include days in the D20 format
      # or just numerically as 20.
      # when in the D20 format, normal sorting will place D100 first
      # so we give it an order
      if (sum(!unique(plotdata$Day) %in% c("D20", "D100", "D60")) == 0) {
        plotdata$Day <- factor(plotdata$Day, levels = c("D20", "D60", "D100"))
      }
      
      # special case if only one selected variable
      if (length(unique(curMeta$Genotype)) == 1) {
        ggplot(plotdata, aes(y=TPM, x=factor(Day))) +
          geom_boxplot(width = 0.4, size = 0.4, staplewidth = 0.3,
                       position = position_dodge(0.8)
          ) + 
          geom_dotplot(
            aes(fill = Batch, color = Batch), 
            binaxis='y', stackdir='center', dotsize = 0.6,
            position = position_dodge(0.2)
          ) +
          labs(title = paste(gene, unique(curMeta$Genotype)), x="Day")
      } else {
        ggplot(plotdata, aes(y=TPM, x=factor(.data[[values$xAxisBy]]),
                             fill=factor(.data[[values$colorBy]]))) +
          geom_boxplot(width = 0.4, size = 0.4, staplewidth = 0.3,
                       position = position_dodge(0.6)) +
          geom_dotplot(
            aes(fill = factor(.data[[values$colorBy]])), 
            binaxis='y', stackdir='center', dotsize = 0.6,
            position = position_dodge(0.6)
          ) +
          labs(title = gene, x=values$xAxisBy, fill = values$colorBy)
      }
    } else {
      return(NULL)
    }
  })
  
  # a few genes have some saved information about them to display
  output$markerInfo <- renderPlot({
    req(input$gene)
    gene <- input$gene
    if (gene %in% markers) {
      t1 <- tableGrob(t(markerTables[[gene]]), theme = ttheme_default(base_size=10)) 
      #t2 <- markerTables[[gene]]
      grid.arrange(t1)
    }
  })
  
  # WT DE results table for the current gene
  output$DERes <- renderPlot({
    req(input$gene)
    gene <- input$gene
    if (gene %in% rownames(sigDEs)) {
      t2 <- tableGrob(t(sigDEs[gene,]), theme = ttheme_default(base_size=10))
      #t2 <- sigDEs[gene,]
      grid.arrange(t2)
    }
  })
  
  # imaging stuff - not the focus right now
  # pax6, tspan6
  output$IHC20 <- renderImage({
    req(input$gene)
    curGene <- input$gene
    curImg <- imagingData %>% filter(geneName == curGene, day == "D20", genotype == "WT")
    #print(nrow(curImg))
    req(nrow(curImg) == 1)
    filename <- normalizePath(curImg$shinyPath[1])
    return(list(src = filename, alt = "D20 image WT", width=150))
  }, deleteFile = FALSE)
  
  output$IHC60 <- renderImage({
    curGene <- input$gene
    curImg <- imagingData %>% filter(geneName == curGene, day == "D60", genotype == "WT")
    req(nrow(curImg) == 1)
  
    filename <- normalizePath(curImg$shinyPath[1])
    return(list(src = filename, alt = "D60 image WT", width=150))
  }, deleteFile = FALSE)
  
  output$IHC100 <- renderImage({
    curGene <- input$gene
    curImg <- imagingData %>% filter(geneName == curGene, day == "D100", genotype == "WT")
    req(nrow(curImg) == 1)
    filename <- normalizePath(curImg$shinyPath[1])
    return(list(src = filename, alt = "D100 image WT", width=150))
  }, deleteFile = FALSE)
}