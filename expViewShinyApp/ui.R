#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(shinyFiles)
library(magrittr)
library(shinycssloaders)
library(shinydashboard)

#shinyData <- readRDS("shinyData.rds")

fluidPage(
  
  # Application title
  titlePanel("Organoid expression - bulk RNA-seq"),
  
  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(width=3, style = paste0("height: 80vh; overflow: hidden;"),
                 shinyFilesButton("TPMIn", "Upload RNA-seq TPM", title = "Select TPM .csv or .tsv file:", multiple = FALSE, buttonType = "default", class = NULL, style="margin-bottom: 1em;"),
                 shinyFilesButton("metaIn", "Upload RNA-seq meta" , title = "Select TPM .csv or .tsv file:", multiple = FALSE, buttonType = "default", class = NULL, style="margin-bottom: 1em;"),
                 selectizeInput("gene", "Select gene", choices = NULL),
                 selectInput("xAxis", "Select x-axis var", choices = NULL),
                 selectInput("colorBy", "Select colour by var", choices = NULL),
                 uiOutput("xAxisGroups"),
                 uiOutput("colorGroups"),
    ),
    
    mainPanel(width=9, style = paste0("overflow: hidden;"),
              fluidRow(style = "padding: 10px;",
                       tabBox(id = "tabset1", height = "100%", width = "100%",
                              tabPanel("Main Plots", 
                                       fluidRow(
                                         plotOutput("exprPlot", height="80vh") %>% withSpinner()
                                       )
                              )
                              # for the thesis only leave main functionality 
                              #,
                              #tabPanel("DERes Table", 
                              #         fluidRow(splitLayout(cellWidths = c("50%", "50%"), 
                              #                              plotOutput("DERes"),
                              #                              plotOutput("markerInfo"))
                              #         )
                              #)#,
                              #tabPanel("Other Plots", splitLayout(cellWidths = c("30%", "30%", "30%"),
                              #                                    plotOutput("IHC20"),
                              #                                    plotOutput("IHC60"),
                              #                                    plotOutput("IHC100"))
                              #)
                       )
              )
    )
  )
)
