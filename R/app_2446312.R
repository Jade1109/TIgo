#' User Interface for the Trajectory Analysis Module
#'
#' @param id The Shiny module id.
#' @return A Shiny UI for the Trajectory Analysis Module.
#' @export
#' @importFrom shiny NS fluidPage sidebarLayout sidebarPanel mainPanel h4 fileInput selectInput actionButton textAreaInput numericInput
#' @importFrom shiny uiOutput verbatimTextOutput downloadButton plotlyOutput tableOutput fluidRow div
#' @importFrom plotly ggplotly
#' @importFrom dplyr select
#' @importFrom ggplot2 ggsave
#' @importFrom tools file_ext
#' @importFrom zip zip
#' @importFrom SeuratObject readRDS
#' @importFrom slingshot slingshot
#' @importFrom harmony harmony
#' @importFrom tibble tibble
#' @importFrom purrr map
#' @importFrom dynwrap dynwrap
#' @importFrom dyro dyro
#' @importFrom ggplot2 ggplot
#' @importFrom plotly plot_ly
#' @importFrom clusterProfiler enrichGO
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom pheatmap pheatmap
#' @importFrom shinyWidgets shinyWidgets
#' @importFrom monocle3 monocle3
#' @importFrom waiter waiter
#' @importFrom mgcv mgcv
#' @importFrom broom broom
#' @importFrom dplyr dplyr
#' @importFrom tidyr tidyr
#' @importFrom RColorBrewer RColorBrewer
#' @importFrom DT DT
#'
library(slingshot)
library(Seurat)
library(magrittr)
library(harmony)
library(tibble)
library(purrr)
library(dynwrap)
library(dyno)
library(ggplot2)
library(plotly)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)
library(shinyWidgets)
library(monocle3)
library(waiter)
library(mgcv)
library(broom)
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(DT)
# Source the separate app scripts
source("/Users/apple/Desktop/TIgo/R/seurat.R")
source("/Users/apple/Desktop/TIgo/R/TI.R")
source("/Users/apple/Desktop/TIgo/R/go.R")
options(shiny.maxRequestSize = 1000*1024^2)
#' TIgo UI
#'
#' Creates the user interface for the TIgo Shiny app using a `navbarPage` layout. The UI consists of three tabs:
#' - Seurat object creation
#' - Trajectory Analysis
#' - GO Enrichment Analysis
#'
#' @export
ui <- navbarPage(
  "TIgo",
  tabPanel("Seurat object creation", umapui("umap_module")),
  tabPanel("Trajectory Analysis", TIUI("trajectory_module")),
  tabPanel("GO Enrichment Analysis", goUI("go_module"))
)


#' TIgo Server
#'
#' Defines the server logic for the TIgo Shiny app. It includes server functions for each of the three modules corresponding to the tabs in the UI.
#'
#' @param input, output, session Standard Shiny server arguments.
#' @export
server <- function(input, output, session) {
  # Call module server functions for each tab
  callModule(umapserver, "umap_module")
  callModule(TIserver, "trajectory_module")
  callModule(goserver,"go_module")
}
<<<<<<< HEAD
=======


>>>>>>> origin/main
#' Run TIgo Shiny App
#'
#' This function runs the combined Shiny app defined by the UI and server functions.
#'
#' @export
runTIgoApp <- function() {
  shinyApp(ui = ui, server = server)
}

