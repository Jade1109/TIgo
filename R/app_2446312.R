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
source("seurat.R")
source("TI.R")
source("go.R")
options(shiny.maxRequestSize = 1000*1024^2)
# Combine UIs using navbarPage
ui <- navbarPage(
  "TIgo",
  tabPanel("Seurat object creation", umapui("umap_module")),
  tabPanel("Trajectory Analysis", TIUI("trajectory_module")),
  tabPanel("GO Enrichment Analysis", goUI("go_module"))
)


# Combined server function
server <- function(input, output, session) {
  # Call module server functions for each tab
  callModule(umapserver, "umap_module")
  callModule(TIserver, "trajectory_module")
  callModule(goserver,"go_module")
}

# Run the combined Shiny app
shinyApp(ui = ui, server = server)

