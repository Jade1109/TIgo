#' Define the UMAP Module UI
#'
#' This function creates the user interface for the UMAP module of the Shiny app,
#' including file inputs, text inputs, sliders, buttons, and output elements.
#'
#' @param id A unique identifier for the module. It is used to namespace the input and output IDs.
#' @return A `shiny.tag` object representing the UI for the UMAP module.
#' @importFrom shiny fluidPage sidebarLayout sidebarPanel mainPanel uiOutput plotOutput downloadButton
#' @export
umapui <- function(id) {
  ns <- NS(id)  # Create a namespace function using the provided id
  fluidPage(  # Create a fluid page layout
    sidebarLayout(  # Layout with a sidebar and main panel
      sidebarPanel(  # Sidebar for inputs
        fileInput(ns("h5_files"), "Choose .h5 Files", multiple = TRUE, accept = ".h5"),  # Input for selecting .h5 files
        textInput(ns("project_names"), "Enter Group Names (comma-separated)", value = ""),  # Input for project/group names
        sliderInput(ns("resolution"), "Select Clustering Resolution", min = 0, max = 2, value = 0.8, step = 0.1),  # Slider for clustering resolution
        actionButton(ns("process"), "Process UMAP"),  # Button to start UMAP processing
        selectInput(ns("select_clusters"), "Select Clusters to Subset", choices = NULL, multiple = TRUE),  # Select clusters for subsetting
        actionButton(ns("subset_clusters"), "Subset Clusters")  # Button to subset clusters
      ),
      mainPanel(  # Main panel for output
        uiOutput(ns("helpTrajectoryInfo")),  # UI output for additional info (not used here)
        plotOutput(ns("umap"), width = "1280px", height = "840px"),  # UMAP plot output
        downloadButton(ns("download_seurat"), "Download Processed Seurat Object"),  # Download button for processed Seurat object
        downloadButton(ns("download_umap"), "Download UMAP Plot"),  # Download button for UMAP plot
        downloadButton(ns("download_subset_seurat"), "Download Subsetted Seurat Object")  # Download button for subsetted Seurat object
      )
    )
  )
}

#' Define Server Logic for the UMAP Module
#'
#' This function defines the server logic for the UMAP module in a Shiny app.
#' It handles file processing, UMAP computation, clustering, subsetting, and provides
#' download handlers for Seurat objects and UMAP plots.
#'
#' @param input Shiny input object.
#' @param output Shiny output object.
#' @param session Shiny session object.
#' @export
umapserver <- function(input, output, session) {
  seurat_list <- reactiveVal(list())  # Reactive value to store Seurat objects
  merged <- reactiveVal(NULL)  # Reactive value to store merged Seurat object
  umap <- reactiveVal(NULL)  # Reactive value to store UMAP plot
  subset <- reactiveVal(NULL)  # Reactive value to store subsetted Seurat object

  # Observe the event when the "Process UMAP" button is clicked
  observeEvent(input$process, {
    req(input$h5_files)  # Ensure .h5 files are uploaded
    req(input$project_names)  # Ensure project names are provided

    tryCatch({
      # Get the list of files and project names
      h5_files <- input$h5_files$datapath
      project_names <- strsplit(input$project_names, ",\\s*")[[1]]

      # Check if the number of files matches the number of project names
      if(length(h5_files) != length(project_names)) {
        showNotification("Number of project names must match the number of files", type = "error")
        return(NULL)
      }

      seurat_objs <- list()  # Initialize list to store Seurat objects
      for (i in seq_along(h5_files)) {
        # Read each .h5 file and create a Seurat object
        counts <- Read10X_h5(h5_files[i], use.names = TRUE, unique.features = TRUE)
        seurat_obj <- CreateSeuratObject(counts = counts, project = project_names[i])
        seurat_objs[[project_names[i]]] <- seurat_obj
      }
      seurat_list(seurat_objs)  # Store the list of Seurat objects

      print("Seurat Object is Created.")

      # Perform quality control on each Seurat object
      seurat_objs <- lapply(seurat_objs, perform_quality_control)
      print("QC done")

      # Create a vector of cell identifiers for each Seurat object
      cell_ids <- paste0("Sample", seq_along(seurat_objs))
      print(cell_ids)

      # Merge the Seurat objects into one
      merged_obj <- merge(
        x = seurat_objs[[1]],
        y = seurat_objs[-1],
        add.cell.ids = cell_ids,
        project = "TIgo",
        merge.data = TRUE
      )
      # Run normalization and other processing steps
      merged_obj <- NormalizeData(merged_obj) %>%
        FindVariableFeatures() %>%
        ScaleData() %>%
        RunPCA(verbose = TRUE) %>%
        RunHarmony(group.by.vars = "orig.ident", max.iter.harmony = 25, plot_convergence = FALSE) %>%
        FindNeighbors(reduction = "harmony", dims = 1:20) %>%
        FindClusters(resolution = input$resolution) %>%
        RunUMAP(reduction = "harmony", dims = 1:20)

      merged(merged_obj)  # Save the processed Seurat object
      print("Seurat Object is Processed.")

      # Generate UMAP plot
      umap_plot <- DimPlot(merged_obj, reduction = "umap", split.by = "orig.ident")

      # Get cluster levels for subsetting
      cluster_levels <- levels(Idents(merged_obj))

      # Update selectInput choices with available cluster levels
      updateSelectInput(session, "select_clusters", choices = cluster_levels, selected = cluster_levels)

      # Render UMAP plot
      output$umap <- renderPlot({
        print(umap_plot)
      })
    }, error = function(e) {
      # Show error notification if something goes wrong
      showNotification(paste("Error processing files:", e$message), type = "error")
    })
  })

  # Observe the event when the "Subset Clusters" button is clicked
  observeEvent(input$subset_clusters, {
    req(merged())  # Ensure the merged object is available
    req(input$select_clusters)  # Ensure clusters are selected for subsetting

    tryCatch({
      select_clusters <- input$select_clusters  # Get selected clusters

      # Use the cached merged object for subsetting
      merged_obj <- merged()
      Idents(merged_obj) <- "seurat_clusters"  # Ensure 'seurat_clusters' is the identity class

      # Manually filter the cells based on the selected clusters
      cells_to_keep <- WhichCells(merged_obj, idents = select_clusters)
      subsetted_obj <- merged_obj[, cells_to_keep]

      # Run additional processing steps on the subsetted object
      subsetted_obj <- ScaleData(subsetted_obj) %>%
        RunPCA(verbose = TRUE) %>%
        RunHarmony(group.by.vars = "orig.ident", max.iter.harmony = 25) %>%
        FindNeighbors(reduction = "harmony", dims = 1:20) %>%
        FindClusters() %>%
        RunUMAP(reduction = "harmony", dims = 1:20)

      subset(subsetted_obj)  # Save the subsetted Seurat object

      # Generate UMAP plot for the subsetted object
      umap_plot <- DimPlot(subsetted_obj, reduction = "umap")
      umap(umap_plot)  # Update the reactive UMAP plot

      # Render UMAP plot
      output$umap <- renderPlot({
        print(umap_plot)
      })
    }, error = function(e) {
      # Show error notification if something goes wrong
      showNotification(paste("Error subsetting clusters:", e$message), type = "error")
    })
  })

  # Download handler for the processed Seurat object
  output$download_seurat <- downloadHandler(
    filename = function() {
      paste("processed_seurat_", Sys.Date(), ".rds", sep = "")  # Define filename with current date
    },
    content = function(file) {
      tryCatch({
        saveRDS(merged(), file)  # Save the merged Seurat object
      }, error = function(e) {
        showNotification(paste("Error saving Seurat object:", e$message), type = "error")
      })
    }
  )

  # Download handler for the subsetted Seurat object
  output$download_subset_seurat <- downloadHandler(
    filename = function() {
      paste("subsetted_seurat_", Sys.Date(), ".rds", sep = "")  # Define filename with current date
    },
    content = function(file) {
      tryCatch({
        saveRDS(subset(), file)  # Save the subsetted Seurat object
      }, error = function(e) {
        showNotification(paste("Error saving subsetted Seurat object:", e$message), type = "error")
      })
    }
  )

  # Download handler for the UMAP plot
  output$download_umap <- downloadHandler(
    filename = function() {
      paste("umap_plot_", Sys.Date(), ".png", sep = "")  # Define filename with current date
    },
    content = function(file) {
      tryCatch({
        png(file)  # Start the PNG device
        print(umap())  # Print the UMAP plot
        dev.off()  # Close the device
      }, error = function(e) {
        showNotification(paste("Error saving UMAP plot:", e$message), type = "error")
      })
    }
  )
}

#' Add Percentage of Mitochondrial Genes
#'
#' Adds a new column to the Seurat object that contains the percentage of mitochondrial genes.
#' This is used for quality control purposes.
#'
#' @param seurat_obj A Seurat object.
#' @return The updated Seurat object with an additional column for mitochondrial gene percentage.
#' @export
add_mito_percentage <- function(seurat_obj) {
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")  # Add mitochondrial gene percentage to Seurat object
  return(seurat_obj)  # Return the updated Seurat object
}

#' Perform Quality Control on Seurat Object
#'
#' Applies quality control steps to the Seurat object, including adding mitochondrial gene percentages,
#' and filtering cells based on feature counts and mitochondrial gene percentage.
#'
#' @param seurat_obj A Seurat object.
#' @return The filtered Seurat object after applying quality control.
#' @export
perform_quality_control <- function(seurat_obj) {
  seurat_obj <- add_mito_percentage(seurat_obj)  # Add mitochondrial percentage

  # Plot QC metrics before subsetting
  VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

  # Subset cells based on QC criteria
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

  return(seurat_obj)  # Return the filtered Seurat object
}
