#' Trajectory Analysis Module UI
#'
#' This function defines the UI for the Trajectory Analysis module in a Shiny app.
#'
#' @param id A string representing the namespace ID for the module.
#'
#' @return A Shiny UI definition for the Trajectory Analysis module.
#' @export
#'
#' @importFrom shiny NS fluidPage sidebarLayout sidebarPanel mainPanel h4 fileInput selectInput actionButton textAreaInput numericInput
#' @importFrom shinyWidgets uiOutput verbatimTextOutput downloadButton plotlyOutput tableOutput fluidRow div
#' @importFrom plotly ggplotly
#'
TIUI <- function(id) {
    ns <- NS(id)  # Create a namespace function using the module's id to ensure unique IDs in the UI

    # Create the UI using a fluid page layout
    fluidPage(
        sidebarLayout(
            sidebarPanel(  # Sidebar panel contains input elements for the user
                h4("Trajectory Analysis"),  # A header for the sidebar panel

                # File input for uploading a Seurat object in RDS format
                fileInput(ns("seurat_object_file"), "Upload Seurat Object (RDS format)"),

                # Select input for choosing the trajectory inference method
                selectInput(ns("TI_method"), "Method of trajectory inference analysis",
                            choices = c("Monocle3", "slingshot")),

                # Button to trigger preprocessing of data
                actionButton(ns("runproprecessingdata"), "Preprocess Data"),

                # Placeholder for method-specific UI elements
                uiOutput(ns("method_specific_ui")),

                # Button to trigger trajectory inference analysis
                actionButton(ns("runanalysis"), "Run Trajectory Inference"),

                # Text area input for entering genes of interest, one per line
                textAreaInput(ns("interested_genes"), "Interested Genes (One per line)", value = ""),

                # Button to suggest genes based on analysis
                actionButton(ns("showGenes"), "Suggested Genes"),

                # Numeric input for setting the q-value threshold for differential expression analysis
                numericInput(
                    inputId = ns("qValueThreshold"),
                    label = "Q-value Threshold",
                    value = 0.05,
                    min = 0.00001,
                    max = 0.5,
                    step = 0.00001
                ),

                # Button to run differential expression analysis
                actionButton(ns("runDEG"), "Run DEG")
            ),
            mainPanel(  # Main panel contains output elements for displaying results
                uiOutput(ns("helpTrajectoryInfo")),  # Placeholder for additional trajectory information (if any)

                fluidRow(  # Organize the following elements in a fluid row
                    div(id = "overlap-container",  # A container for the plot and data table

                        # Container for the trajectory plot
                        div(id = "trajectoryPlot-container",
                            plotlyOutput(ns("trajectoryPlot"), height = "600px")  # Output for a Plotly plot
                        ),

                        # Container for the data table displaying analysis results
                        div(id = "dataTable-container",
                            tableOutput(ns("dataTable"))  # Output for a table
                        )
                    )
                ),

                # Button to display the next plot (useful if multiple plots are generated)
                actionButton(ns("nextPlot"), "Next Plot"),

                # Output for displaying trajectory-related textual information
                verbatimTextOutput(ns("trajectoryText"), placeholder = FALSE),

                # Button to download the trajectory plot as a file
                downloadButton(ns("downloadTrajectoryPlot"), "Download Plot"),

                # Button to download differential expression results as a file
                downloadButton(ns("downloadDEG"), "Download DEG File")
            )
        )
    )
}

#' Trajectory Analysis Module Server
#'
#' This function defines the server logic for the Trajectory Analysis module in a Shiny app.
#'
#' @param input,output,session Standard Shiny server function parameters.
#'
#' @return None. This function is used to handle reactive expressions and events in the module.
#' @export
#'
#' @importFrom shiny reactiveValues observeEvent req renderUI renderText renderPlotly renderTable showNotification downloadHandler
#' @importFrom SeuratObject readRDS
#' @importFrom tools file_ext
#' @importFrom zip zip
#' @importFrom ggplot2 ggsave
#' @importFrom dplyr select
#'

TIserver <- function(input, output, session) {
    ns <- session$ns
    # Initialize reactive values
    rv <- reactiveValues()

    rv$lastButtonPressed <- NULL

    # Render method-specific UI
    output$method_specific_ui <- renderUI({
        req(input$TI_method)  # Ensure TI_method is selected

        if (input$TI_method == "slingshot") {
            tagList(
                # Input for the start cluster ID
                textInput(ns("start_cluster"), "Start Cluster", value = "0"),
                # Input for the end cluster ID
                textInput(ns("end_cluster"), "End Cluster", value = ""),
                # Input for the cluster label (used in the Seurat object)
                textInput(ns("cluster_label"), "Cluster Label", value = ""),
                # Input for specifying a column name in metadata for coloring cells
                textInput(ns("color_cells"), "Color By", value = "")
            )
        } else if (input$TI_method == "Monocle3") {
            tagList(
                # Input for the start cluster ID
                textInput(ns("start_cluster"), "Start Cluster", value = "0"),
                # Input for the cluster label (used in the Seurat object)
                textInput(ns("cluster_label"), "Cluster Label", value = ""),
                # Input for specifying a column name in metadata for coloring cells
                textInput(ns("color_cells"), "Color By", value = ""),
                # Input for grouping cells in Monocle3 analysis
                textInput(ns("group_input"), "Group By", value = ""),
                # Checkbox to indicate if leaves should be labeled on the trajectory plot
                checkboxInput(ns("label_leaves_input"), "Label leaves:", value = FALSE),
                # Checkbox to indicate if branch points should be labeled on the trajectory plot
                checkboxInput(ns("label_branch_points_input"), "Label branch points:", value = TRUE)
            )
        }
    })


    observeEvent(input$runproprecessingdata, {
        req(input$TI_method)  # Ensure TI_method is selected

        print("Preprocessing Start")

        if (input$TI_method == "slingshot") {
            if (!is.null(input$seurat_object_file)) {
                # Load Seurat object
                seurat_object <- readRDS(input$seurat_object_file$datapath)

                # Access metadata
                metadata <- seurat_object@meta.data

                # Access gene expression data
                gene_expression <- LayerData(seurat_object, assay = "RNA", layer = "counts")

                # Convert to matrix if needed
                gene_expression <- as.matrix(gene_expression)

                # Store in reactive values
                rv$gene_expression <- gene_expression
                rv$cell_metadata <- metadata
                print("Slingshot data preprocessed")
            } else {
                stop("For slingshot method, please provide a Seurat object file.")
            }

        } else if (input$TI_method == "Monocle3") {
            if (!is.null(input$seurat_object_file)) {
                seurat_object <- readRDS(input$seurat_object_file$datapath)
                processed_data <- processInputData(seurat_object = seurat_object)

                # Store results in reactiveValues
                rv$gene_expression <- processed_data$gene_expression
                rv$cell_metadata <- processed_data$cell_metadata
                rv$gene_annotation <- processed_data$gene_annotation

                print("Monocle3 data preprocessed")
            } else {
                stop("For Monocle3 method, please provide a Seurat object.")
            }
        }

        # Update the output text after processing
        output$trajectoryText <- renderText({
            req(rv$cell_metadata)  # Ensure cell_metadata is available
            colnames_text <- paste(colnames(rv$cell_metadata), collapse = ", ")
            paste("Column names of metadata are:", colnames_text)
        })
    })

    observeEvent(input$runanalysis, {
        # Ensure that TI_method input is available
        req(input$TI_method)

        # Reset output variables
        rv$trajectoryPlot <- NULL
        rv$output_data <- NULL
        rv$results <- NULL

        if (input$TI_method == "slingshot") {
            print("Starting Slingshot trajectory analysis...")

            # Run the Slingshot trajectory analysis and handle any errors
            rv$output_data <- tryCatch({
                runslingshotTrajectory(
                    seurat_object_file = input$seurat_object_file,  # Input Seurat object file
                    start_cluster = input$start_cluster,            # Starting cluster for trajectory
                    end_cluster = input$end_cluster,                # Ending cluster for trajectory
                    seurat_cluster = input$cluster_label,           # Cluster label from Seurat object
                    color_cells = input$color_cells,                # Optional cell color mapping
                    session = session                               # Shiny session object
                )
            }, error = function(e) {
                print(paste("Error during trajectory analysis:", e$message))
                showNotification(paste("Error during trajectory analysis:", e$message), type = "error")
                return(NULL)
            })

            # Check if the output data is not NULL and assign the trajectory plot
            if (!is.null(rv$output_data)) {
                rv$trajectoryPlot <- rv$output_data[[2]]

                # Render the trajectory plot in the UI
                output$trajectoryPlot <- renderPlotly({
                    req(rv$trajectoryPlot)  # Ensure the plot is available
                    ggplotly(rv$trajectoryPlot)
                })

                # Define the download handler for the trajectory plot
                output$downloadTrajectoryPlot <- downloadHandler(
                    filename = function() {
                        "trajectory_plot_slingshot.png"  # Filename for the downloaded plot
                    },
                    content = function(file) {
                        print("Downloading trajectory plot")
                        ggsave(
                            file,
                            plot = rv$trajectoryPlot,  # Plot to save
                            device = "png",              # File format
                            width = 16,                  # Plot width
                            height = 9,                  # Plot height
                            units = "in",                # Units for dimensions
                            dpi = 300                    # Resolution
                        )
                    }
                )
            } else {
                # Clear plot output if data is NULL
                output$trajectoryPlot <- renderPlotly(NULL)
                output$downloadTrajectoryPlot <- downloadHandler(NULL)
            }

        } else if (input$TI_method == "Monocle3") {
            print("Starting Monocle3 trajectory analysis...")

            # Run the Monocle3 trajectory analysis and handle any errors
            rv$results <- tryCatch({
                runmonocle3Trajectory(
                    gene_expression = rv$gene_expression,  # Gene expression data
                    cell_metadata = rv$cell_metadata,      # Cell metadata
                    gene_annotation = rv$gene_annotation,  # Gene annotation data
                    seurat_cluster = input$cluster_label,  # Cluster label for trajectory
                    start_cluster = input$start_cluster,   # Starting cluster for trajectory
                    label_groups_by_cluster = FALSE,       # Whether to label groups by cluster
                    color_cells_by = input$color_cells,    # Optional cell color mapping
                    label_leaves_input = input$label_leaves_input,  # Input for labeling leaves
                    label_branch_points_input = input$label_branch_points_input,  # Input for branch points
                    group_input = input$group_input,       # Grouping input
                    session = session                     # Shiny session object
                )
            }, error = function(e) {
                print(paste("Error during trajectory analysis:", e$message))
                showNotification(paste("Error during trajectory analysis:", e$message), type = "error")
                return(NULL)
            })

            # Check if the results are not NULL and assign the trajectory plot
            if (!is.null(rv$results)) {
                rv$trajectoryPlot <- rv$results[[2]]

                # Render the trajectory plot in the UI
                output$trajectoryPlot <- renderPlotly({
                    req(rv$trajectoryPlot)  # Ensure the plot is available
                    ggplotly(rv$trajectoryPlot)
                })

                # Define the download handler for the trajectory plot
                output$downloadTrajectoryPlot <- downloadHandler(
                    filename = function() {
                        "trajectory_plot_monocle.png"  # Filename for the downloaded plot
                    },
                    content = function(file) {
                        print("Downloading trajectory plot")
                        ggsave(
                            file,
                            plot = rv$trajectoryPlot,  # Plot to save
                            device = "png",         # File format
                            width = 16,             # Plot width
                            height = 9,             # Plot height
                            units = "in",           # Units for dimensions
                            dpi = 300               # Resolution
                        )
                    }
                )
            } else {
                # Clear plot output if data is NULL
                output$trajectoryPlot <- renderPlotly(NULL)
                output$downloadTrajectoryPlot <- downloadHandler(NULL)
            }
        }

        print("Trajectory Analysis completed")
    })

    observeEvent(input$runDEG, {
        req(input$TI_method)

        # Reset output variables
        rv$DEG <- NULL

        print(paste("TI_method selected:", input$TI_method))

        if (input$TI_method == "slingshot") {
            print("DEG analysis started for Slingshot.")
            req(rv$output_data)
            req(input$seurat_object_file)

            # Initialize and check file extension
            inFile <- input$seurat_object_file
            ext <- tools::file_ext(inFile$datapath)

            # Validate the file type and read the file accordingly
            seurat_object <- tryCatch({
                if (ext == "rds") {
                    readRDS(inFile$datapath)  # Read RDS file
                } else {
                    stop("Invalid file type. Please upload an RDS file.")
                }
            }, error = function(e) {
                print(paste("Error reading file:", e$message))
                showNotification(paste("Error reading file:", e$message), type = "error")
                return(NULL)
            })

            # Check if the loaded object is a Seurat object
            if (!inherits(seurat_object, "Seurat")) {
                print("The uploaded RDS file does not contain a Seurat object.")
                showNotification("The uploaded RDS file does not contain a Seurat object.", type = "error")
                return(NULL)
            }

            # Perform DEG analysis
            rv$DEG <- tryCatch({
                runDEGtrajectory_S(
                    output_slingshot = rv$output_data,
                    Tcells = seurat_object,  # Pass the Seurat object
                    q = input$q_value_threshold,
                    session = session
                )
            }, error = function(e) {
                print(paste("Error during DEG analysis:", e$message))
                showNotification(paste("Error during DEG analysis:", e$message), type = "error")
                return(NULL)
            })

            # Set up download handler if DEG analysis is successful
            if (!is.null(rv$DEG)) {
                output$downloadDEG <- downloadHandler(
                    filename = function() {
                        paste("DEG_results_", Sys.Date(), ".zip", sep = "")
                    },
                    content = function(file) {
                        tempDir <- tempdir()  # Get the path to the temp directory

                        # Create temporary files
                        deg_csv <- file.path(tempDir, "DEG_table.csv")
                        plot_png1 <- file.path(tempDir, "TrajectoryPlot1.png")


                        # Save files to tempDir
                        if (!is.null(rv$DEG[[1]])) {
                            write.csv(rv$DEG[[1]], deg_csv)
                        }

                        if (!is.null(rv$DEG[[2]])) {
                            ggsave(plot_png1, plot = rv$DEG[[2]])
                        }

                        # Zip files
                        zip::zip(zipfile = file, files = c(deg_csv, plot_png1))
                    },
                    contentType = "application/zip"
                )
            }

        } else if (input$TI_method == "Monocle3") {
            print("DEG analysis started for Monocle3.")
            req(rv$results)

            rv$DEG <- tryCatch({
                runDEGtrajectory_M(
                    cds = rv$results[[1]],
                    interested_genes = strsplit(input$interested_genes, "\n", fixed = TRUE)[[1]],
                    seurat_cluster = input$cluster_label,
                    start_cluster = input$start_cluster,
                    q = input$qValueThreshold,
                    session = session
                )
            }, error = function(e) {
                print(paste("Error during DEG analysis:", e$message))
                showNotification(paste("Error during DEG analysis:", e$message), type = "error")
                return(NULL)
            })

            # Set up download handler if DEG analysis is successful
            if (!is.null(rv$DEG)) {
                output$downloadDEG <- downloadHandler(
                    filename = function() {
                        paste("DEG_results_", Sys.Date(), ".zip", sep = "")
                    },
                    content = function(file) {
                        tempDir <- tempdir()  # Get the path to the temp directory

                        # Create temporary files
                        deg_csv <- file.path(tempDir, "DEG_table.csv")
                        plot_png1 <- file.path(tempDir, "TrajectoryPlot1.png")
                        plot_png2 <- file.path(tempDir, "TrajectoryPlot2.png")

                        # Save files to tempDir
                        if (!is.null(rv$DEG[[1]])) {
                            write.csv(rv$DEG[[1]], deg_csv)
                        }

                        if (!is.null(rv$DEG[[2]])) {
                            ggsave(plot_png1, plot = rv$DEG[[2]])
                        }

                        if (!is.null(rv$DEG[[3]])) {
                            ggsave(plot_png2, plot = rv$DEG[[3]])
                        }

                        # Zip files
                        zip::zip(zipfile = file, files = c(deg_csv, plot_png1, plot_png2))
                    },
                    contentType = "application/zip"
                )
            }
        }

        # If DEG analysis fails, clear the download handler
        if (is.null(rv$DEG)) {
            output$downloadDEG <- downloadHandler(NULL, NULL)
        }
    })



    observeEvent(input$nextPlot, {
        if (!is.null(rv$DEG)) {
            rv$currentPlotIndex <- rv$currentPlotIndex %% 2 + 1  # Toggle plot index between 1 and 2
            print(paste("Plot index changed to:", rv$currentPlotIndex))  # Debug message
        }
    })

    # Render the table when the button is clicked
    observeEvent(input$showGenes, {
        output$dataTable <- renderTable({
            fixedTable
        }, striped = TRUE, bordered = TRUE, hover = TRUE)  # Optional styling
    })
}

#' Process Input Data for Trajectory Analysis
#'
#' This function processes the input data for trajectory analysis. It accepts either a Seurat object or raw data and extracts the gene expression matrix, cell metadata, and gene annotations.
#'
#' @param gene_expression A gene expression matrix (default is NULL).
#' @param cell_metadata A data frame containing cell metadata (default is NULL).
#' @param gene_annotation A data frame containing gene annotations (default is NULL).
#' @param seurat_object A Seurat object containing RNA assay data (default is NULL).
#'
#' @return A list containing gene expression, cell metadata, and gene annotation.
#' @export
processInputData <- function(gene_expression = NULL, cell_metadata = NULL, gene_annotation = NULL, seurat_object = NULL) {
    if (!is.null(seurat_object)) {
        if (!inherits(seurat_object, "Seurat")) {
            stop("The provided object is not a Seurat object.")
        }
        #generate gene expression data
        gene_expression <- LayerData(seurat_object, assay = "RNA", layer = "counts")
        #generate metadata and subset
        cell_metadata <- seurat_object@meta.data
        cells_in_expression <- colnames(gene_expression)
        cell_metadata <- cell_metadata[cells_in_expression, , drop = FALSE]

        gene_annotation <- data.frame(
            gene_short_name = rownames(gene_expression),
            row.names = rownames(gene_expression)
        )
    }
    return(list(gene_expression = gene_expression, cell_metadata = cell_metadata, gene_annotation = gene_annotation))
}

#' Run Slingshot Trajectory Inference
#'
#' This function runs the Slingshot trajectory inference on a Seurat object. It generates pseudotime values and integrates them into the Seurat metadata.
#'
#' @param seurat_object_file The path to the Seurat object file.
#' @param start_cluster The starting cluster for the trajectory.
#' @param end_cluster The ending cluster for the trajectory.
#' @param seurat_cluster The name of the cluster column in the Seurat object.
#' @param color_cells A vector of colors for cells.
#' @param session The Shiny session object.
#'
#' @return A list containing the Slingshot output, plot, reduced dimension data, and pseudotime values.
#' @export
runslingshotTrajectory <- function(seurat_object_file, start_cluster, end_cluster, seurat_cluster, color_cells, session) {
    tryCatch({
        checkpoints <- list(method_afterpreproc = as.numeric(Sys.time()))

        print("Running Slingshot on the reduced dimension data and clustering...")
        seurat_object <- readRDS(seurat_object_file$datapath)
        rd <- seurat_object@reductions$umap@cell.embeddings
        cl <- as.factor(seurat_object@meta.data[[seurat_cluster]])

        # Verify dimensions
        if (nrow(rd) != length(cl)) {
            stop("The length of clusterLabels does not match the number of rows in the reduced dimension data.")
        }

        sds <- slingshot::slingshot(
            rd,
            clusterLabels = cl,
            start.clus = start_cluster,
            end.clus = end_cluster,
            reducedDim = 'UMAP',
            shrink = 1L,
            reweight = TRUE,
            reassign = TRUE,
            maxit = 10L,
            smoother = "smooth.spline"
        )

        print("Integrating pseudotime into metadata...")
        metadata <- seurat_object@meta.data
        pseudotime <- slingPseudotime(sds, na = FALSE)
        metadata[["pseudotime"]] <- rowMeans(pseudotime, na.rm = TRUE)

        print("Generating trajectory information...")
        lineages <- slingLineages(sds)
        lineage_ctrl <- slingParams(sds)
        cluster_network <- lineages %>%
            map_df(~ tibble(from = .[-length(.)], to = .[-1])) %>%
            unique() %>%
            mutate(
                length = lineage_ctrl$dist[cbind(from, to)],
                directed = TRUE,
                length = 1
            )

        dimred <- slingshot::slingReducedDim(sds)
        cluster <- slingClusterLabels(sds)
        lin_assign <- apply(slingCurveWeights(sds), 1, which.max)

        progressions <- map_df(seq_along(lineages), function(l) {
            ind <- lin_assign == l
            lin <- lineages[[l]]
            pst.full <- slingPseudotime(sds, na = FALSE)[, l]
            pst <- pst.full[ind]
            means <- sapply(lin, function(clID) {
                stats::weighted.mean(pst.full, cluster[, clID])
            })
            # Ensure breaks are unique
            non_ends <- unique(means[-c(1, length(means))])
            unique_breaks <- unique(c(-Inf, non_ends, Inf))
            edgeID.l <- as.numeric(cut(pst, breaks = unique_breaks))
            from.l <- lineages[[l]][edgeID.l]
            to.l <- lineages[[l]][edgeID.l + 1]
            m.from <- means[from.l]
            m.to <- means[to.l]

            pct <- (pst - m.from) / (m.to - m.from)
            pct[pct < 0] <- 0
            pct[pct > 1] <- 1

            tibble(cell_id = names(which(ind)), from = from.l, to = to.l, percentage = pct)
        })

        output_slingshot <-
            dynwrap::wrap_data(
                cell_ids = rownames(metadata)
            ) %>%
            dynwrap::add_trajectory(
                milestone_network = cluster_network,
                progressions = progressions,
                lineages = lineages
            ) %>%
            dynwrap::add_dimred(
                dimred = dimred
            ) %>%
            dynwrap::add_timings(checkpoints)

        plot <- dynplot::plot_dimred(
            output_slingshot,
            color_cells = "pseudotime"
        )

        output_slingshot = list(output_slingshot, plot, rd, metadata[["pseudotime"]])

        print("Trajectory analysis and data preparation completed.")
        return(output_slingshot)
    }, error = function(e) {
        print(paste("Error in runslingshotTrajectory:", e))
        sendSweetAlert(
            session = session,
            title = "Trajectory Inference Error",
            text = paste("An error occurred during trajectory inference:", e),
            type = "error"
        )
        return(NULL)
    })
}

#' Run Monocle3 Trajectory Inference
#'
#' This function runs the Monocle3 trajectory inference on processed input data. It creates a cell data set (CDS) and performs preprocessing, dimension reduction, clustering, and graph learning.
#'
#' @param gene_expression A gene expression matrix.
#' @param cell_metadata A data frame containing cell metadata.
#' @param gene_annotation A data frame containing gene annotations.
#' @param seurat_cluster The name of the cluster column in the Seurat object.
#' @param start_cluster The starting cluster for the trajectory.
#' @param label_groups_by_cluster A logical value indicating whether to label groups by cluster.
#' @param color_cells_by The method used to color cells in the plot.
#' @param label_branch_points_input A logical value indicating whether to label branch points.
#' @param session The Shiny session object.
#' @param label_leaves_input A logical value indicating whether to label leaves.
#' @param group_input A grouping variable for alignment.
#'
#' @return A list containing the Monocle3 CDS object and a trajectory plot.
#' @export
runmonocle3Trajectory <- function(gene_expression, cell_metadata, gene_annotation, seurat_cluster, start_cluster, label_groups_by_cluster, color_cells_by, label_branch_points_input, session, label_leaves_input, group_input) {
    tryCatch({
        if (!inherits(gene_expression, "dgCMatrix")) {
            library(Matrix)
            gene_expression <- as(gene_expression, "dgCMatrix")
        }

        if (!is.data.frame(cell_metadata)) {
            cell_metadata <- as.data.frame(cell_metadata)
        }

        if (!is.data.frame(gene_annotation)) {
            gene_annotation <- as.data.frame(gene_annotation)
        }

        print("Running Monocle3 on the reduced dimension data and clustering...")
        print("gene_expression class:")
        print(class(gene_expression))
        print("cell_metadata class:")
        print(class(cell_metadata))
        print("gene_annotation class:")
        print(class(gene_annotation))
        print("group_input:")
        print(group_input)

        cds <- new_cell_data_set(gene_expression,
                                 cell_metadata = cell_metadata,
                                 gene_metadata = gene_annotation)
        print("seurat object created")
        cds <- preprocess_cds(cds, num_dim = 50)
        print("preprocessing done")
        cds <- align_cds(cds, alignment_group = group_input)
        print("align done")
        cds <- reduce_dimension(cds, reduction_method = "UMAP")
        print("reduction done")
        cds <- cluster_cells(cds)
        print("clustering done")
        cds <- learn_graph(cds)
        print("lineage done")

        # Check if get_earliest_principal_node returns a valid value
        root_nodes <- get_earliest_principal_node(cds, seurat_cluster, start_cluster)
        print("root_nodes:")
        print(root_nodes)
        if (is.null(root_nodes)) {
            stop("get_earliest_principal_node returned NULL. Please check the function and its inputs.")
        }

        cds <- order_cells(cds, root_pr_nodes = root_nodes)
        print("Ordering of cells completed")

        # Visualize the trajectory inference
        plot <- plot_cells(cds,
                           color_cells_by = color_cells_by,
                           label_leaves = label_leaves_input,
                           label_branch_points = label_branch_points_input,
                           label_groups_by_cluster = TRUE) +
            theme(
                text = element_text(size = 20), # General text size
                plot.title = element_text(size = 20, face = "bold"), # Title text size
                legend.title = element_text(size = 15), # Legend title text size
                legend.text = element_text(size = 10), # Legend text size
                axis.text = element_text(size = 10), # Axis text size
                axis.title = element_text(size = 13),# Axis title text size
                strip.text = element_text(size = 13)  # Adjust font size for group labels
            )

        out <- list(cds, plot)
        return(out)
    }, error = function(e) {
        print(paste("Error in runmonocle3Trajectory:", e$message))
        sendSweetAlert(
            session = session,
            title = "Trajectory Inference Error",
            text = paste("An error occurred during trajectory inference:", e$message),
            type = "error"
        )
        return(NULL)
    })
}

#' Run Differential Expression Analysis Using Monocle3
#'
#' This function performs differential expression analysis along a trajectory using Monocle3. It identifies differentially expressed genes (DEGs) and generates a heatmap and pseudotime plot.
#'
#' @param cds A Monocle3 cell data set (CDS) object.
#' @param interested_genes A vector of gene names of interest.
#' @param seurat_cluster The name of the cluster column in the Seurat object.
#' @param start_cluster The starting cluster for the trajectory.
#' @param q The threshold for the q-value to identify significant genes.
#' @param session The Shiny session object.
#'
#' @return A list containing the aggregated expression matrix, heatmap, pseudotime plot, and DEG results.
#' @export
runDEGtrajectory_M <- function(cds, interested_genes, seurat_cluster, start_cluster, q, session) {
    tryCatch({
        # Perform DEG analysis
        cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=4)
        print("calculation done")

        results_df <- as.data.frame(cds_pr_test_res)

        pr_deg_ids <- row.names(subset(cds_pr_test_res, q_value < q))  # Use user-defined q_value


        # Process for heatmap
        gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=c(10^seq(-6, -1)))
        cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), cell_group=colData(cds)[[seurat_cluster]])
        agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
        row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))

        # Create heatmap
        plot1 <- pheatmap::pheatmap(agg_mat, cluster_rows = TRUE, cluster_cols = TRUE, scale = "column", clustering_method = "ward.D2", fontsize = 6)

        # Process for pseudotime plot
        interested_lineage_cds <- cds[rowData(cds)$gene_short_name %in% interested_genes,]
        interested_lineage_cds <- order_cells(interested_lineage_cds, root_pr_nodes=get_earliest_principal_node(interested_lineage_cds, seurat_cluster, start_cluster))
        plot2 <- plot_genes_in_pseudotime(interested_lineage_cds, color_cells_by="pseudotime", min_expr=0.5)

        output <- list(agg_mat, plot1, plot2, results_df)
        return(output)
    }, error = function(e) {
        print(paste("Error in runDEGtrajectory:", e))
        sendSweetAlert(session = session, title = "DEG Error", text = paste("An error occurred during DEG analysis:", e), type = "error")
        return(NULL)
    })
}

#' Run Differential Expression Analysis Using Slingshot
#'
#' This function performs differential expression analysis along a trajectory using Slingshot. It fits a GAM model for each gene and generates a heatmap of the expression patterns.
#'
#' @param output_slingshot The output from the Slingshot trajectory analysis.
#' @param Tcells A Seurat object containing RNA assay data.
#' @param q The threshold for the q-value to identify significant genes.
#' @param session The Shiny session object.
#'
#' @return A list containing the DEG results, heatmap, and smooth terms.
#' @export
runDEGtrajectory_S <- function(output_slingshot, Tcells, q, session) {
    tryCatch({
        # Extract data from output_slingshot
        metadata <- output_slingshot$metadata
        dimred <- output_slingshot$dimred
        milestone_network <- output_slingshot$milestone_network
        progressions <- output_slingshot$progressions

        # Extract raw counts
        counts_matrix <- LayerData(Tcells, assay = "RNA", layer = "counts")


        # Convert counts_matrix to Matrix if needed
        counts_matrix <- as(counts_matrix, "Matrix")

        print("DATA prprocessed")

        # Ensure column names match cell IDs
        cell_ids <- colnames(counts_matrix)

        # Ensure pseudotime is named correctly
        cell_ids_pseudotime <- names(pseudotime)
        common_cell_ids <- intersect(cell_ids, cell_ids_pseudotime)


        # Subset pseudotime to match cell IDs in counts matrix
        pseudotime_aligned <- pseudotime[cell_ids]

        # Extract expression matrix and pseudotime
        expression_long <- as.data.frame(counts_matrix) %>%
            rownames_to_column("gene") %>%
            pivot_longer(cols = -gene, names_to = "cell_id", values_to = "expression")

        expression_long$pseudotime <- pseudotime_aligned[expression_long$cell_id]

        # Fit a GAM model for each gene and extract smooth terms
        print("Calculation started")
        results <- expression_long %>%
            group_by(gene) %>%
            do({
                model <- tryCatch({
                    gam(expression ~ s(pseudotime), data = .)
                }, error = function(e) {
                    warning(sprintf("GAM model fitting failed for gene %s: %s", unique(.$gene), e$message))
                    return(NULL)
                })
                if (!is.null(model)) {
                    tidy(model) %>%
                        filter(term == "s(pseudotime)") %>%
                        mutate(gene = unique(.$gene))
                } else {
                    data.frame(gene = unique(.$gene), term = NA, p.value = NA, q_value = NA)
                }
            })

        # Adjust p-values for multiple testing (e.g., using the Benjamini-Hochberg method)
        results <- results %>%
            mutate(q_value = p.adjust(p.value, method = "BH"))

        print("filtering.....")
        # Filter genes with q-value <= q_value_threshold
        significant_genes <- results %>%
            filter(q_value < q) %>%
            pull(gene)

        if (length(significant_genes) == 0) {
            stop("No significant genes found with the given q_value_threshold.")
        }

        # Extract smooth terms for each gene
        smooth_terms <- expression_long %>%
            filter(gene %in% significant_genes) %>%
            group_by(gene) %>%
            do({
                model <- tryCatch({
                    gam(expression ~ s(pseudotime), data = .)
                }, error = function(e) {
                    warning(sprintf("GAM model fitting failed for gene %s: %s", unique(.$gene), e$message))
                    return(NULL)
                })
                if (!is.null(model)) {
                    data.frame(pseudotime = seq(min(.$pseudotime), max(.$pseudotime), length.out = 100),
                               fitted = predict(model, newdata = data.frame(pseudotime = seq(min(.$pseudotime), max(.$pseudotime), length.out = 100))))
                } else {
                    data.frame(pseudotime = seq(min(.$pseudotime), max(.$pseudotime), length.out = 100),
                               fitted = NA)
                }
            }) %>%
            ungroup()

        # Reshape data for heatmap
        heatmap_data <- smooth_terms %>%
            pivot_wider(names_from = pseudotime, values_from = fitted) %>%
            column_to_rownames(var = "gene")

        # Create a data frame for column annotations with numeric pseudotime values
        pseudotime_values <- as.numeric(colnames(heatmap_data))  # Convert column names to numeric values
        annotation_df <- data.frame(pseudotime = pseudotime_values)
        rownames(annotation_df) <- colnames(heatmap_data)

        # Create a custom color scale based on the pseudotime range
        pseudotime_breaks <- seq(min(pseudotime_values), max(pseudotime_values), length.out = 100)
        pseudotime_colors <- colorRampPalette(c("yellow", "white", "green"))(100)
        pseudotime_color_map <- setNames(pseudotime_colors, pseudotime_breaks)

        print("Drawing heatmap.....")
        # Create the heatmap
        heatmap <- pheatmap(
            as.matrix(heatmap_data),
            cluster_rows = TRUE,
            cluster_cols = TRUE,
            scale = "row",
            color = colorRampPalette(c("blue", "white", "red"))(100),
            annotation_col = annotation_df,
            annotation_colors = list(
                pseudotime = pseudotime_color_map
            ),
            show_rownames = FALSE,
            show_colnames = FALSE
        )
        print("DEG analysis completed.")

        result <- list(results = results, heatmap = heatmap, smooth_terms = smooth_terms, heatmap_data = heatmap_data)

    }, error = function(e) {
        print(paste("Error in runDEGtrajectory:", e$message))
        sendSweetAlert(session = session, title = "DEG Error", text = paste("An error occurred during DEG analysis:", e$message), type = "error")
        return(NULL)
    })
}

#' Get Earliest Principal Node
#'
#' This function identifies the earliest principal node for ordering cells in a Monocle3 trajectory. It calculates the root principal node based on the start cluster.
#'
#' @param cds A Monocle3 cell data set (CDS) object.
#' @param seurat_cluster The name of the cluster column in the Seurat object.
#' @param start_cluster The starting cluster for the trajectory.
#'
#' @return The earliest principal node in the trajectory graph.
#' @export
get_earliest_principal_node <- function(cds, seurat_cluster, start_cluster){
    # Identify the cell IDs that belong to the start cluster
    cell_ids <- which(colData(cds)[[seurat_cluster]] == start_cluster)
    if (length(cell_ids) == 0) {
        stop("No cells found for the specified start cluster.")
    }
    print(paste("Number of cells in start cluster:", length(cell_ids)))

    # Get the closest vertex matrix
    closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
    if (is.null(closest_vertex)) {
        stop("closest_vertex is NULL.")
    }
    closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
    print(dim(closest_vertex))

    # Calculate the root principal nodes
    vertex_table <- table(closest_vertex[cell_ids,])
    if (length(vertex_table) == 0) {
        stop("Vertex table is empty.")
    }
    print(vertex_table)
    root_pr_node <- names(which.max(vertex_table))
    if (length(root_pr_node) == 0) {
        stop("No root principal node found.")
    }
    root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(root_pr_node)]

    root_pr_nodes
}

# Define column data
isg <- c("ADAR", "B2M", "BATF2", "BST2", "C1S", "CASP1", "CASP8", "CCRL2",
         "CD47", "CD74", "CMPK2", "CMTR1", "CNP", "CSF1", "CXCL10", "CXCL11", "DDX60",
         "DHX58", "EIF2AK2", "ELF1", "EPSTI1", "GBP2", "GBP4", "GMPR", "HELZ2", "HERC6",
         "HLA-C", "IFI27", "IFI30", "IFI35", "IFI44", "IFI44L", "IFIH1", "IFIT2", "IFIT3",
         "IFITM1", "IFITM2", "IFITM3", "IL15", "IL4R", "IL7", "IRF1", "IRF2", "IRF7",
         "IRF9", "ISG15", "ISG20", "LAMP3", "LAP3", "LGALS3BP", "LPAR6", "LY6E", "MOV10",
         "MVB12A", "MX1", "NCOA7", "NMI", "NUB1", "OAS1", "OASL", "OGFR", "PARP12", "PARP14",
         "PARP9", "PLSCR1", "PNPT1", "PROCR", "PSMA3", "PSMB8", "PSMB9", "PSME1", "PSME2",
         "RIPK2", "RNF31", "RSAD2", "RTP4", "SAMD9", "SAMD9L", "SELL", "SLC25A28", "SP110",
         "STAT2", "TAP1", "TDRD7", "TENT5A", "TMEM140", "TRAFD1", "TRIM14", "TRIM21",
         "TRIM25", "TRIM26", "TRIM5", "TXNIP", "UBA7", "UBE2L6", "USP18", "WARS1")

nk_cytokine_expr <- c("CCL3", "CCL4", "IFNG", "TNF")

t_cell_exhaustion <- c(
    "KIAA1671", "ABI3", "ADAM19", "AKNA", "ARHGAP9", "ARL6IP1", "ARMC7", "BCL2A1",
    "CBX4", "CCL3", "CCL4", "CCL5", "CD160", "CD164", "CD27", "CD3E", "CD3G", "CD7",
    "CD82", "CD8A", "CST7", "CXCR6", "DAPK2", "DTX1", "DUSP2", "EFHD2", "EIF4A2",
    "FAM189B", "FASLG", "FOXN3", "FYN", "GIMAP1", "GIMAP6", "GIMAP7", "GLRX",
    "GNG2", "GRAMD1A", "GZMA", "GZMB", "GZMK", "HLA-A", "HCST", "HSPA5", "ID2",
    "IFIT3", "IL21R", "ISG15", "ITK", "ITPKB", "LAG3", "LAX1", "LRRK1", "MBNL1",
    "MXD4", "NR4A2", "PDCD1", "PFDN5", "PLA2G16", "PLAC8", "PRDX5", "PRKCH",
    "PSMB10", "PSMB8", "PSME1", "PTGER4", "PTPN18", "PTPN22", "RGS1", "RGS2",
    "RGS3", "RTP4", "RUNX3", "SH2D2A", "SHISA5", "SIPA1", "SLC3A2", "STAT1",
    "STK17B", "TAP2", "TAPBP", "TAPBPL", "TNFRSF1B", "TOX", "UCP2", "VMP1", "ZBP1"
)

hla_class_ii <- c(
    "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1",
    "HLA-DQA2", "HLA-DQB1", "HLA-DQB1-AS1", "HLA-DQB2", "HLA-DRA", "HLA-DRB1",
    "HLA-DRB5"
)

nk_cell_exhaustion <- c("PDCD1","LAG3","HAVCR2")

detected_cytokines <- c(
    "CSF1", "CSF2", "CSF3", "IFNA1", "IFNA13", "IFNA17", "IFNA2", "IFNA21",
    "IFNA4", "IFNA5", "IFNA6", "IFNA7", "IFNB1", "IFNE", "IFNG", "IFNK", "IFNL1",
    "IFNL2", "IFNL3", "IFNW1", "IFNA10", "IFNA14", "IFNA16", "IFNA8", "IL10",
    "IL11", "IL12A", "IL12B", "IL13", "IL15", "IL16", "IL17A", "IL17B", "IL17C",
    "IL17D", "IL17F", "IL18", "IL19", "IL1A", "IL1B", "IL1F10", "IL2", "IL20",
    "IL22", "IL23A", "IL24", "IL25", "IL26", "IL27", "IL3", "IL31", "IL32", "IL33",
    "IL34", "IL36B", "IL36G", "IL37", "IL4", "IL5", "IL6", "IL7", "IL8", "IL9",
    "IL21", "IL36A", "TNFSF10", "TNFSF11", "TNFSF12", "TNFSF13", "TNFSF13B", "TNFSF14",
    "TNFSF15", "TNFSF18", "TNFSF4", "TNFSF8", "TNFSF9", "CXCL1", "CXCL10", "CXCL11",
    "CXCL12", "CXCL13", "CXCL14", "CXCL16", "CXCL17", "CXCL2", "CXCL3", "CXCL5",
    "CXCL6", "CXCL9", "CCL11", "CCL13", "CCL14", "CCL15", "CCL16", "CCL17", "CCL18",
    "CCL19", "CCL2", "CCL20", "CCL21", "CCL22", "CCL23", "CCL24", "CCL25", "CCL26",
    "CCL28", "CCL3", "CCL3L1", "CCL3L3", "CCL4", "CCL4L1", "CCL4L2", "CCL5", "CCL7",
    "CCL8", "CCL1", "TSLP", "LIF", "OSM", "TNF", "LTA", "LTB", "CD40L", "FASL",
    "CD70", "TGFB1", "MIF", "CX3CL1"
)

# Define maximum length
max_length <- max(length(isg), length(nk_cytokine_expr), length(t_cell_exhaustion), length(hla_class_ii), length(nk_cell_exhaustion), length(detected_cytokines))

# Extend columns to the maximum length with empty strings or NA
fixedTable <- data.frame(
    ISG = c(isg, rep("", max_length - length(isg))),
    `NK or CD8+ T cell pro-inflammatory cytokine expression` = c(nk_cytokine_expr, rep("", max_length - length(nk_cytokine_expr))),
    `T cell exhaustion` = c(t_cell_exhaustion, rep("", max_length - length(t_cell_exhaustion))),
    `HLA Class II` = c(hla_class_ii, rep("", max_length - length(hla_class_ii))),
    `NK cell exhaustion` = c(nk_cell_exhaustion, rep("", max_length - length(nk_cell_exhaustion))),
    `Detected Cytokines` = c(detected_cytokines, rep("", max_length - length(detected_cytokines)))
)

