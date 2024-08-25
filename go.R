# Define the GO Enrichment module UI
goUI <- function(id) {
    ns <- NS(id)  # Create a namespace function using the module's id to ensure unique IDs in the UI

    # Create the UI using a fluid page layout
    fluidPage(
        titlePanel("GO Enrichment Analysis"),  # Title panel for the app

        sidebarLayout(  # Layout with a sidebar and main panel
            sidebarPanel(  # Sidebar contains input elements for the user
                # File input for uploading a gene list
                fileInput(ns("file1"), "Upload Gene List",
                          accept = c(".txt", ".csv")),  # Accepts text or CSV files

                # Text input for specifying the column name that contains the gene IDs
                textInput(ns("geneColumn"), "Gene Column Name", value = "GeneID"),

                # Button to trigger the GO Enrichment Analysis
                actionButton(ns("goButton"), "Run GO Enrichment Analysis")
            ),

            mainPanel(  # Main panel contains output elements for displaying results
                # Output for displaying the enrichment plot
                plotOutput(ns("enrichmentPlot"), width = "100%", height = "840px"),

                # Button to download the enrichment results as a CSV file
                downloadButton(ns("downloadResults"), "Download Enrichment Results"),

                # Button to download the enrichment plot as a PNG file
                downloadButton(ns("downloadPlot"), "Download Enrichment Plot")
            )
        )
    )
}



# Define the GO Enrichment module server logic
goserver <- function(input, output, session) {
     print("GO enrichment started.")
    # Reactive expression to process and store the uploaded gene list
    geneList <- reactive({
        req(input$file1)  # Ensure that a file is uploaded

        inFile <- input$file1  # Get the uploaded file
        ext <- tools::file_ext(inFile$datapath)  # Get the file extension

        # Read the file based on its extension
        if (ext == "csv") {
            data <- read.csv(inFile$datapath)  # Read CSV file
        } else if (ext == "txt") {
            data <- read.table(inFile$datapath, header = TRUE)  # Read text file with a header
        } else {
            stop("Invalid file type")  # Stop if the file type is not supported
        }

        # Check if the specified gene column exists in the data
        if (!(input$geneColumn %in% colnames(data))) {
            stop("Gene column name does not match any column in the data")
        }

        # Extract the gene IDs from the specified column
        gene_ids <- data[[input$geneColumn]]
        gene_ids
    })

    # Reactive expression to perform the GO enrichment analysis
    enrichmentResults <- eventReactive(input$goButton, {
        req(geneList())  # Ensure the gene list is available

        gene_ids <- geneList()  # Get the gene IDs

        # Perform GO enrichment analysis using the enrichGO function
        enrich_result <- enrichGO(gene = gene_ids,
                                  OrgDb = org.Hs.eg.db,  # Specify the organism database (human in this case)
                                  keyType = "SYMBOL",  # The type of gene identifiers (e.g., SYMBOL)
                                  ont = "BP",  # Ontology: Biological Process (BP)
                                  pAdjustMethod = "BH",  # p-value adjustment method (Benjamini-Hochberg)
                                  qvalueCutoff = 0.05)  # Q-value cutoff for significance

        # Handle cases where no significant GO terms are found
        if (nrow(enrich_result) == 0) {
            stop("No significant GO terms found.")
        }

        return(enrich_result)
        print("GO enrichment completed.")
    })

    # Render the GO enrichment plot
    output$enrichmentPlot <- renderPlot({
        req(enrichmentResults())  # Ensure enrichment results are available

        enrich_result <- enrichmentResults()  # Get the enrichment results

        # Generate a bar plot for the top GO enrichment categories
        plot <- barplot(enrich_result, showCategory = 10) +  # Show top 10 categories
            ggtitle("GO Enrichment Analysis Results")  # Add a title to the plot

        return(plot)
    })

    # Download handler for the enrichment results
    output$downloadResults <- downloadHandler(
        filename = function() {
            paste("go_enrichment_results.csv")  # Define the filename for the downloaded results
        },
        content = function(file) {
            req(enrichmentResults())  # Ensure enrichment results are available

            enrich_result <- enrichmentResults()  # Get the enrichment results
            write.csv(as.data.frame(enrich_result), file, row.names = FALSE)  # Save the results as a CSV file
        }
    )

    # Download handler for the enrichment plot
    output$downloadPlot <- downloadHandler(
        filename = function() {
            paste("go_enrichment_plot.png")  # Define the filename for the downloaded plot
        },
        content = function(file) {
            req(enrichmentResults())  # Ensure enrichment results are available

            plot <- renderPlot({  # Re-render the plot for saving
                enrich_result <- enrichmentResults()  # Get the enrichment results

                # Generate the plot
                barplot(enrich_result, showCategory = 10) +
                    ggtitle("GO Enrichment Analysis Results")
            })

            # Save the plot as a PNG file
            ggsave(file, plot = plot, width = 10, height = 6)
        }
    )
}
