#' @export
rescorClass <- R6::R6Class(
    "rescorClass",
    inherit = rescorBase,
    private = list(
        .run = function() {
            
            # Get the data
            if (is.null(self$options$vars) || length(self$options$vars) == 0)
                return()
            
            data <- self$data
            vars <- self$options$vars
            iterations <- self$options$iterations
            cores <- self$options$cores
            seed <- self$options$seed
            
            # Select only the specified variables
            df <- data[, vars]
            
            # Convert to numeric if needed
            for (col in names(df)) {
                df[[col]] <- as.numeric(df[[col]])
            }
            
            # Run the analysis
            tryCatch({
                # Get residual correlation cutoff
                rescor <- RIgetResidCor(df, iterations = iterations, cpu = cores, seed = seed)
                
                # Get the correlation matrix and cutoff
                resid_results <- RIresidcorr(df, rescor)
                resid_matrix <- resid_results$resmat
                cutoff_value <- resid_results$cutoff
                cutoff_value99 <- resid_results$cutoff99
                
                # Set cutoff information as HTML
                cutoffHtml <- self$results$cutoffText
                cutoff_html <- paste0(
                    "<p>Simulation-based relative cutoff <strong>(99.5th percentile):</strong> ", 
                    round(cutoff_value, 3), 
                    "</p>",
                    "<p>Simulation-based relative cutoff <strong>(99th percentile):</strong> ", 
                    round(cutoff_value99, 3), 
                    "</p>",
                    "<p>Residual correlations above the cutoff value may indicate local dependence.</p>"
                )
                cutoffHtml$setContent(cutoff_html)
                
                # Get the table object
                table <- self$results$residcorr
                
                # Get item names (column names from the matrix)
                item_names <- colnames(resid_matrix)
                
                # Add columns dynamically based on item names
                # First column is the row label
                table$addColumn(name = "rowname", title = "", type = "text")
                
                # Add a column for each item
                for (item in item_names) {
                    table$addColumn(name = item, title = item, type = "text")
                }
                
                # Add rows with data
                for (i in seq_len(nrow(resid_matrix))) {
                    row_values <- list(rowname = item_names[i])
                    for (j in seq_along(item_names)) {
                        row_values[[item_names[j]]] <- as.character(resid_matrix[i, j])
                    }
                    table$addRow(rowKey = i, values = row_values)
                }

            }, error = function(e) {
                # Handle errors gracefully
                stop(paste("Error in residual correlation analysis:", e$message))
            })
        }
    )
)
