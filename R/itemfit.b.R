#' @export
itemfitClass <- R6::R6Class(
    "itemfitClass",
    inherit = itemfitBase,
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
                # Get fit statistics using simulations
                results <- rasch_itemfit(df, iterations = iterations, cores = cores, seed = seed)
                
                # Get the table object
                table <- self$results$fitTable
                
                # Populate the table with results
                for (i in seq_len(nrow(results))) {
                    table$setRow(rowNo = i, values = list(
                        item = results$Item[i],
                        infitmsq = results$InfitMSQ[i],
                        thresholds = results$`Infit thresholds`[i],
                        infitdiff = results$`Infit diff`[i],
                        location = results$`Relative location`[i]
                    ))
                }

            }, error = function(e) {
                # Handle errors gracefully
                stop(paste("Error in item fit analysis:", e$message))
            })
        }
    )
)
