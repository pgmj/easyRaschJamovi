
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
            
            # Select only the specified variables
            df <- data[, vars]
            
            # Convert to numeric if needed
            for (col in names(df)) {
                df[[col]] <- as.numeric(df[[col]])
            }
            
            # Run the analysis
            tryCatch({
                # Get fit statistics
                gf <- RIgetfit(df, iterations = iterations, cpu = cores)
                fit_results <- RIitemfit(df, gf)
                
                # Get item parameters for thresholds
                item_params <- RIitemparams(df, output = "dataframe")
                
                # Get the fit table
                table <- self$results$fitTable
                
                # Populate the table
                for (i in seq_along(vars)) {
                    item_name <- vars[i]
                    
                    # Get infit MSQ value
                    infitmsq <- fit_results[fit_results$Item == item_name, "InfitMSQ"]
                    
                    # Get relative location
                    location <- fit_results[fit_results$Item == item_name, "Relative location"]
                    
                    # Get thresholds as text
                    threshold_cols <- grep("^Threshold", names(item_params), value = TRUE)
                    thresholds <- item_params[i, threshold_cols]
                    threshold_text <- paste(round(na.omit(as.numeric(thresholds)), 2), collapse = ", ")
                    
                    # Add row to table
                    table$setRow(rowNo = i, values = list(
                        item = item_name,
                        infitmsq = infitmsq,
                        location = location,
                        thresholds = threshold_text
                    ))
                }
                
            }, error = function(e) {
                # Handle errors gracefully
                stop(paste("Error in item fit analysis:", e$message))
            })
        }
    )
)
