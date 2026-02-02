
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
                # Get fit statistics using simulations
                gf <- RIgetfit(df, iterations = iterations, cpu = cores)
                
                # Get item fit table (recreate logic from RIitemfit without kable output)
                if(max(as.matrix(df), na.rm = TRUE) == 1) {
                    erm_out <- eRm::RM(df)
                    item_avg_locations <- coef(erm_out, "beta") * -1
                    person_avg_locations <- eRm::person.parameter(erm_out)[["theta.table"]][["Person Parameter"]] %>%
                        mean(na.rm = TRUE)
                    relative_item_avg_locations <- item_avg_locations - person_avg_locations
                } else if(max(as.matrix(df), na.rm = TRUE) > 1) {
                    erm_out <- eRm::PCM(df)
                    item_avg_locations <- RIitemparams(df, output = "dataframe") %>%
                        pull(Location)
                    person_avg_locations <- RIestThetasOLD(df) %>%
                        mean(na.rm = TRUE)
                    relative_item_avg_locations <- item_avg_locations - person_avg_locations
                }
                
                # Get conditional MSQ
                cfit <- iarm::out_infit(erm_out)
                
                # Get item parameters for thresholds
                item_params <- RIitemparams(df, output = "dataframe")
                
                # Get the fit table
                table <- self$results$fitTable
                
                # Populate the table
                for (i in seq_along(vars)) {
                    item_name <- vars[i]
                    
                    # Get infit MSQ value
                    infitmsq <- round(cfit$Infit[i], 3)
                    
                    # Get relative location
                    location <- round(relative_item_avg_locations[i], 2)
                    
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
