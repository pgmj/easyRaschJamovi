
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
            
            # Select only the specified variables
            df <- data[, vars]
            
            # Convert to numeric if needed
            for (col in names(df)) {
                df[[col]] <- as.numeric(df[[col]])
            }
            
            # Run the analysis
            tryCatch({
                # Get residual correlation cutoff
                rescor <- RIgetResidCor(df, iterations = iterations, cpu = cores)
                
                # Get the correlation matrix as HTML
                resid_table <- RIresidcorr(df, rescor)
                
                # Convert kable object to HTML string
                html_output <- as.character(resid_table)
                
                # Set the HTML content
                html <- self$results$residcorr
                html$setContent(html_output)
                
            }, error = function(e) {
                # Handle errors gracefully
                stop(paste("Error in residual correlation analysis:", e$message))
            })
        }
    )
)
