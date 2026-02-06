#' @export
rescorClass <- R6::R6Class(
  "rescorClass",
  inherit = rescorBase,
  private = list(
    .run = function() {
      # Get the data
      if (is.null(self$options$vars) || length(self$options$vars) == 0) {
        return()
      }

      # Check minimum number of variables
      if (length(self$options$vars) < 2) {
        stop("You need at least two variables to run an analysis.")
      }

      data <- self$data
      vars <- self$options$vars
      iterations <- self$options$iterations
      cores <- self$options$cores
      seed <- self$options$seed

      # Select only the specified variables (drop = FALSE keeps it as dataframe)
      df <- data[, vars, drop = FALSE]

      # Convert to numeric - handle factors properly
      for (col in names(df)) {
        if (is.factor(df[[col]])) {
          # For factors, convert via character to preserve actual values
          df[[col]] <- as.numeric(as.character(df[[col]]))
        } else {
          df[[col]] <- as.numeric(df[[col]])
        }
      }

      # Check for duplicate/identical variables using correlation
      check_identical_vars <- function(df_numeric) {
        n_vars <- ncol(df_numeric)
        identical_pairs <- list()

        for (i in 1:(n_vars - 1)) {
          for (j in (i + 1):n_vars) {
            # Use correlation - identical columns have r = 1
            cor_val <- cor(df_numeric[[i]], df_numeric[[j]], use = "complete.obs")
            if (!is.na(cor_val) && cor_val == 1) {
              identical_pairs <- append(identical_pairs,
                                        list(c(names(df_numeric)[i], names(df_numeric)[j])))
            }
          }
        }
        return(identical_pairs)
      }

      # Check AFTER conversion
      identical_pairs <- check_identical_vars(df)

      if (length(identical_pairs) > 0) {
        # Format the pairs for the message
        pair_strings <- sapply(identical_pairs, function(p) paste0("'", p[1], "' and '", p[2], "'"))
        pair_msg <- paste(pair_strings, collapse = "; ")

        if (ncol(df) == 2) {
          # Only two variables and they're identical - stop
          stop(paste("The two selected items are identical:", pair_msg,
                     "- please select different items."))
        } else {
          # More than two variables but some are identical - warning
          jmvcore::reject(
            "Warning: Some items appear to be identical ({pairs}). This may affect results.",
            pairs = pair_msg
          )
        }
      }

      # Check for all-NA columns (variables with no valid data)
      all_na_cols <- sapply(df, function(x) all(is.na(x)))
      if (any(all_na_cols)) {
        bad_vars <- names(df)[all_na_cols]
        stop(paste(
          "The following variables contain no valid numeric data:",
          paste(bad_vars, collapse = ", ")
        ))
      }

      # Check minimum response category coding (must start at 0)
      min_val <- min(as.matrix(df), na.rm = TRUE)
      if (min_val > 0) {
        stop("Variables/items need to be continuous (numeric integers) and the lowest response category must be coded as 0. Please recode your data.")
      }
      if (min_val < 0) {
        stop("Response categories cannot be negative. Please recode your data.")
      }

      # Check for sufficient complete cases
      n_complete <- sum(complete.cases(df))
      if (n_complete == 0) {
        stop("No complete cases found in the data. Each row must have responses for all selected items.")
      }
      if (n_complete < 30) {
        jmvcore::reject("Warning: Only {n} complete cases found. Results may be unreliable with small samples.",
          n = n_complete
        )
      }

      # Check for sufficient response variation per item
      for (col in names(df)) {
        unique_vals <- length(unique(na.omit(df[[col]])))
        if (unique_vals < 2) {
          stop(paste0(
            "Item '", col, "' has no variation in responses. ",
            "Each item needs at least two different response values."
          ))
        }
      }

      # Run the analysis
      tryCatch(
        {
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
        },
        error = function(e) {
          # Handle errors gracefully
          stop(paste("Error in residual correlation analysis:", e$message))
        }
      )
    }
  )
)
