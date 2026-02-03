#' Get simulation based cutoff values for Yen's Q^3^ residual correlations
#'
#' Based on Christensen et al. (2017, DOI: 10.1177/0146621616677520),
#' uses package `mirt` (Chalmers, 2012) to assess local dependence.
#'
#' Uses a dataframe with response data to simulate residual correlation values
#' across n simulations based on estimated item & person locations.
#'
#' Results include mean, max and difference between the mean and max for each
#' iteration. Also, 95th, 99th, 99.5th and 99.9th percentile values are
#' calculated for use with `RIresidcorr()` as cutoff value, since the max value
#' may be spurious and dependent on number of iterations.
#'
#' Uses multi-core processing. To find how many cores you have on your computer,
#' use `parallel::detectCores()`. Remember to keep 1-2 cores free.
#'
#' @param data Dataframe with response data
#' @param iterations Number of simulation iterations (needed)
#' @param cpu Number of CPU cores to use (4 is default)
#' @param seed Optional random seed for reproducibility (default NULL)
#' @export
RIgetResidCor <- function(data, iterations = 150, cpu = 4, seed = NULL) {
  # set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # get sample size
  sample_n <- nrow(data)

  if (min(as.matrix(data), na.rm = T) > 0) {
    stop("The lowest response category needs to coded as 0. Please recode your data.")
  } else if (max(as.matrix(data), na.rm = T) > 1 && min(as.matrix(data), na.rm = T) == 0) {
    # get item threshold locations for response data
    item_locations <- RIitemparams(data, output = "dataframe") %>%
      dplyr::select(!Location) %>%
      janitor::clean_names() %>%
      as.matrix()

    n_items <- nrow(item_locations)

    # item threshold locations in list format for simulation function
    itemlist <- list()
    for (i in 1:n_items) {
      itemlist[[i]] <- list(na.omit(item_locations[i, ]))
    }

    # get number of response categories for each item for later use in checking complete responses
    itemlength <- list()
    for (i in 1:n_items) {
      itemlength[i] <- length(na.omit(item_locations[i, ]))
      names(itemlength)[i] <- names(data)[i]
    }

    # estimate theta values in response data
    thetas <- RIestThetasOLD(data)

    # create object to store results from multicore loop
    residcor <- list()
    residcor <- foreach(i = 1:iterations,
      .options.future = list(seed = TRUE)
    ) %dofuture% {
      # resampled vector of theta values (based on sample properties)
      inputThetas <- sample(thetas, size = sample_n, replace = TRUE)

      # simulate response data based on thetas and items above
      testData <- SimPartialScore(
        deltaslist = itemlist,
        thetavec = inputThetas
      ) %>%
        as.data.frame()

      names(testData) <- names(data)

      # check that simulated dataset has responses in all categories
      data_check <- testData %>%
        # make factor to not drop any consequtive response categories with 0 responses
        mutate(across(everything(), ~ factor(.x, levels = c(0:itemlength[[as.character(expression(.x))]])))) %>%
        pivot_longer(everything()) %>% # screws up factor levels, which makes the next step necessary
        dplyr::count(name, value, .drop = FALSE) %>%
        pivot_wider(
          names_from = "name",
          values_from = "n"
        ) %>%
        dplyr::select(!value) %>%
        # mark missing cells with NA for later logical examination
        mutate(across(everything(), ~ car::recode(.x, "0=NA", as.factor = FALSE))) %>%
        as.data.frame() %>%
        dplyr::select(all_of(names(data))) # get item sorting correct

      # match response data generated with itemlength
      item_ccount <- list()
      for (i in 1:n_items) {
        item_ccount[i] <- list(data_check[c(1:itemlength[[i]]), i])
      }

      # check if any item has 0 responses in a response category that should have data
      if (any(is.na(unlist(item_ccount)))) {
        return("Missing cells in generated data.")
      }

      # create Yen's Q3 residual correlation matrix
      mirt.rasch <- mirt(testData, model = 1, itemtype = "Rasch", verbose = FALSE, accelerate = "squarem")
      resid <- residuals(mirt.rasch, type = "Q3", digits = 2, verbose = FALSE)
      diag(resid) <- NA

      data.frame(
        mean = mean(resid, na.rm = TRUE),
        max = max(resid, na.rm = TRUE)
      )
    }
  } else if (max(as.matrix(data), na.rm = T) == 1 && min(as.matrix(data), na.rm = T) == 0) {
    # estimate item threshold locations from data
    erm_out <- eRm::RM(data)
    item_locations <- erm_out$betapar * -1
    names(item_locations) <- names(data)

    # estimate theta values from data using MLE (temporary fix for RM with missing data)
    thetas <- eRm::person.parameter(erm_out)[["theta.table"]][["Person Parameter"]]

    # create object to store results from multicore loop
    residcor <- list()
    residcor <- foreach(i = 1:iterations) %dopar% {
      # resample vector of theta values (based on sample properties)
      inputThetas <- sample(thetas, size = sample_n, replace = TRUE)

      # simulate response data based on thetas and items above
      testData <-
        psychotools::rrm(inputThetas, item_locations, return_setting = FALSE) %>%
        as.data.frame()

      # TEMPORARY FIX START
      # check that all items have at least 8 positive responses, otherwise eRm::RM() fails
      n_resp <-
        testData %>%
        as.matrix() %>%
        colSums2() %>%
        t() %>%
        as.vector()

      if (min(n_resp, na.rm = TRUE) < 8) {
        return("Missing cells in generated data.")
      }
      # END TEMP FIX

      # create Yen's Q3 residual correlation matrix
      mirt.rasch <- mirt(testData, model = 1, itemtype = "Rasch", verbose = FALSE, accelerate = "squarem")
      resid <- residuals(mirt.rasch, type = "Q3", digits = 2, verbose = FALSE)
      diag(resid) <- NA

      data.frame(
        mean = mean(resid, na.rm = TRUE),
        max = max(resid, na.rm = TRUE)
      )
    }
  }

  # identify datasets with inappropriate missingness
  nodata <- lapply(residcor, is.character) %>% unlist()
  iterations_nodata <- which(nodata)

  actual_iterations <- iterations - length(iterations_nodata)

  # get all results to a dataframe
  if (actual_iterations == iterations) {
    results <-
      bind_rows(residcor) %>%
      mutate(diff = max - mean)
  } else {
    results <-
      bind_rows(residcor[-iterations_nodata]) %>%
      mutate(diff = max - mean)
  }

  out <- list()
  out$results <- results
  out$actual_iterations <- actual_iterations

  out$sample_n <- sample_n
  out$sample_summary <- summary(thetas)

  out$max_diff <- max(results$diff)
  out$sd_diff <- sd(results$diff)
  out$p95 <- quantile(results$diff, .95)
  out$p99 <- quantile(results$diff, .99)
  out$p995 <- quantile(results$diff, .995)
  out$p999 <- quantile(results$diff, .999)

  return(out)
}

#' Correlation matrix of Rasch residuals using Yen's Q^3^
#'
#' Uses package `mirt` (Chalmers, 2012).
#'
#' Mandatory option to set relative cutoff-value over the average of all
#' item residual correlations. It is strongly recommended to use the function
#' `RIgetResidCor()` to retrieve an appropriate cutoff value for your data.
#'
#' @param data Dataframe with item data only
#' @param cutoff Relative value above the average of all item residual correlations
#' @param output Default HTML table, optional "quarto" or "dataframe"
#' @export
RIresidcorr <- function(data, cutoff) {

  cutoff_used <- cutoff$p995
  cutoff_used99 <- cutoff$p99
  iter <- cutoff$actual_iterations

  mirt.rasch <- mirt(data, model = 1, itemtype = "Rasch", verbose = FALSE, accelerate = "squarem") # unidimensional Rasch model
  resid <- residuals(mirt.rasch, type = "Q3", digits = 2, verbose = FALSE) # get Q3 residuals

  diag(resid) <- NA # make the diagonal of correlation matrix NA instead of 1
  resid <- as.data.frame(resid)

  mean.resid <- resid %>%
    as.matrix() %>%
    mean(na.rm = TRUE)

  resid[upper.tri(resid)] <- NA

  resid <- resid %>%
    mutate_if(is.character, as.numeric) %>%
    mutate(across(everything(), ~ round(.x, 3))) %>% 
    rowwise() %>%
    mutate(
      `LD indicated in row` = case_when(
        any(c_across(everything()) > cutoff_used) ~ "x",
        TRUE ~ NA_character_
      )
    ) %>%
    ungroup()
  
  #ld_indicated <- resid$ld_indicated
  #resid$ld_indicated <- NULL
  
  #resid[upper.tri(resid)] <- "" # remove values in upper right triangle to clean up table
  #diag(resid) <- "" # same for diagonal

  rescor_results <- list(resmat = resid,
                         cutoff = cutoff_used,
                         cutoff99 = cutoff_used99)
  return(rescor_results)
}

