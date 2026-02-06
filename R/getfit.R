#' Get simulation based cutoff values for item fit values
#'
#' This function uses your response data to simulate datasets that fit the
#' Rasch model to find a credible range of item fit values. The function
#' outputs an object that is strongly recommended to save to an object, since it
#' takes some time to run this function when using many iterations/simulations.
#'
#' The output is a list object, which can in turn be used with two different
#' functions. Most importantly, you can use it with `RIitemfit()` to get
#' conditional highlighting of cutoff values based on your sample size and
#' item parameters. Each item gets its own cutoff thresholds.
#'
#' The function `RIgetfitPlot()` uses the package `ggdist` to plot the
#' distribution of fit values from the simulation results.
#'
#' Uses multi-core processing. To find how many cores you have on your computer,
#' use `parallel::detectCores()`. Remember to keep 1-2 cores free.
#'
#' Since version 0.2.4.2, the default is to only use complete cases in the
#' simulations, since this is what the conditional item fit function uses and
#' numbers should be more comparable using this method.
#'
#' @param data Dataframe with response data
#' @param iterations Number of simulation iterations (use 200-400)
#' @param cpu Number of CPU cores to use
#' @param na.omit Defaults to TRUE to produce conditional fit comparable values
#' @param seed Optional random seed for reproducibility (default NULL)
#' @export
RIgetfit <- function(data, iterations = 150, cpu = 4, na.omit = TRUE, seed = NULL) {
  # since we want comparable values to conditional item fit, which only uses
  # complete cases, we remove any missing responses by default
  if (na.omit == TRUE) {
    data <- na.omit(data)
  }
  sample_n <- nrow(data)

  # set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  } else {
    # ensure .Random.seed exists in .GlobalEnv
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      runif(1)
    }
  }
  
  # generate vector of random seeds for reproducible simulations
  seeds <- sample.int(.Machine$integer.max, iterations)

  registerDoParallel(cores = cpu)

    if (is.null(ncol(data)) || ncol(data) < 2) {
    stop("You need at least two variables to run an analysis.")
  }

  if (min(as.matrix(data), na.rm = T) > 0) {
    stop("Variables need to be numeric and the lowest response category coded as 0. Please recode your data.")
  } else if (max(as.matrix(data), na.rm = T) == 1 && min(as.matrix(data), na.rm = T) == 0) {
    # estimate item threshold locations from data
    erm_out <- eRm::RM(data)
    item_locations <- erm_out$betapar * -1
    names(item_locations) <- names(data)

    # estimate theta values from data using WLE
    mirt_out <- mirt(data, itemtype = "Rasch", verbose = FALSE)
    thetas <- mirt::fscores(mirt_out, method = "WLE", verbose = FALSE)

    fitstats <- list()
    #registerDoRNG(seeds[17])
    fitstats <- foreach(i = 1:iterations) %dopar% {
      # reproducible seed
      set.seed(seeds[i])
      # resampled vector of theta values (based on sample properties)
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

      # get conditional MSQ
      rm_out <- eRm::RM(testData, se = FALSE)
      cfit <- iarm::out_infit(rm_out)

      # create dataframe
      item.fit.table <- data.frame(
        InfitMSQ = cfit$Infit,
        OutfitMSQ = cfit$Outfit
      ) %>%
        round(3) %>%
        rownames_to_column("Item")
    }
  } else if (max(as.matrix(data), na.rm = T) > 1 && min(as.matrix(data), na.rm = T) == 0) {
    # estimate item threshold locations from data
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

    # estimate theta values from data using WLE
    thetas <- RIestThetasOLD(data)

    fitstats <- list()
    fitstats <- foreach(i = 1:iterations) %dopar% {
      # reproducible seed
      set.seed(seeds[i])
      # resampled vector of theta values (based on sample properties)
      inputThetas <- sample(thetas, size = sample_n, replace = TRUE)

      # simulate response data based on thetas and items above
      testData <- SimPartialScore(
        deltaslist = itemlist,
        thetavec = inputThetas
      ) %>%
        as.data.frame()

      names(testData) <- names(data)

      # check that data has responses in all categories
      data_check <- testData %>%
        # make factor to not drop any consecutive response categories with 0 responses
        mutate(across(everything(), ~ factor(.x, levels = c(0:itemlength[[as.character(expression(.x))]])))) %>%
        pivot_longer(everything()) %>% # screws up factor levels, which makes the next step necessary
        dplyr::count(name, value, .drop = FALSE) %>%
        pivot_wider(
          names_from = "name",
          values_from = "n"
        ) %>%
        dplyr::select(!value) %>%
        # mark missing cells with NA for later logical examination with if(is.na)
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

      # get conditional MSQ
      pcm_out <- psychotools::pcmodel(testData, hessian = FALSE)
      cfit <- iarm::out_infit(pcm_out)

      # create dataframe
      item.fit.table <- data.frame(
        InfitMSQ = cfit$Infit,
        OutfitMSQ = cfit$Outfit
      ) %>%
        round(3) %>%
        rownames_to_column("Item")
    }
  }

  fitstats$sample_n <- sample_n
  fitstats$sample_summary <- summary(thetas)

  return(fitstats)
}

#' Calculate conditional infit MSQ statistics
#'
#' Automatically uses RM (dichotomous data) or PCM (polytomous data) depending
#' on data structure.
#'
#' Uses `iarm::out_infit()` to calculate conditional mean square fit statistics
#' for all items. See Müller (2020, DOI: 10.1186/s40488-020-00108-7) for details.
#' Note: only uses complete cases! This is explicitly mentioned in the automatic
#' table caption text.
#'
#' Simulated datasets that have zero responses in any response category that
#' should have data will automatically be removed/skipped from analysis,
#' which means that final set of iterations may be lower than specified by user.
#'
#'
#' @param data Dataframe with response data
#' @param simcut Object output from `RIgetfit()`
#' @param output Optional "dataframe" or "quarto"
#' @param cutoff Default `c(.001,.999)`
#' @export
RIitemfit <- function(data, simcut, cutoff = c(.001,.999)) {

  if(min(as.matrix(data), na.rm = T) > 0) {
    stop("The lowest response category needs to coded as 0. Please recode your data.")
  } else if(na.omit(data) %>% nrow() == 0) {
    stop("No complete cases in data.")
  } else if(max(as.matrix(data), na.rm = T) == 1) {
    erm_out <- eRm::RM(data)
    item_avg_locations <- coef(erm_out, "beta")*-1 # item coefficients
    person_avg_locations <- eRm::person.parameter(erm_out)[["theta.table"]][["Person Parameter"]] %>%
      mean(na.rm = TRUE)
    relative_item_avg_locations <- item_avg_locations - person_avg_locations
  } else if(max(as.matrix(data), na.rm = T) > 1) {
    erm_out <- eRm::PCM(data)
    item_avg_locations <- RIitemparams(data, output = "dataframe") %>%
      pull(Location)
    person_avg_locations <- RIestThetasOLD(data) %>%
      mean(na.rm = TRUE)
    relative_item_avg_locations <- item_avg_locations - person_avg_locations
  }

  # get conditional MSQ
  cfit <- iarm::out_infit(erm_out)
  # get count of complete cases
  n_complete <- nrow(na.omit(data))

  # create dataframe
  item.fit.table <- data.frame(InfitMSQ = cfit$Infit) %>%
    round(3) %>%
    rownames_to_column("Item") %>%
    add_column(`Relative location` = round(relative_item_avg_locations,2))

  if (!missing(simcut)) {

    # get number of iterations used to get simulation based cutoff values
    iterations <- length(simcut) - 2

    nodata <- lapply(simcut, is.character) %>% unlist()
    iterations_nodata <- which(nodata)

    actual_iterations <- iterations - length(iterations_nodata)

    # summarise simulations and set cutoff values
    if (actual_iterations == iterations) {
      lo_hi <-
        bind_rows(simcut[1:iterations]) %>%
        group_by(Item) %>%
        summarise(min_infit_msq = quantile(InfitMSQ, cutoff[1]),
                  max_infit_msq = quantile(InfitMSQ, cutoff[2])
        )
    } else {
      lo_hi <-
        bind_rows(simcut[1:iterations][-iterations_nodata]) %>%
        group_by(Item) %>%
        summarise(min_infit_msq = quantile(InfitMSQ, cutoff[1]),
                  max_infit_msq = quantile(InfitMSQ, cutoff[2])
        )
    }

    lo_hi$Item <- names(data)

    # get upper/lower values into a dataframe
    if (actual_iterations == iterations) {
      fit_table <-
        bind_rows(simcut[1:iterations]) %>%
        group_by(Item) %>%
        summarise(inf_thresh = paste0("[",round(quantile(InfitMSQ, cutoff[1]),3),", ",round(quantile(InfitMSQ, cutoff[2]),3),"]")
        )
    } else {
      fit_table <-
        bind_rows(simcut[1:iterations][-iterations_nodata]) %>%
        group_by(Item) %>%
        summarise(inf_thresh = paste0("[",round(quantile(InfitMSQ, cutoff[1]),3),", ",round(quantile(InfitMSQ, cutoff[2]),3),"]")
        )
    }
    # add thresholds to dataframe and calculate differences between thresholds and observed values
    item.fit.table <-
      item.fit.table %>%
      add_column(`Infit thresholds` = fit_table$inf_thresh, .after = "InfitMSQ") %>%
      left_join(lo_hi, by = "Item") %>%
      mutate(infit_lo = abs(InfitMSQ - min_infit_msq),
             infit_hi = abs(InfitMSQ - max_infit_msq),
             `Infit diff` = round(pmin(infit_lo,infit_hi),3)
      ) %>%
      mutate(`Infit diff` = ifelse(yes = "no misfit", no = `Infit diff`, InfitMSQ > min_infit_msq & InfitMSQ < max_infit_msq)) %>%
      dplyr::select(!contains(c("lo","hi","min","max"))) %>%
      add_column(`Relative location` = round(relative_item_avg_locations,2))

    return(item.fit.table)
  }
}
