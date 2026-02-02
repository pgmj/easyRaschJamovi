#' Helper functions for easyRaschJamovi
#'
#' @import iarm
#' @import eRm
#' @import mirt
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import janitor
#' @import doParallel
#' @import purrr
#' @import catR
#' @import psychotools
#' @import readxl
#' @import readr
#' @import matrixStats
#' @import knitr
#' @import doFuture
#' @importFrom magrittr %>%
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom doFuture %dofuture%
#' @importFrom doFuture registerDoFuture
NULL

options(rgl.useNULL = TRUE)

RIcheckdata <- function(data, n = 3) {

  if (min(as.matrix(data), na.rm = T) > 0) {
    stop("The lowest response category needs to coded as 0. Please recode your data.")
  }

  # count number of responses per cell (item & response category)
  allcells_n <- data %>%
    pivot_longer(everything()) %>%
    dplyr::count(name,value, .drop = FALSE) %>%
    pull(n)

  # check whether any cells have fewer than 3 responses and return logical TRUE/FALSE
  return(any(allcells_n < n))
}

RIitemparams <- function(dfin, fontsize = 15, output = "dataframe",
                         detail = "thresholds", filename = "item_params.csv",
                         tbl_width = 90) {

  if(RIcheckdata(dfin) == TRUE) {
    warning("Warning! Your data has less than 3 responses in some response categories. Results may not be reliable.")
    # mirt is less unreliable than eRm in this situation
    mirt_out <- mirt(dfin, model=1, itemtype='Rasch', verbose = FALSE)
    item.locations <- coef(mirt_out, simplify = TRUE, IRTpars = TRUE)$items %>%
      as.data.frame() %>%
      dplyr::select(!a) %>%
      as.matrix()
    item.locations <- item.locations - mean(item.locations, na.rm = TRUE)
    maxcat <- dfin %>%
      pivot_longer(everything()) %>%
      dplyr::count(name,value) %>%
      pull(value) %>%
      max(na.rm = TRUE)
    item.locations <- item.locations %>%
      as.data.frame() %>%
      set_names(paste0("Threshold ", 1:maxcat))

  } else if(RIcheckdata(dfin) == FALSE) {
    erm_out <- PCM(dfin)
    item.locations <- as.data.frame(thresholds(erm_out)[[3]][[1]][, -1] - mean(thresholds(erm_out)[[3]][[1]][, -1], na.rm=T))
  }

  item_difficulty <- item.locations %>%
    mutate(Location = rowMeans(., na.rm = TRUE), .before = `Threshold 1`) %>%
    mutate(across(where(is.numeric), ~ round(.x, 3)))

  # detailed df
  item_params <- item_difficulty %>%
    mutate(all_item_avg = mean(Location)) %>%
    mutate(relative_avg_loc = Location - all_item_avg) %>%
    mutate(relative_lowest_tloc = `Threshold 1` - all_item_avg) %>%
    rownames_to_column("itemnr")

  # get the highest threshold value from each item - since the number of thresholds can vary, this needs special treatment
  highest_loc <- item_params %>%
    pivot_longer(cols = starts_with("Threshold"),
                 names_to = "threshold",
                 values_to = "t_location") %>%
    group_by(itemnr) %>%
    dplyr::filter(t_location == max(t_location, na.rm = TRUE)) %>%
    ungroup() %>%
    dplyr::select(itemnr, t_location) %>%
    dplyr::rename(highest_tloc = t_location)

  # join the highest threshold location to the item_params df
  item_params <- item_params %>%
    left_join(highest_loc, by = "itemnr") %>%
    mutate(relative_highest_tloc = highest_tloc - all_item_avg) %>%
    dplyr::relocate(all_item_avg, .after = relative_highest_tloc) %>%
    dplyr::select(-highest_tloc) %>%
    as.data.frame()

  if (output == "file" & detail == "thresholds") {
    item_difficulty %>%
      dplyr::select(!Location) %>%
      set_names(paste0("threshold_", 1:ncol(.))) %>%
      write_csv(., file = filename)
  }
  else if (output == "file" & detail == "all") {
    item_params %>%
      write_csv(., file = filename)
  }
  else if (output == "dataframe" & detail == "thresholds") {
    return(item_difficulty)
  }
  else if (output == "dataframe" & detail == "all") {
    return(item_params)
  }
  else if (output == "table" & detail == "thresholds") {
    item_difficulty %>%
      mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
      dplyr::relocate(Location, .after = last_col()) %>%
      dplyr::rename('Item location' = Location) %>%
      knitr::kable()
  }
  else if (output == "table" & detail == "all") {
    item_params %>%
      mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
      dplyr::rename('Item location' = Location,
                    'Relative item location' = relative_avg_loc,
                    'Relative lowest threshold' = relative_lowest_tloc,
                    'Relative highest threshold' = relative_highest_tloc) %>%
      dplyr::select(!all_item_avg)%>%
      knitr::kable()
  }

}

RIestThetas <- function(data) {
  # Estimate person parameters using WLE for dichotomous Rasch model data
  if (min(as.matrix(data), na.rm = T) > 0) {
    stop("The lowest response category needs to coded as 0. Please recode your data.")
  }
  
  if (max(as.matrix(data), na.rm = T) != 1) {
    stop("RIestThetas is only for dichotomous data. Use RIestThetasOLD for polytomous data.")
  }
  
  df.erm <- eRm::RM(data)
  person_params <- eRm::person.parameter(df.erm)
  
  # Return a list with WLE values
  list(WLE = person_params[["theta.table"]][["Person Parameter"]])
}

RIestThetasOLD <- function(data, itemParams, method = "WL",
                           theta_range = c(-10,10)) {

  if (min(as.matrix(data), na.rm = T) > 0) {
    stop("The lowest response category needs to coded as 0. Please recode your data.")

  } else if (max(as.matrix(data), na.rm = T) == 1) {
    model <- "RM"
  } else if (max(as.matrix(data), na.rm = T) > 1) {
    model <- "PCM"
  }

  # define function to call from purrr::map_dbl later.
  estTheta <- function(personResponse, itemParameters = itemParams, rmod = model,
                       est = method, rtheta = theta_range) {
    thetaEst(itemParameters, as.numeric(as.vector(personResponse)), model = rmod,
             method = est, range = rtheta)
  }
  # if no itemParams are given, calculate them based on input dataframe
  if (missing(itemParams) & model == "PCM") {
    erm_out <- PCM(data)
    itemParams <- thresholds(erm_out)[[3]][[1]][, -1] - mean(thresholds(erm_out)[[3]][[1]][, -1], na.rm=T)

    # Transpose dataframe to make persons to columns, then output a vector with thetas
    data %>%
      t() %>%
      as.data.frame() %>%
      map_dbl(., estTheta)

  } else if (missing(itemParams) & model == "RM") {
    df.erm <- RM(data)
    itemParams <- as.matrix(coef(df.erm, "beta")*-1)

    # Transpose dataframe to make persons to columns, then output a vector with thetas
    data %>%
      t() %>%
      as.data.frame() %>%
      map_dbl(., ~ estTheta(.x, rmod = NULL))

  } else if (!missing(itemParams) & model == "PCM") {

    # Transpose dataframe to make persons to columns, then output a vector with thetas
    data %>%
      t() %>%
      as.data.frame() %>%
      map_dbl(., estTheta)

  } else if (!missing(itemParams) & model == "RM") {
    df.erm <- RM(data)
    itemParams <- as.matrix(coef(df.erm, "beta")*-1)

    # Transpose dataframe to make persons to columns, then output a vector with thetas
    data %>%
      t() %>%
      as.data.frame() %>%
      map_dbl(., ~ estTheta(.x, rmod = NULL))

  }
}

### The two functions below were written by nicklas.korsell@ri.se

SimPolyItem <- function(x, thetas) {
  k <- length(x) + 1 # Antal kategorier för denna item
  n <- length(thetas) # Antal respondenter

  # Tar diffen mellan varje theta och varje delta
  Y <- outer(thetas, x, "-")

  # Kumulativa summan (theta -delta1) + (theta - delta2) + ...
  cumsums <-
    t(
      rbind(
        rep(0, times = n), # Lägg på en rad med 0:or
        apply(X = Y, MARGIN = 1, FUN = cumsum)
      )
    )

  # Exponentiera
  expcumsums <- exp(cumsums)

  # Beräkna nämnaren i normeringen
  norms <- apply(X = expcumsums, MARGIN = 1, FUN = sum)

  # Gör själva normeringen
  z <- expcumsums / norms

  # Nu gör vi 1 simulering för denna item för alla individer
  vapply(
    X = 1:n,
    FUN = function(x) {
      sample(x = 0:(k - 1), size = 1, replace = TRUE, prob = z[x, ])
    },
    FUN.VALUE = 1
  )
}

# Simulate Rasch partial credit model data
SimPartialScore <- function(deltaslist, thetavec) {
  # Gör om varje element i deltalist till en numerisk vector
  deltaslist <- lapply(X = deltaslist, FUN = unlist)

  # Anropa SimPolyItem för varje item
  sapply(
    X = deltaslist,
    FUN = SimPolyItem,
    thetas = thetavec
  )
}
