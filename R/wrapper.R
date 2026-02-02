# These are the two main functions

rasch_itemfit <- function(data, iterations = 150, cores = 4, seed = NULL) {
  gf <- RIgetfit(data, iterations = iterations, cpu = cores, seed = seed)
  RIitemfit(data,gf)
}

rasch_rescor <- function(data, iterations = 250, cores = 4, seed = NULL) {
  rescor <- RIgetResidCor(data, iterations = iterations, cpu = cores, seed = seed)
  RIresidcorr(data,rescor)
}
