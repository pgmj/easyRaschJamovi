# These are the two main functions

rasch_itemfit <- function(data, iterations = 150, cores = 4) {
  gf <- RIgetfit(data, iterations = iterations, cpu = cores)
  RIitemfit(data,gf)
}

rasch_rescor <- function(data, iterations = 250, cores = 4) {
  rescor <- RIgetResidCor(data, iterations = iterations, cpu = cores)
  RIresidcorr(data,rescor)
}
