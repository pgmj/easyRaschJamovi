# manually load functions to test them locally
source("R/helpers.R")
source("R/rescor.R")
source("R/getfit.R")
source("R/wrapper.R")

# test underlying functions for item fit
gf <- RIgetfit(pcmdat2)
RIitemfit(pcmdat2, gf)
# test wrapper that does both steps in one
rasch_itemfit(pcmdat2)


# test underlying functions for residual correlations
rescor <- RIgetResidCor(pcmdat2)
RIresidcorr(pcmdat2,rescor)
# test wrapper that does both steps in one
rasch_rescor(pcmdat2)

