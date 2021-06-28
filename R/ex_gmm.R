#' Simulated Data from Gaussian Mixture
#'
#' Example data for \code{DSC_weights}, \code{DSC_bc}, \code{DSC_plot}, \code{DSC_CI}.
#' There are 30 control units, each has a distribution of a mixture of 3 Gaussian distributions.
#' The target has a distribution of a mixture of 4 Gaussian distributions.
#' 
#' @export
#' @return
#' \item{\code{target}}{a vector.}
#' \item{\code{control}}{a matrix.}


ex_gmm=function(){
  ##### Things the researcher needs to choose
  # Required parameters for both examples
  numdraws <- 1000 # number of draws for each distribution
  sd.mult <- 5 # seed_multiplier to change the randomness in a controlled way when generating the
  # controls. The number used is immaterial, it is just used to get the same
  # random seed each time
  sd.target <- 1860 # seed for the target

  num.con <- 30 # number of control variables
  

  # Mixture of 3 Gaussians 
  con <- matrix(0, nrow=num.con, ncol=numdraws)
  for (ii in 1:num.con){
    set.seed(sd.mult*ii)
    # generating uniformly distributed weights in the unit simplex
    a1 <- matrix(runif(3), nrow=1)
    a1 <- sweep(a1, 1, rowSums(a1), FUN="/")
    components <- sample(1:3,prob=a1,size=numdraws,replace=TRUE)
    mus <- runif(3,-10,10)
    sigmas <- runif(3,0.5,6)
    sigmas <- (sigmas + t(sigmas))/2
    con[ii,] <- rnorm(numdraws)*sigmas[components]+mus[components]
  }
  
  
  # generating the target distribution as a mixture of 4 Normals
  set.seed(sd.target)
  a1 <- matrix(runif(4), nrow=1)
  a1 <- sweep(a1, 1, rowSums(a1), FUN="/")
  components <- sample(1:4,prob=a1,size=numdraws,replace=TRUE)
  mus <- runif(4,-10,10)
  sigmas <- runif(4,0.5,6)
  sigmas <- (sigmas + t(sigmas))/2
  treat <- rnorm(numdraws)*sigmas[components]+mus[components]
  
  
  con <- t(con)
  target <- treat

  return(list(target=target, control=con))
}





























