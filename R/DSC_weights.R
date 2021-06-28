#' Weights in the DSC Method
#'
#' Function for obtaining the weights in the DSC method at a specific time period.
#' @param controls a list or matrix:
#' in list form, each list element contains a vector of observations for the given control unit;
#' in matrix form, each column corresponds to one unit and each row is one observation.
#' @param target a vector of observations.
#' @param M an integer specifying the number of draws from the uniform distribution for approximating the integral.
#' @param choose_solver a solver for the optimization problem; see \code{\link[CVXR]{installed_solvers}} in CVXR for more options.
#' @return \code{DSC_weights} returns a list containing the following components:
#' \item{\code{weights}}{a matrix of the optimal weights.}
#' \item{\code{distance}}{a float of the squared Wasserstein distance between the target and the corresponding barycenter.}
#' @examples
#' #simulated data from Gaussian mixture
#' #ex_gmm() calls the simulated data
#' #detail can be found by ??ex_gmm
#' DSC_weights(ex_gmm()$control,ex_gmm()$target, M=100)
#' @export
#'
#'
#' @import CVXR

DSC_weights <- function(controls,target, M = 500, choose_solver = "SCS"){
  # the control distributions can be given in list form, where each list element contains a
  # vector of observations for the given control unit, in matrix form;
  # in matrix- each column corresponds to one unit and each row is one observation.
  # The list-form is useful, because the number of draws for each control group can be different.
  # The target must be given as a vector
  if (!is.vector(target)){
    stop("Target needs to be given as a vector.")
  }
  if (!is.list(controls) && !is.matrix(controls)){
    stop ("Controls need to be given in either list-form or matrix form.")
  }
  # if the controls are given in matrix form we turn them into a list
  if (is.matrix(controls)) {
    controls.h <- list()
    for (ii in 1:ncol(controls)){
      controls.h[[ii]] <- as.vector(controls[,ii])
    }
    controls <- controls.h
  }
  # M is the number of draws from the uniform distribution for approximating the integral

  ## Creating a function for the empirical quantile function
  myquant <- function(X,q){
    # sort if unsorted
    if (is.unsorted(X)) X <- sort(X)
    # compute empirical CDF
    X.cdf <- 1:length(X) / length(X)
    # obtain the corresponding empirical quantile
    return(X[which(X.cdf >= q)[1]])
  }

  ## Sampling from this quantile function M times
  Mvec <- runif(M, min = 0, max = 1)
  controls.s <- matrix(0,nrow = M, ncol = length(controls))
  for (jj in 1:length(controls)){
    controls.s[,jj] <- mapply(myquant, Mvec, MoreArgs = list(X=controls[[jj]]))
  }

  target.s <- matrix(0, nrow = M, ncol=1)
  target.s[,1] <- mapply(myquant, Mvec, MoreArgs = list(X=target))

  ## Solving the optimization using CVXR
  # the variable we are after
  theweights <- Variable(length(controls))
  # the objective function
  objective <- sum_squares(controls.s %*% theweights - target.s)
  # the constraints for the unit simplex
  constraints <- list(theweights>=0, sum_entries(theweights) == 1)
  # the optimization problem
  problem <- Problem(Minimize(objective),constraints)
  # solving the optimization problem
  results <- solve(problem, solver = choose_solver)

  # returning the optimal weights and the value function which provides the
  # squared Wasserstein distance between the target and the corresponding barycenter
  theweights.opt <- results$getValue(theweights)
  thedistance.opt <- results$value/M

  if (results$status != "optimal") {
    stop(paste('Convex optimization with CVXR for obtaining the optimal weights did ',
               'not find a solution. Take a smaller value for  M.', sep=""))
  }

  return(list(weights=theweights.opt, distance=thedistance.opt))
}
