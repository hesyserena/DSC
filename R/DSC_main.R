#' Main Ouputs of DSC
#'
#' Function for computing pre-treatment weights, optimal weights, post-treatment barycenters, post-treatment confidence intervals.
#' @param c_df a list of list: each sublist contains observations of every control units at a specific time period.
#' e.g. \code{c_df[[1]][[2]]} contains the observations of the second control unit at the first time period.
#' @param t_df a list of observations. e.g. \code{t_df[[1]]} contains the observations of the target at the first time period.
#' @param T0 an integer specifying the end of pre-treatment period.
#' @param M an integer specifying the number of draws from the uniform distribution for approximating the integral.
#' @param solver a solver for the optimization problem; see \code{\link[CVXR]{installed_solvers}} in CVXR for more options.
#' @param ww 0 or a vector. By default, i.e. 0, arithmetic mean is used for calculating optimal weights.
#' Otherwise, a vector of specific weights is used.
#' @param evgrid a vector of gridpoints used to evaluate quantile functions.
#' @param cl a float specifying the confidence level.
#' @param num.redraws an integer specifying the number of redraws used in the bootstrap approach.
#' @return \code{DSC_main} returns a list containing the following components:
#' \item{\code{weight_t}}{a list of weights in each pre-treatment period.}
#' \item{\code{weights}}{a vector of optimal weights.}
#' \item{\code{bc_post}}{a list of post-treatment barycenters.}
#' \item{\code{CI_post}}{a list of post-treatment confidence intervals.}
#' @examples
#' #simulated data from independent normal distributions
#' #ex_normal() calls the simulated data
#' #detail can be found by ??ex_normal
#' DSC_main(c_df=ex_normal()$control, t_df=ex_normal()$target, T0=5, M=100, solver="SCS", ww="arithmetic", evgrid=seq(from=0, to=1, length.out=101), cl=0.99, num.redraws=500)
#' @export
#'
#'



DSC_main=function(c_df, t_df, T0, M=100, solver="SCS", ww=0, evgrid=seq(from=0, to=1, length.out=101), cl=0.99, num.redraws=500){

  #calculate lambda_t for t<=T0
  lambda_t=list()

  for (t in 1:T0){
    lambda_t[[t]]=DSC_weights(c_df[[t]],as.vector(t_df[[t]]), M, choose_solver = solver)$weights

  }

  #calculate the average optimal lambda
  if (length(ww)==1){
    w_t=rep(1/T0, T0)
    lambda.opt=matrix(unlist(lambda_t),ncol=T0)%*%w_t
  }else{
    lambda.opt=matrix(unlist(lambda_t),ncol=T0)%*%ww
  }


  #calculate the barycenters for post-treatment period, t>T0
  bc_post=list()

  for (t in (T0+1):length(c_df)){
    bc_post[[t-T0]]=DSC_bc(c_df[[t]],lambda.opt,evgrid)$bc

  }

  #the confidence intervals for post-treatment period barycenters
  CI=list()
  for (t in (T0+1):length(c_df)){
    cat(paste('post-treatment CI, t =',t), '\n')
    CI[[t-T0]]=DSC_CI(df_control=c_df[[t]], df_target=as.vector(t_df[[t]]), M, solver="SCS", w_t=lambda.opt, cl, num.redraws, evgrid,
                      graph=FALSE, y_name='y', x_name='x')
  }

  #return pre-treatment weights, optimal weights, post-treatment barycenters, post-treatment confidence interval
  return(list(weight_t=lambda_t, weights=lambda.opt, bc_post=bc_post, CI_post=CI))
}










