#' Create a Plot for the Target Distribution and Barycenter
#'
#' Create a plot for the target distribution (in black) and barycenter (in red).
#' @param evgrid a vector of gridpoints used to evaluate quantile functions.
#' @param bc a vector of evaluated barycenter; can be calculated by \code{\link[DSC]{DSC_bc}}.
#' @param target a vector of observations.
#' @param y_name a string for the title of the y-axis.
#' @param x_name a string for the title of the x-axis.
#' @export
#' @examples
#' #simulated data from Gaussian mixture
#' #ex_gmm() calls the simulated data
#' #detail can be found by ??ex_gmm
#' #obtaining the optimal weights
#' DSC_res=DSC_weights(ex_gmm()$control,ex_gmm()$target, M=100)
#' weights=DSC_res[[1]]
#' DSC_res2=DSC_bc(ex_gmm()$control,weights,seq(from=0,to=1,length.out=1001))
#' DSC_plot(bc=DSC_res2[[2]],target=ex_gmm()$target)


DSC_plot=function(evgrid = seq(from=0, to=1, length.out=1001), bc, target,
                  y_name='y', x_name='x'){

  # function to compute quantile function
  myquant <- function(X,q){
    # sort if unsorted
    if (is.unsorted(X)) X <- sort(X)
    # compute empirical CDF
    X.cdf <- 1:length(X) / length(X)
    # obtain the corresponding empirical quantile
    return(X[which(X.cdf >= q)[1]])
  }

  # computing the target quantile function
  target.q <- mapply(myquant, evgrid, MoreArgs = list(X=target))

  # plot
  plot(evgrid, bc, xlab='',ylab='',
       type='l', col="red", lwd=4)
  lines(evgrid, target.q,lwd=4)
  legend("bottomright",legend = c("Target", "DSC"),
         col=c("black", "red"),
         lty= c(1,1), lwd = c(4,4), cex = 1.5)
  title(ylab=y_name, line=2, cex.lab=1.5)
  title(xlab=x_name, line=2, cex.lab=1.5)
}











