#' Simulated data from independent normal distributions
#'
#' Example data for \code{DSC_main}, \code{DSC_per}.
#' There are 10 periods. Post-treatment periods are 1:5.
#' There are 30 control units. 
#' In each period, each control unit has a distribution of normal distribution with randomly generated mean in unif(0,1), 
#' ramdomly gnerated variance in unif(0.8,1.2).
#' The target has a distribution of N(0,1) in t=1:5, and N(2,1) in t=6:10.
#' 
#' @export
#' @return
#' \item{\code{target}}{a list.}
#' \item{\code{control}}{a list.}



ex_normal=function(){
  
  target=list()
  
  for (i in 1:5){
    set.seed(666+i)
    target[[i]]=rnorm(1000,0,1)
  }
  
  
  for (i in 6:10){
    set.seed(666+i)
    target[[i]]=rnorm(1000,2,1)
  }
  
  
  control=list()
  for (i in 1:10){
    control[[i]]=list()
  }
  
  
  for (j in 1:30){
    set.seed(888+j)
    m=runif(1,0,1)
    sigma=runif(1,0.8,1.2)
    for (i in 1:10){
      control[[i]][[j]]=rnorm(1000,m,sigma)
    }
  }
  
  return(list(target=target, control=control))
}

