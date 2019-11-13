#' Mean
#' 
#' @export
rs.mean <- function(rsobj, mode=c("intrinsic","extrinsic")){
  ###############################################
  # Check and Rearrange as 3d array
  if (!rs_checker(rsobj)){
    stop("* rs.mean : input 'rsobj' is invalid.")
  }
  spobj = rs_transform(rsobj, to.sphere=TRUE) # 3d array, each row is sphere
  dimsp = dim(spobj)
  
  ###############################################
  # Parameters
  N = dimsp[1] # nrow
  P = dimsp[2] # ncol
  M = dimsp[3] # number of copies
  
  mode   = match.arg(mode)
  output = array(0,c(N,P))
  for (n in 1:N){
    ttgt = t(spobj[n,,]) # (M x P) matrix; rows are sphere vectors
    if (all(mode=="intrinsic")){
      mvec = as.vector(RiemSphere::mle.spnorm(ttgt)$mu)
    } else {
      etgt = base::colMeans(ttgt) # extrinsic mean
      mvec = as.vector(etgt/sqrt(sum(etgt^2)))
    }
    output[n,] = as.vector(mvec^2)
  }
  
  
  ###############################################
  # Return
  return(output)
}
