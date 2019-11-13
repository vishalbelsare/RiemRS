#' Pairwise Distance
#' 
#' 
#' @export
rs.pdist <- function(rsobj, as.dist=FALSE){
  ###############################################
  # Check and Rearrange as 3d array
  if (!rs_checker(rsobj)){
    stop("* rs.pdist : input 'rsobj' is invalid.")
  }
  spobj = rs_transform(rsobj, to.sphere=TRUE) # 3d array, each row is sphere
  
  ###############################################
  # Parameters
  N = dim(spobj)[3]  # number of slices
  outmat = array(0,c(N,N))
  for (i in 1:(N-1)){
    spmat1 = spobj[,,i]
    for (j in (i+1):N){
      spmat2 = spobj[,,j]
      outmat[i,j] <- outmat[j,i] <- rs.pdist.bi(spmat1, spmat2)
    }
  }
  
  ###############################################
  # Return
  if (as.dist){
    return(stats::as.dist(outmat))
  } else {
    return(outmat)
  }
}


#   -----------------------------------------------------------------------
#' @keywords internal
#' @noRd
rs.pdist.bi <- function(mat1, mat2){
  N = nrow(mat1)
  P = ncol(mat1)
  
  bivec = rep(0,N)
  for (n in 1:N){
    bivec[n] = ((base::acos(sum(as.vector(mat1[n,])*as.vector(mat2[n,]))))^2)
  }
  return(sqrt(sum(bivec)))
}
