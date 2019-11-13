## Auxiliary Functions
# (1) rs_checker   : either a list or 3d array
# (2) rs_transform : let's make it a nicer way




# (1) rs_checker ----------------------------------------------------------
#' @keywords internal
#' @noRd
rs_checker <- function(rsinput){
  if (is.array(rsinput)){       # 1. 3d array
    if (length(dim(rsinput))!=3){
      stop("* RiemRS : input array should be 3-dimensional.")
    }
    if (all(apply(rsinput, 3, rs_check_single)==TRUE)){
      return(TRUE)
    } else {
      return(FALSE)
    }
  } else if (is.list(rsinput)){ # 2. list of matrices
    if (length(rsinput)<2){
      stop("* RiemRS : input list should be of length at least 2.")
    }
    matss = unlist(base::lapply(rsinput, base::is.matrix))
    if (any(matss==FALSE)){
      stop("* RiemRS : all elements in a given list should be matrix type.")
    }
    nrows = unlist(base::lapply(rsinput, nrow))
    ncols = unlist(base::lapply(rsinput, ncol))
    if (length(unique(nrows))!=1){
      stop("* RiemRS : all matrices in a list should have same number of rows.")
    }
    if (length(unique(ncols))!=1){
      stop("* RiemRS : all matrices in a list should have same number of columns.")
    }
    if (all(unlist(base::lapply(rsinput, rs_check_single))==TRUE)){
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
}
#' @keywords internal
#' @noRd
rs_check_single <- function(rsmat){
  rsum = rowSums(rsmat)
  
  cond1 = (is.matrix(rsmat))
  cond2 = all(abs(rsum-rep(1,length(rsum))) < sqrt(.Machine$double.eps))
  cond3 = all(rsmat >= 0)
  if (cond1&&cond2&&cond3){
    return(TRUE)
  } else {
    return(FALSE)
  }
}


# (2) rs_transform --------------------------------------------------------
#' @keywords internal
#' @noRd
rs_transform <- function(rsinput, to.sphere=FALSE){
  # stack as 3d array
  if (is.list(rsinput)){
    m = nrow(rsinput[[1]])
    n = ncol(rsinput[[1]])
    p = length(rsinput)
    tmp = array(0,c(m,n,p))
    for (i in 1:p){
      tmp[,,i] = rsinput[[i]]
    }
  } else {
    tmp  = rsinput
    dims = dim(rsinput)
    m = dims[1]
    n = dims[2]
    p = dims[3]
  }
  
  # for each slice, let's do the transformation
  output = array(0,c(m,n,p))
  for (i in 1:p){
    output[,,i] = rs_trf_single(tmp[,,i])
  }
  
  # return
  if (to.sphere){
    return(sqrt(output))
  } else {
    return(output)
  }
}
#' @keywords internal
#' @noRd
rs_trf_single <- function(themat){
  sqrtval = sqrt(.Machine$double.eps)
  n = nrow(themat)
  p = ncol(themat)
  vec1p = 1:p
  outmat = array(0,c(n,p))
  for (i in 1:n){
    tgtvec = as.vector(themat[i,]) 
    if (any(tgtvec==0)){
      id0 = which(tgtvec==0)
      id1 = setdiff(vec1p, id0)
      tgtvec[id1] = tgtvec[id1] - rep(sqrtval/length(id1), length(id1))
      tgtvec[id0] = tgtvec[id0] + rep(sqrtval/length(id0), length(id0))
    }
    outmat[i,] = tgtvec/sum(tgtvec)
  }
  return(outmat)
}


