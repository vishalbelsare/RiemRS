#' k-Medoids Clustering for Row-Stochastic Matrices
#' 
#' 
#' @export
rs.kmedoids <- function(rsobj, k=2){
  ###############################################
  # compure pairwise distance matrix
  pdist = RiemRS::rs.pdist(rsobj, as.dist=TRUE)
  myk   = round(k)
  
  ###############################################
  # k-medoids algorithm based on RiemBaseExt
  tmpout = RiemBaseExt::rclust.kmedoids(pdist, k=myk)
    
  ###############################################
  # return output
  output = list()
  output$clustering = tmpout$clustering
  return(output)
}