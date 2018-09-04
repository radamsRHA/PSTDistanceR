#' Pipeline.SpeciesTreeDistances.ScaleAllBranches: function to measure phylogenetic distances across a range of branch length scaling between two species tree models
#'
#' This function returns a matrix containing the 5 distance metics (3 probabilistic + Robinson-Foulds + BHV)
#' @param vector.ScaleFactor Vector to containing the values to scale the branches of the second species tree model 
#' @param string.PathParentDir Path to a parent directory that will be used for subdirectories
#' @param string.PathHybridCoal Path to HYBRID-COAL executable
#' @param handle.SpeciesTree.Model1 Species Tree Model 1
#' @param handle.SpeciesTree.Model2 Species Tree Model 2
#' @keywords multispecies coalescent, model distance, incomplete lineage sorting, hybridization
#' @return XXX\cr \cr 
#' @export
#' @examples
#' 
#' library(phangorn)
#' library(PSTD)
#' library(distory)
#' 
#' string.SpeciesTree.1a <- "(((A:1,B:1):1,C:2):1,D:3);"
#' handle.SpeciesTree.1a <- read.tree(text = string.SpeciesTree.1a)
#' 
#' vector.Gamma <- seq(0.01, 10, 0.01)
#' handle.Figure1d <- Pipeline.SpeciesTreeDistances.ScaleAllBranches(vector.ScaleFactor = vector.Gamma, 
#'    string.PathParentDir = '~/Desktop/', 
#'    string.PathHybridCoal = '/Applications/hybrid-coal-v0.2.1-beta/hybrid-coal', 
#'    handle.SpeciesTree.Model1 = handle.SpeciesTree.1a, 
#'    handle.SpeciesTree.Model2 = handle.SpeciesTree.1a)
#'    
#' plot(handle.Figure1d[,1], ylim = c(0, 0.75), col = "darkblue", type = "l", xaxt='n', yaxt="n", lty = 1)
#' axis(1,at=seq(from = 0, to = 1000, by = 50),labels=seq(from = 0, to = 10, by = 0.50))
#' axis(2,at=seq(from = 0, to = 0.75, by = 0.10),labels=seq(from = 0, to = 0.75, by = 0.10))
#' lines(handle.Figure1d[,2], col = "blue", lty = 2)
#' lines(handle.Figure1d[,3], col = "lightblue", lty = 4)
#' lines(handle.Figure1d[,5]*0.025, col = "red", lty = 3)

##################################################
# Pipeline.SpeciesTreeDistances.ScaleAllBranches #
##################################################
Pipeline.SpeciesTreeDistances.ScaleAllBranches <- function(vector.ScaleFactor, string.PathParentDir, string.PathHybridCoal, handle.SpeciesTree.Model1, handle.SpeciesTree.Model2){
  
  ##################################################
  # Create directory for analysis with hybrid-coal #
  ##################################################
  string.ChildDir0 = paste(string.PathParentDir, '/Pipeline_PSTD_ScaleAllBranches_', Sys.Date(), sep = "")
  unlink(string.ChildDir0, recursive = T)
  dir.create(string.ChildDir0, showWarnings = T, recursive = T)
  matrix.DistanceResults <- matrix(nrow = length(vector.ScaleFactor), ncol = 5)
  rownames(matrix.DistanceResults) <- vector.ScaleFactor
  colnames(matrix.DistanceResults) <- c("Hellinger", "Kullback_Leibler", "Jensen_Shannon", "Robinson_Foulds", "BHV")
  
  ##############################
  # Loop through scale factors #
  ##############################
  for (i in 1:length(vector.ScaleFactor)){
    
    numeric.Scale.i <- vector.ScaleFactor[i]
    
    ##########################################
    # Create subdir for scale-based analyses #
    ##########################################
    string.ChildDir1 <- paste(string.ChildDir0, '/ScaleAllBranches_XXX', sep = "")
    string.ChildDir1 <- gsub(pattern = "XXX", replacement = numeric.Scale.i, x = string.ChildDir1)
    unlink(string.ChildDir1, recursive = T)
    dir.create(string.ChildDir1, showWarnings = T, recursive = T)
    
    ################################################
    # scale branch lengths of species tree model 2 #
    ################################################
    handle.SpeciesTree.Model2.i <- handle.SpeciesTree.Model2
    handle.SpeciesTree.Model2.i$edge.length <- handle.SpeciesTree.Model2.i$edge.length*numeric.Scale.i
    
    ############################
    # Conduct PSTD analyses... #
    ############################
    vector.Probabilistic.SpeciesTree.Distances <- Compute.Probabilistic.SpeciesTree.Distances(handle.SpeciesTree.Model1 = handle.SpeciesTree.Model1, 
                                                                                              handle.SpeciesTree.Model2 = handle.SpeciesTree.Model2.i, 
                                                                                              string.PathParentDir = string.ChildDir1, 
                                                                                              string.PathHybridCoal = string.PathHybridCoal)
    
    numeric.Hellinger <- vector.Probabilistic.SpeciesTree.Distances[1]
    numeric.Kullback_Leibler <- vector.Probabilistic.SpeciesTree.Distances[2]
    numeric.Jensen_Shannon <- vector.Probabilistic.SpeciesTree.Distances[3]
    
    ############################
    # Get BHV and RF distances #
    ############################
    numeric.Robinson_Foulds <- RF.dist(tree1 = handle.SpeciesTree.Model1, tree2 = handle.SpeciesTree.Model2.i, normalize = T, rooted = T)
    numeric.BHV <- dist.multiPhylo(x = c(handle.SpeciesTree.Model1, handle.SpeciesTree.Model2.i), method = "geodesic")
    
    ############################
    # Append to results matrix #
    ############################
    matrix.DistanceResults[i,] <- c(numeric.Hellinger, numeric.Kullback_Leibler, numeric.Jensen_Shannon, numeric.Robinson_Foulds, numeric.BHV)
    
  }
    
  return(matrix.DistanceResults = matrix.DistanceResults)
  
}