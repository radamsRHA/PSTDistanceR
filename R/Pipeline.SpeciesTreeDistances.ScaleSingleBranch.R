#' Pipeline.SpeciesTreeDistances.ScaleSingleBranch: function to measure phylogenetic distances across a range of branch length scaling (for a single branch on the second model) between two species tree models
#'
#' This function returns a matrix containing the 5 distance metics (3 probabilistic + Robinson-Foulds + BHV) across a range of scaler values
#' @param vector.ScaleFactor Vector to containing the values to scale the branches of the second species tree model 
#' @param string.PathParentDir Path to a parent directory that will be used for subdirectories
#' @param string.PathHybridCoal Path to HYBRID-COAL executable
#' @param handle.SpeciesTree.Model1 Species Tree Model 1
#' @param handle.SpeciesTree.Model2 Species Tree Model 2
#' @param numeric.BranchToScale Number indicating the branch to scale on species tree model 2
#' @keywords multispecies coalescent, model distance, incomplete lineage sorting, hybridization
#' @return XXX\cr \cr 
#' @export
#' @examples
#' 


###################################################
# Pipeline.SpeciesTreeDistances.ScaleSingleBranch #
###################################################
Pipeline.SpeciesTreeDistances.ScaleSingleBranch <- function(vector.ScaleFactor, string.PathParentDir, string.PathHybridCoal, handle.SpeciesTree.Model1, handle.SpeciesTree.Model2, numeric.BranchToScale){
 
  ##################################################
  # Create directory for analysis with hybrid-coal #
  ##################################################
  string.ChildDir0 = paste(string.PathParentDir, '/Pipeline_PSTD_ScaleSingleBranch_', numeric.BranchToScale, '_', Sys.Date(), sep = "")
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
    string.ChildDir1 <- paste(string.ChildDir0, '/ScaleSingleBranches_YYY_XXX', sep = "")
    string.ChildDir1 <- gsub(pattern = "XXX", replacement = numeric.Scale.i, x = string.ChildDir1)
    string.ChildDir1 <- gsub(pattern = "YYY", replacement = numeric.BranchToScale, x = string.ChildDir1)
    unlink(string.ChildDir1, recursive = T)
    dir.create(string.ChildDir1, showWarnings = T, recursive = T)
    
    ################################################
    # scale branch lengths of species tree model 2 #
    ################################################
    handle.SpeciesTree.Model2.i <- handle.SpeciesTree.Model2
    handle.SpeciesTree.Model2.i$edge.length[numeric.BranchToScale] <- handle.SpeciesTree.Model2.i$edge.length[numeric.BranchToScale] *numeric.Scale.i
    
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
    numeric.BHV <-dist.multiPhylo(x = c(handle.SpeciesTree.Model1, handle.SpeciesTree.Model2.i), method = "geodesic")
    
    ############################
    # Append to results matrix #
    ############################
    matrix.DistanceResults[i,] <- c(numeric.Hellinger, numeric.Kullback_Leibler, numeric.Jensen_Shannon, numeric.Robinson_Foulds, numeric.BHV)
    
  }  
  
  return(matrix.DistanceResults = matrix.DistanceResults)
}