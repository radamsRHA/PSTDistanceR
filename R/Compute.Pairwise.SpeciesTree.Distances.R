#' Compute.Pairwise.SpeciesTree.Distances: function to compute a matrix of pairwise distances among a set of species trees
#'
#' This function returns a list of pairwise distance matrix between all provided specie trees
#' @param list.SpeciesTrees List of input species trees
#' @param vector.SpeciesTreeNames Vector of species tree names
#' @param string.PathParentDir String for path to parent directory for analyses
#' @param string.PathHybridCoal String for path to hybrid-coal executable
#' @keywords multispecies coalescent, model distance, incomplete lineage sorting, hybridization
#' @return vector.Probabilistic.SpeciesTree.Distances Vector containing the three PSTD measures\cr \cr 
#' @export
#' @examples
#' 

##########################################
# Compute.Pairwise.SpeciesTree.Distances #
##########################################
Compute.Pairwise.SpeciesTree.Distances <- function(list.SpeciesTrees, vector.SpeciesTreeNames, string.PathParentDir, string.PathHybridCoal){
  
  ##########################################################
  # Specify matrices used for computing pairwise distances #
  ##########################################################
  numeric.NumberOfTrees <- length(list.SpeciesTrees)
  matrix.Hellinger_Distance <- matrix(nrow = numeric.NumberOfTrees, ncol = numeric.NumberOfTrees)
  matrix.Kullback_Leibler_Distance <- matrix(nrow = numeric.NumberOfTrees, ncol = numeric.NumberOfTrees)
  matrix.Jensen_Shannon_Distance <- matrix(nrow = numeric.NumberOfTrees, ncol = numeric.NumberOfTrees)
  matrix.RobinsonFoulds_Distance <- matrix(nrow = numeric.NumberOfTrees, ncol = numeric.NumberOfTrees)
  matrix.BHV_Distance <- matrix(nrow = numeric.NumberOfTrees, ncol = numeric.NumberOfTrees)
  rownames(matrix.Hellinger_Distance) <- colnames(matrix.Hellinger_Distance) <- rownames(matrix.Kullback_Leibler_Distance) <- colnames(matrix.Kullback_Leibler_Distance) <- rownames(matrix.Jensen_Shannon_Distance) <- colnames(matrix.Jensen_Shannon_Distance) <- rownames(matrix.RobinsonFoulds_Distance) <- colnames(matrix.RobinsonFoulds_Distance) <- rownames(matrix.BHV_Distance) <- colnames(matrix.BHV_Distance) <- vector.SpeciesTreeNames
  
  ########################################
  # Create parent directory for analyses #
  ########################################
  string.ChildDir0 = paste(string.PathParentDir, '/Pairwise_PSTD_Analyses', Sys.Date(), sep = "")
  unlink(string.ChildDir0, recursive = T)
  dir.create(string.ChildDir0, showWarnings = T, recursive = T)
  
  #############################
  # Conduct pairwise analyses #
  #############################
  for (i in 1:numeric.NumberOfTrees){
    
    ############################
    # set first species tree i #
    ############################
    handle.SpeciesTreeModel.i <- list.SpeciesTrees[[i]]
    string.SpeciesTree.i.Name <- vector.SpeciesTreeNames[i]
    
    ####################################
    # loop through second species tree #
    ####################################
    for (j in 1:numeric.NumberOfTrees){
      
      #############################
      # set second species tree j #
      #############################
      if (j != i){
        handle.SpeciesTreeModel.j <- list.SpeciesTrees[[j]]
        string.SpeciesTree.j.Name <- vector.SpeciesTreeNames[j]
        string.PairwiseDir <- paste("SpeciesTree_", string.SpeciesTree.i.Name, "__SpeciesTree_", string.SpeciesTree.j.Name, sep = "")
        string.Path.PairwiseDir <- paste(string.ChildDir0, '/', string.PairwiseDir, sep = "")
        handle.PSTD.RESULTS <- Compute.Probabilistic.SpeciesTree.Distances(handle.SpeciesTree.Model1 = handle.SpeciesTreeModel.i, 
                                                                           handle.SpeciesTree.Model2 = handle.SpeciesTreeModel.j, 
                                                                           string.PathParentDir = string.Path.PairwiseDir, 
                                                                           string.PathHybridCoal = string.PathHybridCoal)
      
        #####################
        # extract distances #
        #####################
        numeric.Hellinger <- handle.PSTD.RESULTS[1]
        numeric.Kullback_Leibler <- handle.PSTD.RESULTS[2]
        numeric.Jensen_Shannon <- handle.PSTD.RESULTS[3]
        numeric.RobinsonFoulds <- RF.dist(tree1 = handle.SpeciesTreeModel.i, tree2 = handle.SpeciesTreeModel.j, normalize = T, rooted = T)
        numeric.BHV <- dist.multiPhylo(x = c(handle.SpeciesTreeModel.i, handle.SpeciesTreeModel.j), method = "geodesic")
        
        #################################
        # Append to respective matrices #
        #################################
        matrix.Hellinger_Distance[i,j] <- numeric.Hellinger
        matrix.Kullback_Leibler_Distance[i,j] <- numeric.Kullback_Leibler
        matrix.Jensen_Shannon_Distance[i,j] <- numeric.Jensen_Shannon
        matrix.RobinsonFoulds_Distance[i,j] <- numeric.RobinsonFoulds
        matrix.BHV_Distance[i,j] <- numeric.BHV
        
      }
      if (j == i){
        matrix.Hellinger_Distance[i,j] <- 0
        matrix.Kullback_Leibler_Distance[i,j] <- 0
        matrix.Jensen_Shannon_Distance[i,j] <- 0
        matrix.RobinsonFoulds_Distance[i,j] <- 0
        matrix.BHV_Distance[i,j] <- 0
      }
    }
  }
  
  return(list(matrix.Hellinger_Distance = matrix.Hellinger_Distance, 
              matrix.Kullback_Leibler_Distance = matrix.Kullback_Leibler_Distance, 
              matrix.Jensen_Shannon_Distance = matrix.Jensen_Shannon_Distance, 
              matrix.RobinsonFoulds_Distance = matrix.RobinsonFoulds_Distance, 
              matrix.BHV_Distance = matrix.BHV_Distance))
}
