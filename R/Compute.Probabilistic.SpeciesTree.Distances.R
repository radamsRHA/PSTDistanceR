#' Compute.Probabilistic.SpeciesTree.Distances: function to compute probabilistic species tree distances between two species tree models
#'
#' This function returns a vector of length (3): representing the "Hellinger", "Kullback_Leibler", and "Jensen_Shannon" distances between to trees
#' @param handle.SpeciesTree.Model1 Species Tree Model 1
#' @param handle.SpeciesTree.Model2 Species Tree Model 2
#' @param string.PathParentDir Path to a parent directory that will be used for subdirectories
#' @param string.PathHybridCoal Path to HYBRID-COAL executable
#' @keywords multispecies coalescent, model distance, incomplete lineage sorting, hybridization
#' @return vector.Probabilistic.SpeciesTree.Distances Vector containing the three PSTD measures\cr \cr 
#' @export
#' @examples
#' 
#' string.SpeciesTree.1a <- "(((A:1,B:1):1,C:2):1,D:3);"
#' handle.SpeciesTree.1a <- read.tree(text = string.SpeciesTree.1a)
#' 
#' string.SpeciesTree.1b <- "(((A:1,C:1):1,B:2):1,D:3);"
#' handle.SpeciesTree.1b <- read.tree(text = string.SpeciesTree.1b)
#' 
#' 
#' Compute.Probabilistic.SpeciesTree.Distances(handle.SpeciesTree.Model1 = handle.SpeciesTree.1a, 
#'    handle.SpeciesTree.Model2 = handle.SpeciesTree.1b, 
#'    string.PathParentDir = '~/Desktop/EXAMPLE/SmallTree/', 
#'    string.PathHybridCoal = '/Applications/hybrid-coal-v0.2.1-beta/hybrid-coal')

###############################################
# Compute.Probabilistic.SpeciesTree.Distances #
###############################################
Compute.Probabilistic.SpeciesTree.Distances <- function(handle.SpeciesTree.Model1, handle.SpeciesTree.Model2, string.PathParentDir, string.PathHybridCoal){
  
  ##############################################
  # Create top-level subdirectory for analysis #
  ##############################################
  string.ChildDir0 = paste(string.PathParentDir, '/PSTD_', Sys.Date(), sep = "")
  unlink(string.ChildDir0, recursive = T)
  dir.create(string.ChildDir0, showWarnings = T, recursive = T)
  
  ###################
  # Run HYBRID-COAL #
  ###################
  vector.GeneTreeProbs.Model.1 <- RunHybridCoal(string.PathHybridCoal = string.PathHybridCoal, 
                                                handle.SpeciesTree = handle.SpeciesTree.Model1, 
                                                string.PathParentDir = string.ChildDir0, 
                                                string.SpeciesTree.ModelName = "Model1")
  
  vector.GeneTreeProbs.Model.2 <- RunHybridCoal(string.PathHybridCoal = string.PathHybridCoal, 
                                                handle.SpeciesTree = handle.SpeciesTree.Model2, 
                                                string.PathParentDir = string.ChildDir0, 
                                                string.SpeciesTree.ModelName = "Model2")
  ###############
  # Compute PST #
  ###############
  vector.Probabilistic.SpeciesTree.Distances <- ComputePSTD(vectorM1 = vector.GeneTreeProbs.Model.1, vectorM2 = vector.GeneTreeProbs.Model.2)
  names(vector.Probabilistic.SpeciesTree.Distances) <- c("Hellinger", "Kullback_Leibler", "Jensen_Shannon")
  
  return(vector.Probabilistic.SpeciesTree.Distances = vector.Probabilistic.SpeciesTree.Distances)
}