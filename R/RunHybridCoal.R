#' RunHybridCoal: function to run HYBRID-COAL on a given species tree model, and return a vector of gene tree probabilities for the species tree 
#'
#' This function returns vector.GeneTreeProbabilities: a vector of gene tree probabilities for the species model (sum ~1.0)
#' @param string.PathHybridCoal Path to HYBRID-COAL executable
#' @param handle.SpeciesTree Species tree
#' @param string.PathParentDir Path to a parent directory that will be used for subdirectories
#' @param string.SpeciesTree.ModelName String giving the name of the species model 
#' @keywords multispecies coalescent, model distance, incomplete lineage sorting, hybridization
#' @return vector.GeneTreeProbabilities: vector containing a set of gene tree probabilities for the species tree model\cr 
#' @export
#' @examples
#' 
#' string.SpeciesTree.1a <- "(((A:1,B:1):1,C:2):1,D:3);"
#' handle.SpeciesTree.1a <- read.tree(text = string.SpeciesTree.1a)
#' 
#' RunHybridCoal(string.PathHybridCoal = '/Applications/hybrid-coal-v0.2.1-beta/hybrid-coal', 
#'    handle.SpeciesTree = handle.SpeciesTree.1a, 
#'    string.PathParentDir = '~/Desktop/EXAMPLE/SmallTree/', 
#'    string.SpeciesTree.ModelName = "Model1a")
#'    

#################
# RunHybridCoal #
#################
RunHybridCoal <- function(string.PathHybridCoal, handle.SpeciesTree, string.PathParentDir, string.SpeciesTree.ModelName){
  
  #####################################################
  # Create subdirectory for analysis with hybrid-coal #
  #####################################################
  string.ChildDir0 = paste(string.PathParentDir, '/SpeciesTreeModel_', string.SpeciesTree.ModelName, '_', Sys.Date(), '/',sep = "")
  unlink(string.ChildDir0, recursive = T)
  dir.create(string.ChildDir0, showWarnings = T, recursive = T)
  
  #######################################
  # Set working dir to Parent directory #
  #######################################
  setwd(dir = string.ChildDir0)
  string.PathToSpeciesTree <- paste(string.ChildDir0, '/SpeciesTreeModel_', string.SpeciesTree.ModelName, '.tree',  sep = "")
  write.tree(phy = handle.SpeciesTree, file = string.PathToSpeciesTree)
  string.HybridCoal.Cmd <- gsub(pattern = "YYY", replacement = string.PathHybridCoal, x = "YYY -o XXX -sp ZZZ")
  string.HybridCoal.Cmd <- gsub(pattern = "XXX", replacement = string.SpeciesTree.ModelName, x = string.HybridCoal.Cmd)
  string.HybridCoal.Cmd <- gsub(pattern = "ZZZ", replacement = paste('./SpeciesTreeModel_', string.SpeciesTree.ModelName, '.tree', sep = ""), x = string.HybridCoal.Cmd)
  
  ###################
  # run HYBRID-COAL #
  ###################
  system(string.HybridCoal.Cmd)
  string.ToGeneTreeProbsFile <- paste(string.ChildDir0, '/', string.SpeciesTree.ModelName, '.prob', sep = "")
  
  ###################################
  # Extract gene tree probabilities #
  ###################################
  vector.GeneTreeProbabilities <- Read.GeneTree.Probabilities(string.PathToGeneTreeProbs = string.ToGeneTreeProbsFile)
  
  return(vector.GeneTreeProbabilities = vector.GeneTreeProbabilities)
  
}