#' Read.GeneTree.Probabilities: function to extract the gene tree probabilities from the output of HYBRID-COAL
#'
#' This function returns a vector of containing gene tree probabilities (assumed to be in the same order when comparing two gene tree probability vectors)
#' @param string.PathToGeneTreeProbs Vector containing a set of likelihoods under model 1
#' @keywords multispecies coalescent, model distance, incomplete lineage sorting, hybridization
#' @return vector.GeneTreeProbabilities: vector containing a set of gene tree probabilities\cr 
#' @export
#' @examples
#' vector.GeneTreeProbs.M1 <- Read.GeneTree.Probabilities(string.PathToGeneTreeProbs = '~/Desktop/EXAMPLE/SmallTree/M1.prob')
#' 

###############################
# Read.GeneTree.Probabilities #
###############################
Read.GeneTree.Probabilities <- function(string.PathToGeneTreeProbs){
  
  ###############
  # Read infile #
  ###############
  handle.GeneTreeProbFile <- readLines(string.PathToGeneTreeProbs)
  numeric.NumberOfGeneTrees <- length(handle.GeneTreeProbFile)
  vector.GeneTreeProbabilities <- rep(NA, numeric.NumberOfGeneTrees)
  
  #########################################
  # Loop through gene trees and get probs #
  #########################################
  for (i in 1:numeric.NumberOfGeneTrees){
    
    #####################
    # read line by line #
    #####################
    string.GeneTreeProb.Line <- handle.GeneTreeProbFile[i]
    handle.StringSplit.GeneTreeProb.Line <- strsplit(x = string.GeneTreeProb.Line, split = '\t')[[1]]
    vector.GeneTreeProbabilities[i] <- as.numeric(handle.StringSplit.GeneTreeProb.Line[2])
    
  }
  return(vector.GeneTreeProbabilities = vector.GeneTreeProbabilities)
}