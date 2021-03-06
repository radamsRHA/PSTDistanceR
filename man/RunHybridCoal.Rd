% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RunHybridCoal.R
\name{RunHybridCoal}
\alias{RunHybridCoal}
\title{RunHybridCoal: function to run HYBRID-COAL on a given species tree model, and return a vector of gene tree probabilities for the species tree}
\usage{
RunHybridCoal(string.PathHybridCoal, handle.SpeciesTree,
  string.PathParentDir, string.SpeciesTree.ModelName)
}
\arguments{
\item{string.PathHybridCoal}{Path to HYBRID-COAL executable}

\item{handle.SpeciesTree}{Species tree}

\item{string.PathParentDir}{Path to a parent directory that will be used for subdirectories}

\item{string.SpeciesTree.ModelName}{String giving the name of the species model}
}
\value{
vector.GeneTreeProbabilities: vector containing a set of gene tree probabilities for the species tree model\cr
}
\description{
This function returns vector.GeneTreeProbabilities: a vector of gene tree probabilities for the species model (sum ~1.0)
}
\examples{

string.SpeciesTree.1a <- "(((A:1,B:1):1,C:2):1,D:3);"
handle.SpeciesTree.1a <- read.tree(text = string.SpeciesTree.1a)

RunHybridCoal(string.PathHybridCoal = '/Applications/hybrid-coal-v0.2.1-beta/hybrid-coal', 
   handle.SpeciesTree = handle.SpeciesTree.1a, 
   string.PathParentDir = '~/Desktop/EXAMPLE/SmallTree/', 
   string.SpeciesTree.ModelName = "Model1a")
   
}
\keyword{coalescent,}
\keyword{distance,}
\keyword{hybridization}
\keyword{incomplete}
\keyword{lineage}
\keyword{model}
\keyword{multispecies}
\keyword{sorting,}
