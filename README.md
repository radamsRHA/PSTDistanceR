# PSTDistanceR
PSTDistanceR: R package for computing probabilistic distances between species tree models under the multispecies coalscent

The R package PSTDistanceR is freely available to download and distribute from github <https://github.com/radamsRHA/PSTDistanceR/>. To install and load PSTDistanceR, you must first install the R packages `phangorn`, `distory`, `Rcpp`, and the program `hybrid-coal-v0.2.1-beta`. 

```
install.packages("devtools")
install.packages("distory")
install.packages("phangorn")
```
Now using devtools we can install `PSTDistanceR` from github:

```
library(devtools)
install_github("radamsRHA/PSTDistanceR")
library(PSTDistanceR) # Load package ThetaMater
```

Here's an example usage of PSTDistanceR for computing the distance between two species trees:

```
library(PSTDistanceR)
library(phangorn)

string.SpeciesTree.1a <- "(((A:1,B:1):1,C:2):1,D:3);"
handle.SpeciesTree.1a <- read.tree(text = string.SpeciesTree.1a)

string.SpeciesTree.1b <- "(((A:1,C:1):1,B:2):1,D:3);"
handle.SpeciesTree.1b <- read.tree(text = string.SpeciesTree.1b)


Compute.Probabilistic.SpeciesTree.Distances(handle.SpeciesTree.Model1 = handle.SpeciesTree.1a, 
   handle.SpeciesTree.Model2 = handle.SpeciesTree.1b, 
   string.PathParentDir = '~/Desktop/', 
   string.PathHybridCoal = '/Applications/hybrid-coal-v0.2.1-beta/hybrid-coal')

```


