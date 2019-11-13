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
