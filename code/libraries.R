## A script to install all necessary libraries
## written by: Andee Kaplan

# Load and install libraries ----
if(!require(ggplot2)) install.packages("ggplot2"); 
if(!require(dplyr)) install.packages("dplyr"); 
if(!require(tidyr)) install.packages("tidyr"); 
if(!require(rootSolve)) install.packages("rootSolve");
if(!require(devtools)) install.packages("devtools")
if(!require(conclique)) devtools::install_github("andeek/conclique"); 
if(!require(Rcpp)) install.packages("Rcpp")
if(!require(RcppArmadillo)) install.packages("RcppArmadillo") 
if(!require(igraph)) install.packages("igraph") 
if(!require(tidyverse)) install.packages("tidyverse") 
if(!require(purrr)) install.packages("purrr")
if(!require(purrrlyr)) install.packages("purrrlyr") 
if(!require(LaplacesDemon)) install.packages("LaplacesDemon")

