## A script to cinstall and load all necessary libraries
## written by: Andee Kaplan

# Load and install libraries ----
if(!require(knitr)) install.packages("knitr"); library(knitr)
if(!require(ggplot2)) install.packages("ggplot2"); library(ggplot2)
if(!require(dplyr)) install.packages("dplyr"); library(dplyr)
if(!require(tidyr)) install.packages("tidyr"); library(tidyr)
if(!require(scales)) install.packages("scales"); library(scales)
if(!require(rootSolve)) install.packages("rootSolve"); library(rootSolve)
if(!require(agridat)) install.packages("agridat"); library(agridat)
if(!require(devtools)) install.packages("devtools")
if(!require(conclique)) devtools::install_github("andeek/conclique"); library(conclique)
if(!require(xtable)) install.packages("xtable"); library(xtable)
if(!require(Rcpp)) install.packages("Rcpp")
if(!require(RcppArmadillo)) install.packages("RcppArmadillo") 
if(!require(igraph)) install.packages("igraph") 
if(!require(tidyverse)) install.packages("tidyverse") 
if(!require(purrr)) install.packages("purrr")
if(!require(purrrlyr)) install.packages("purrrlyr") 

