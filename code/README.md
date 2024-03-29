# Code appendix for "Simulating Markov random fields with a conclique-based Gibbs sampler"

This is the code appendix accompanying "Simulating Markov random fields with a conclique-based Gibbs sampler", a paper published in the Journal of Computational and Graphical Statistics.

## Authors

Andee Kaplan (andee.kaplan@colostate.edu), Colorado State University  
Mark S. Kaiser, Iowa State University  
Soumendra N. Lahiri, Washington University in St. Louis  
Daniel J. Nordman, Iowa State University 

## Contents

This folder contains all necessary R and C++ scripts for generating the figures in the paper. Code in the folder is written such that the folder is the working directory and everything is run relative to that. Contents include:

1. libraries.R -- A script to install and load all necessary libraries in R.
2. plot.R -- A script to create all the figures in the paper. Depends on code in the helpers/ folder and stored .RData objects (saved from other helper scripts) in the written_results/ folder. Some results have been cached because the amount of time necessary to run these scripts.
3. gibbs_efficiency.R -- A script to create all the Gibbs computational and algorithmic timing results.
4. sw_efficiency.R -- A script to create the Swendsen-Wang computational and algorithmic timing results.
5. dsatur_network_example.R -- A script to create the concliques using DSatur and perform the analysis of Section 5.2.
6. helpers/ -- A folder containing functions used in the main scripts as well as scripts for generating store results.
7. written_results/ -- A folder containing cached timing results for the 4nn Gaussian MRF and the network timing examples.
8. README.md -- this instructional README document.

