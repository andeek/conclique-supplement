# Data for "Simulating Markov random fields with a conclique-based Gibbs sampler"

This is the data used in "Simulating Markov random fields with a conclique-based Gibbs sampler", a paper published in the Journal of Computational and Graphical Statistics.

## Authors

Andee Kaplan (andee.kaplan@colostate.edu), Colorado State University  
Mark S. Kaiser, Iowa State University  
Soumendra N. Lahiri, Washington University in St. Louis  
Daniel J. Nordman, Iowa State University 

## Contents

We include two data files that correspond to the real data example in Section 5.2. The raw data can be obtained from the Political Analysis Dataverse (https://dataverse.harvard.edu/dataverse/pan), and the citation for its use is

Chyzh, Olga V, and Mark S Kaiser. 2019. "A Local Structure Graph Model: Modeling Formation of Network Edges as a Function of Other Edges." Political Analysis In Press.

### edgedat.txt

This is a file with columns `ccode1`, `ccode2`, `edge`, `defense`, `mil_ratio`, `tot_trade`, `joint_dem`, `year`, `edgelab`, and `nis`
 
`ccode1` and `ccode2` are original country codes from Chyzh and Kaiser (2019) and `edge` is the original edge from the paper  
`defense` is the response (0/1) for presence of edge  
`mil_ratio` is covariate military size ratio  
`tot_trade` is covariate total reciprocal trade (needs to be transformed log(x+1) for model)  
`joint_dem` is covariate for joint democracy (0/1 values but not the response)  
`year` is year 1947 to 2007  
`edgelab` is my new edge index 1 through 45,513  
`nis` is size of neighborhood for that edge  

### nbhds.txt

This is neighborhood matrix size 45,513 by 554.  Each row is a potential edge with the edge index in column 1 (1 through 45,513).  The other entries in the row are indices of that edges neighbors or 0s to fill in the matrix.
 
For example, row 1 is for edge 1 which has 8 neighbors and then zeros to fill in the row. The edge with the most neighbors has 553 neighbors.

