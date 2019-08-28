// helper functions for performing Swendsen-Wang algorithm
// written by: Andee Kaplan

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp; using namespace arma;

// utility function for clustering -- join groups ----
arma::uvec union_cpp(arma::uvec labels, int x, int y) {
  uvec res = labels;
  res(y) = labels(x);
  return(res);
}

// utility function to get cluster labels ----
// Hoshen-Kopelman Algorithm
// https://www.ocf.berkeley.edu/~fricke/projects/hoshenkopelman/hoshenkopelman.html
// [[Rcpp::export]]  
arma::uvec get_cluster_labels_cpp(arma::mat neighbor_pairs, arma::uvec bonds, int n) {
  // setup variables
  uvec clusters(n);
  int largest_label = 0;
  
  for(int i = 0; i < n; ++i) {
    // for each point, find it's bonded neighbors that are "up" and "left"
    // neighbor pairs is such that v1 < v2 always, i.e. v1 will always be up or left of v2
    uvec idx = find((neighbor_pairs.col(1) == (i + 1)) && (bonds == 1));
    mat up_left_full = neighbor_pairs.rows(idx);
    uvec up_left = conv_to<uvec>::from(up_left_full.col(0)) - 1;
    
    if(up_left.n_elem == 0) { // Neither a neighbor above nor to the left
      ++largest_label; // Make a new, as-yet-unused cluster label
      clusters(i) = largest_label;
    } else { // neighbors either left or above
      if(up_left.n_elem > 1) {
        // link them
        for(int j = 1; j < up_left.n_elem; ++j) {
          clusters = union_cpp(clusters, up_left(0), up_left(j));
        }
      }
      // find label
      clusters(i) = clusters(up_left(0));
    }
  }
  return(clusters);
}  

// [[Rcpp::export]]  
arma::mat get_eta_a_cpp(arma::mat neighbor_pairs, int n) {
  mat eta_a(n, n, fill::zeros);
  for(int i = 0; i < n; ++i) {
    for(int j  = 0; j < n; ++j) {
      if(i <= j) {
        uvec match = (neighbor_pairs.col(0) == (i + 1)) && (neighbor_pairs.col(1) == (j + 1));
        if(sum(match) > 0) {
          uvec match_idx = find(match == 1); 
          mat eta_row = neighbor_pairs.rows(match_idx); //should be only 1
          eta_a(i, j) = eta_row(2);
          eta_a(j, i) = eta_row(2);
        }
      }
    }
  }
  return(eta_a);
}