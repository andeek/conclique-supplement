## A script to create the Swendsen-Wang computational and algorithmic timing results
## Depends on code in the helpers/ folder and 
## Results displayed in Table 2
## be sure packages in 0_libraries.R are installed prior to running
## written by: Andee Kaplan

## load libraries -----
library(conclique)

## Reproducibility ------------------------------------------
set.seed(1022)

## Swendsen Wang Functions -----------------------------------------------------------
# utility function to get neighbor pairs from neighbor list 
get_neighbor_pairs <- function(neighbors, eta) {
  # neighbors is a list containing dfs with neighboring vertices
  # eta is a vector with the same length
  stopifnot(length(eta) == length(neighbors))
  
  neighbor_pairs <- data.frame(v1 = NULL, v2 = NULL, eta = NULL)
  for(i in seq_along(neighbors)) {
    neigh <- neighbors[[i]]
    for(j in seq_len(ncol(neigh) - 1) + 1) {
      neighbor_pairs <- rbind(neighbor_pairs, data.frame(v1 = neigh[, 1], v2 = neigh[, j], eta = eta[i]))
    }
  }
  return(neighbor_pairs[neighbor_pairs$v1 < neighbor_pairs$v2, ])
}

# utility function for the swendson wang external field values
external_field_fn <- function(kappa, eta_a) {
  # y is a vector of the actual binary values for the field
  # kappa is either a numeric value or a vector of values of length equal to y
  # eta_a is the adjacency matrix times eta_ij
  stopifnot(length(kappa) == 1 | length(kappa) == dim(eta_a)[1])
  
  if(length(kappa) == 1) kappa <- rep(kappa, dim(eta_a)[1])
  
  log(kappa/(1 - kappa)) - eta_a %*% matrix(kappa - 1/2, ncol = 1) ## centered
}

# swendson wang algorithm 
sw <- function(neighbors, params, inits, n_iter = 100) {
  # neighbors is the result of calling conclique::get_neighbors()
  # params is a named list with the appropriate parameter values
  # init is the initialization values for the lattice
  # n_iter is the number of iterations to run
  
  # neighbor_pairs is a df with three columns
  # first two columns are neighboring locations
  # column 3 is the eta that applies to that location pair
  neighbor_pairs <- get_neighbor_pairs(neighbors, params$eta)
  
  # adjacency from neigbor pairs
  # scaled by eta for each pair of locations
  n <- length(inits)
  eta_a <- get_eta_a_cpp(as.matrix(neighbor_pairs), n)
  
  # store the results in dat
  dat <- matrix(NA, nrow = n_iter + 1, ncol = n)
  dat[1,] <- as.vector(inits)
  
  # external field values
  ext_field <- external_field_fn(params$kappa, eta_a)
  
  # get the samples
  for(i in seq_len(n_iter) + 1) {
    z <- dat[i - 1,]
    dat[i, ] <- sw_sampler(z, neighbor_pairs, params, ext_field, eta_a, n)
  }
  return(dat[-1,])
}

# inner sampler function 
sw_sampler <- function(z, neighbor_pairs, params, ext_field, eta_a, n){
  # get upper_bound for each pair of neighbors
  upper_bound <- exp(neighbor_pairs[,3]/2 * (z[neighbor_pairs[, 1]] == z[neighbor_pairs[, 2]]))
  
  # sample u uniformly, one for each pair of neighbors
  u <- vapply(upper_bound, function(x) runif(1, 0, x), FUN.VALUE = numeric(1))
  
  # is u > 1? these are the bonds
  bonds <- u > 1
  
  # cluster
  clusters <- as.numeric(get_cluster_labels_cpp(as.matrix(neighbor_pairs), bonds, n))
  clust_list <- as.list(unique(clusters))
  names(clust_list) <- unique(clusters)
  
  # compute p for each cluster
  p <- lapply(clust_list, function(clust) 1/(1 + exp(-sum(ext_field[clusters == clust]))))
  
  # sample bern(p) for each cluster
  x <- lapply(p, function(i) rbinom(1, 1, i))
  
  # return values 
  return(unlist(x[as.character(clusters)]))
}

## source additional sw functs 
Rcpp::sourceCpp("helpers/sw_functs.cpp")

#create grid and neighbors
N <- 40
grid <- matrix(1:(N*N), nrow = N)
lattice <- lattice_4nn_torus(dimvec = c(N, N))
neighbors <- get_neighbors(lattice, TRUE, grid)

## Compute SW Algorithmic Efficiency ---------------------------------
# parameters formatted for sw
par0 <- list(eta = c(.5, .5), kappa = .5)
par1 <- list(eta = c(.2, .7), kappa = .5)
reg <- cbind(1, as.vector(col(grid))) %*% matrix(c(-10, .5), ncol = 1)
par2 <- list(eta = c(.2, .7), kappa = exp(reg)/(1 + exp(reg)))
# make the ising model with no external field and close to the critical point
par3 <- list(eta = c(.88, .88), kappa = uniroot(function(k, eta) log(k/(1-k)) - 4*eta*(k - 1/2), interval = c(0, 1), eta = .88, tol = 10e-10)$root)

## repeat 
times <- 10
num_sample <- 10000

## prep for timing 
n <- length(grid)
neighbor_pairs0 <- get_neighbor_pairs(neighbors, par0$eta)
neighbor_pairs1 <- get_neighbor_pairs(neighbors, par1$eta)
neighbor_pairs2 <- get_neighbor_pairs(neighbors, par2$eta)
neighbor_pairs3 <- get_neighbor_pairs(neighbors, par3$eta)
eta_a0 <- get_eta_a_cpp(as.matrix(neighbor_pairs0), n)
eta_a1 <- get_eta_a_cpp(as.matrix(neighbor_pairs1), n)
eta_a2 <- get_eta_a_cpp(as.matrix(neighbor_pairs2), n)
eta_a3 <- get_eta_a_cpp(as.matrix(neighbor_pairs3), n)

## storage 
algo_eff <- list(sw = list(m0 = rep(NA, times), m1 = rep(NA, times), m2 = rep(NA, times)), m3 = rep(NA, times))
samp_sw_m0 <- list()
samp_sw_m1 <- list()
samp_sw_m2 <- list()
samp_sw_m3 <- list()

## do the sampling 
for(i in seq_len(times)) {
  ## inits 
  inits <- matrix(rbinom(n = N*N, size = 1, prob = .5), nrow = N)
  
  ## samples 
  samp_sw_m0[[i]] <- sw(neighbors, par0, inits, num_sample)
  samp_sw_m1[[i]] <- sw(neighbors, par1, inits, num_sample)
  samp_sw_m2[[i]] <- sw(neighbors, par2, inits, num_sample)
  samp_sw_m3[[i]] <- sw(neighbors, par3, inits, num_sample)
  
  ## get A 
  algo_eff$sw$m0[i] <- min(apply(samp_sw_m0[[i]], 2, LaplacesDemon::IAT)^(-1))
  algo_eff$sw$m1[i] <- min(apply(samp_sw_m1[[i]], 2, LaplacesDemon::IAT)^(-1))
  algo_eff$sw$m2[i] <- min(apply(samp_sw_m2[[i]], 2, LaplacesDemon::IAT)^(-1), na.rm = TRUE)
  algo_eff$sw$m3[i] <- min(apply(samp_sw_m3[[i]], 2, LaplacesDemon::IAT)^(-1), na.rm = TRUE)
}

## Compute SW Computational Efficiency ---------------------------------
## only care about the sampler, not the up front costs
## do for 20000 samples, as with conclique/sequential (check this)
times_sw_m0 <- apply(do.call(rbind, samp_sw_m0[1:2]), 1, function(x) {
  ext_field <- external_field_fn(par0$kappa, eta_a0)
  time <- Sys.time()
  tmp <- sw_sampler(x, neighbor_pairs0, par0, ext_field, eta_a0, n)
  time2 <- Sys.time()
  as.numeric(difftime(time2, time, units = "secs"))
})
times_sw_m1 <- apply(do.call(rbind, samp_sw_m1[1:2]), 1, function(x) {
  ext_field <- external_field_fn(par1$kappa, eta_a1)
  time <- Sys.time()
  tmp <- sw_sampler(x, neighbor_pairs1, par1, ext_field, eta_a1, n)
  time2 <- Sys.time()
  as.numeric(difftime(time2, time, units = "secs"))
})
times_sw_m2 <- apply(do.call(rbind, samp_sw_m2[1:2]), 1, function(x) {
  ext_field <- external_field_fn(par2$kappa, eta_a2)
  time <- Sys.time()
  tmp <- sw_sampler(x, neighbor_pairs2, par2, ext_field, eta_a2, n)
  time2 <- Sys.time()
  as.numeric(difftime(time2, time, units = "secs"))
})
times_sw_m3 <- apply(do.call(rbind, samp_sw_m3[1:2]), 1, function(x) {
  ext_field <- external_field_fn(par3$kappa, eta_a3)
  time <- Sys.time()
  tmp <- sw_sampler(x, neighbor_pairs3, par3, ext_field, eta_a3, n)
  time2 <- Sys.time()
  as.numeric(difftime(time2, time, units = "secs"))
})