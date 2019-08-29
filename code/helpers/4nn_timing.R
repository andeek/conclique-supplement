## helper script to generate written_results/4nn_timing.RData
## records timing results for simulating from a spatial 4nn Gaussian MRF
## grid sizes: 5x5, 10x10, 20x20, 30x30, 50x50, 75x75
## number of iterations: 100, 1000, 5000, 10000
## written by: Andee Kaplan

#load libraries ------------------------------------
library(igraph)
library(dplyr)
library(conclique) #devtools::install_github("andeek/conclique")
library(tidyr)
library(purrr)

#setup params for simulation ----------------------
Ns <- c(5, 10, 20, 30, 50, 75)
n_iter <- c(100, 1000, 5000, 10000)
params <- tbl_df(expand.grid(N = Ns,
                             eta = .24, rho = 1, kappa = 0)) 

#create concliques list -------------------------------
concliques_4nn <- function(grid) {
  concliques <- list()
  concliques[[1]] <- grid[(row(grid) %% 2 == 1 & col(grid) %% 2 == 1) | (row(grid) %% 2 == 0 & col(grid) %% 2 == 0)]
  concliques[[2]] <- grid[(row(grid) %% 2 == 0 & col(grid) %% 2 == 1) | (row(grid) %% 2 == 1 & col(grid) %% 2 == 0)]
  class(concliques) <- "conclique_cover"
  return(concliques)
}

#create grid, lattices, neighbors, concliques for MCMC ------------------
params %>%
  mutate(grid = map(N, function(N) matrix(1:(N*N), nrow = N))) %>%
  mutate(conclique = map(grid, concliques_4nn)) %>%
  mutate(lattice = map(N, function(N) lattice_4nn_torus(dimvec = c(N, N)))) %>%
  mutate(inits = map(N, function(N) matrix(0, nrow = N, ncol = N))) %>%
  mutate(neighbors = map(lattice, get_neighbors)) -> lattices

lattices %>%
  left_join(expand.grid(N = Ns,
                        n.iter = n_iter)) -> params

#run mcmc for various param values ---------------------------
times <- 10

#timing functions
time_conclique <- function(concliques, neighbors, inits, sampler, params, n.iter, times) {
  res <- rep(0, times)
  for(i in 1:times) {
    time <- Sys.time()
    run_conclique_gibbs(concliques, neighbors, inits, sampler, params, n.iter)
    res[i] <- as.numeric(difftime(Sys.time(), time, units = "secs"))
  }
  return(res)
}

time_sequential <- function(neighbors, inits, sampler, params, n.iter, times) {
  res <- rep(0, times)
  for(i in 1:times) {
    time <- Sys.time()
    run_sequential_gibbs(neighbors, inits, sampler, params, n.iter)
    res[i] <- as.numeric(difftime(Sys.time(), time, units = "secs"))
    
  }
  return(res)
}

params %>%
  group_by(N, eta, rho, kappa, n.iter) %>%
  do(data.frame(conclique = time_conclique(.$conclique[[1]], .$neighbors[[1]], .$inits[[1]], "gaussian_single_param", list(eta = .$eta, kappa = .$kappa, rho= .$rho), .$n.iter, times),
     sequential = time_sequential(.$neighbors[[1]], .$inits[[1]], "gaussian_single_param", list(eta = .$eta, kappa = .$kappa, rho= .$rho), .$n.iter, times))) -> timings

save(timings, file = "written_results/4nn_timing.RData")