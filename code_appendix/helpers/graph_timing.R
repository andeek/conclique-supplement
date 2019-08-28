## helper script to generate written_results/graph_timing.RData
## records timing results for simulating from a network MRF
## graph sizes: 10, 30, 50, 100, 200
## written by: Andee Kaplan

# libraries -----
library(igraph)
library(tidyverse)
library(conclique)

source("helpers/get_concliques.R")
source("helpers/get_neighbors.R")

params <- expand.grid(n = c(10, 30, 50, 100, 200), 
                     eta = .2,
                     kappa = .5,
                     n_iter = 1000)

#create grid, lattices, neighbors, concliques for MCMC ------------------
params %>% 
  select(n, eta, kappa) %>%
  unique() %>%
  mutate(conclique = map(n, edge_concliques)) %>%
  mutate(neighbors = map(n, edge_neighbors)) %>%
  mutate(inits = map(n, function(n){ matrix(0, choose(n, 2)) })) -> mcmc_pieces

params %>%
  left_join(mcmc_pieces) -> params

times <- 10

#timing functions ----
time_conclique <- function(concliques, neighbors, inits, sampler, params, n_iter, times) {
  res <- rep(0, times)
  for(i in 1:times) {
    time <- Sys.time()
    run_conclique_gibbs(concliques, neighbors, inits, sampler, params, n_iter)
    res[i] <- as.numeric(difftime(Sys.time(), time, units = "secs"))
  }
  return(res)
}
time_sequential <- function(neighbors, inits, sampler, params, n_iter, times) {
  res <- rep(0, times)
  for(i in 1:times) {
    time <- Sys.time()
    run_sequential_gibbs(neighbors, inits, sampler, params, n_iter)
    res[i] <- as.numeric(difftime(Sys.time(), time, units = "secs"))
    
  }
  return(res)
}

params %>%
  group_by(n, eta, kappa, n_iter) %>%
  do(data.frame(conclique = time_conclique(.$conclique[[1]], .$neighbors[[1]], .$inits[[1]], "binary_single_param", list(eta = .$eta, kappa = .$kappa), .$n_iter, times),
                sequential = time_sequential(.$neighbors[[1]], .$inits[[1]], "binary_single_param", list(eta = .$eta, kappa = .$kappa), .$n_iter, times))) -> timings

save(timings, file = "written_results/graph_timing.RData")




