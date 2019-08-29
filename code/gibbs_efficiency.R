## A script to create all the Gibbs computational and algorithmic timing results
## Results displayed in Table 2
## be sure packages in 0_libraries.R are installed prior to running
## written by: Andee Kaplan

## load libraries -----
library(dplyr)
library(conclique)

## Reproducibility ------------------------------------------
set.seed(1022)

## Algorithmic efficiency ------------------------------------------------
#create concliques list 
concliques_4nn <- function(grid) {
  concliques <- list()
  concliques[[1]] <- grid[(row(grid) %% 2 == 1 & col(grid) %% 2 == 1) | (row(grid) %% 2 == 0 & col(grid) %% 2 == 0)]
  concliques[[2]] <- grid[(row(grid) %% 2 == 0 & col(grid) %% 2 == 1) | (row(grid) %% 2 == 1 & col(grid) %% 2 == 0)]
  class(concliques) <- "conclique_cover"
  return(concliques)
}

N <- 40
params0 <- list(eta = .5, kappa = .5)
params1 <- list(eta_1 = .2, eta_2 = .7, kappa = .5)
params2 <- list(eta_1 = .2, eta_2 = .7, beta_0 = -10, beta_1 = .5)
# make the ising model with no external field and close to the critical point
params3 <- list(eta = .88, kappa = uniroot(function(k, eta) log(k/(1-k)) - 4*eta*(k - 1/2), interval = c(0, 1), eta = .88, tol = 10e-10)$root)

num_sample <- 10000

#create grid, lattices, neighbors, concliques for MCMC 
grid <- matrix(1:(N*N), nrow = N)
conclique <- concliques_4nn(grid)
lattice <- lattice_4nn_torus(dimvec = c(N, N))
neighbors <- get_neighbors(lattice, TRUE, grid)

## repeat 
times <- 10
algo_eff <- list(conc = list(m0 = rep(NA, times), m1 = rep(NA, times), m2 = rep(NA, times), m3 = rep(NA, times)),
                 seq = list(m0 = rep(NA, times), m1 = rep(NA, times), m2 = rep(NA, times), m3 = rep(NA, times)))

## do one for later 
inits <- matrix(rbinom(n = N*N, size = 1, prob = .5), nrow = N)
samp_conc_m0 <- run_conclique_gibbs(conclique, neighbors, inits, "binary_single_param", params0, num_sample)
samp_seq_m0 <- run_sequential_gibbs(neighbors, inits, "binary_single_param", params0, num_sample)
samp_conc_m1 <- run_conclique_gibbs(conclique, neighbors, inits, "binary_two_param", params1, num_sample)
samp_seq_m1 <- run_sequential_gibbs(neighbors, inits, "binary_two_param", params1, num_sample)
samp_conc_m2 <- run_conclique_gibbs(conclique, neighbors, inits, "binary_two_param_reg", params2, num_sample)
samp_seq_m2 <- run_sequential_gibbs(neighbors, inits, "binary_two_param_reg", params2, num_sample)
samp_conc_m3 <- run_conclique_gibbs(conclique, neighbors, inits, "binary_single_param", params3, num_sample)
samp_seq_m3 <- run_sequential_gibbs(neighbors, inits, "binary_single_param", params3, num_sample)


## get Alg
algo_eff$conc$m0[1] <- min(apply(samp_conc_m0, 2, LaplacesDemon::IAT)^(-1))
algo_eff$seq$m0[1] <- min(apply(samp_seq_m0, 2, LaplacesDemon::IAT)^(-1))
algo_eff$conc$m1[1] <- min(apply(samp_conc_m1, 2, LaplacesDemon::IAT)^(-1))
algo_eff$seq$m1[1] <- min(apply(samp_seq_m1, 2, LaplacesDemon::IAT)^(-1))
algo_eff$conc$m2[1] <- min(apply(samp_conc_m2, 2, LaplacesDemon::IAT)^(-1), na.rm = TRUE)
algo_eff$seq$m2[1] <- min(apply(samp_seq_m2, 2, LaplacesDemon::IAT)^(-1), na.rm = TRUE)
algo_eff$conc$m3[1] <- min(apply(samp_conc_m3, 2, LaplacesDemon::IAT)^(-1))
algo_eff$seq$m3[1] <- min(apply(samp_seq_m3, 2, LaplacesDemon::IAT)^(-1))

## do the rest 
for(i in 2:times) {
  ## get samples 
  inits <- matrix(rbinom(n = N*N, size = 1, prob = .5), nrow = N)
  samp_conc_m0 <- run_conclique_gibbs(conclique, neighbors, inits, "binary_single_param", params0, num_sample)
  samp_seq_m0 <- run_sequential_gibbs(neighbors, inits, "binary_single_param", params0, num_sample)
  samp_conc_m1 <- run_conclique_gibbs(conclique, neighbors, inits, "binary_two_param", params1, num_sample)
  samp_seq_m1 <- run_sequential_gibbs(neighbors, inits, "binary_two_param", params1, num_sample)
  samp_conc_m2 <- run_conclique_gibbs(conclique, neighbors, inits, "binary_two_param_reg", params2, num_sample)
  samp_seq_m2 <- run_sequential_gibbs(neighbors, inits, "binary_two_param_reg", params2, num_sample)
  samp_conc_m3 <- run_conclique_gibbs(conclique, neighbors, inits, "binary_single_param", params3, num_sample)
  samp_seq_m3 <- run_sequential_gibbs(neighbors, inits, "binary_single_param", params3, num_sample)

  ## get A 
  algo_eff$conc$m0[i] <- min(apply(samp_conc_m0, 2, LaplacesDemon::IAT)^(-1))
  algo_eff$seq$m0[i] <- min(apply(samp_seq_m0, 2, LaplacesDemon::IAT)^(-1))
  algo_eff$conc$m1[i] <- min(apply(samp_conc_m1, 2, LaplacesDemon::IAT)^(-1))
  algo_eff$seq$m1[i] <- min(apply(samp_seq_m1, 2, LaplacesDemon::IAT)^(-1))
  algo_eff$conc$m2[i] <- min(apply(samp_conc_m2, 2, LaplacesDemon::IAT)^(-1), na.rm = TRUE)
  algo_eff$seq$m2[i] <- min(apply(samp_seq_m2, 2, LaplacesDemon::IAT)^(-1), na.rm = TRUE)
  algo_eff$conc$m3[i] <- min(apply(samp_conc_m3, 2, LaplacesDemon::IAT)^(-1))
  algo_eff$seq$m3[i] <- min(apply(samp_seq_m3, 2, LaplacesDemon::IAT)^(-1))
}

## Computational efficiency ----------------------------------------
## formatting functions
get_data <- function(data, neighbors, gibbs = c("conclique", "sequential"), concliques = NULL, N = NULL) {
  stopifnot((gibbs == "conclique" & !is.null(concliques)) | gibbs == "sequential")
  
  nums <- lapply(neighbors, function(neigh) {
    rowSums(!is.na(neigh[, -1]))
  }) 
  sums <- lapply(neighbors, function(neigh) {
    rowSums(matrix(data[neigh[, -1]], ncol = ncol(neigh) - 1, byrow = TRUE))
  }) 
  u <- (0:(length(data) - 1)) %% N + 1
  
  
  if(gibbs == "sequential") {
    res <- lapply(seq_along(data), function(i) list(data = data[i], nums = lapply(nums, function(num) num[i]),
                                                    sums = lapply(sums, function(sum) sum[i]), u = u[i]))
  } else {
    res <- lapply(concliques, function(conc) {
      list(data = data[conc], 
           nums = lapply(nums, function(num) num[conc]),
           sums = lapply(sums, function(sum) sum[conc]), 
           u = u[conc])
    })
  }
  return(res)
}

format_conc_m0 <- lapply(seq_len(nrow(samp_conc_m0)), function(row) {
  get_data(samp_conc_m0[row, ], neighbors, "conclique", conclique, N)
}) %>%
  unlist(recursive = FALSE)
format_seq_m0 <- lapply(1:10, function(row) {
  get_data(samp_seq_m0[row, ], neighbors, "sequential", conclique, N)
}) %>%
  unlist(recursive = FALSE)
format_conc_m1 <- lapply(seq_len(nrow(samp_conc_m1)), function(row) {
  get_data(samp_conc_m1[row, ], neighbors, "conclique", conclique, N)
}) %>%
  unlist(recursive = FALSE)
format_seq_m1 <- lapply(1:10, function(row) {
  get_data(samp_seq_m1[row, ], neighbors, "sequential", conclique, N)
}) %>%
  unlist(recursive = FALSE)
format_conc_m2 <- lapply(seq_len(nrow(samp_conc_m2)), function(row) {
  get_data(samp_conc_m2[row, ], neighbors, "conclique", conclique, N)
}) %>%
  unlist(recursive = FALSE)
format_seq_m2 <- lapply(1:10, function(row) {
  get_data(samp_seq_m2[row, ], neighbors, "sequential", conclique, N)
}) %>%
  unlist(recursive = FALSE)
format_conc_m3 <- lapply(seq_len(nrow(samp_conc_m3)), function(row) {
  get_data(samp_conc_m3[row, ], neighbors, "conclique", conclique, N)
}) %>%
  unlist(recursive = FALSE)
format_seq_m3 <- lapply(1:10, function(row) {
  get_data(samp_seq_m3[row, ], neighbors, "sequential", conclique, N)
}) %>%
  unlist(recursive = FALSE)

## Get computational times
times_conc_m0 <- lapply(format_conc_m0, function(x) {
  time <- Sys.time()
  tmp <- binary_single_param_sampler(x, params0)
  time2 <- Sys.time()
  as.numeric(difftime(time2, time, units = "secs"))
}) %>% unlist()
times_seq_m0 <- lapply(format_seq_m0, function(x) {
  time <- Sys.time()
  tmp <- binary_single_param_sampler(x, params0)
  time2 <- Sys.time()
  as.numeric(difftime(time2, time, units = "secs"))
}) %>% unlist()
times_conc_m1 <- lapply(format_conc_m1, function(x) {
  time <- Sys.time()
  tmp <- binary_two_param_sampler(x, params1)
  time2 <- Sys.time()
  as.numeric(difftime(time2, time, units = "secs"))
}) %>% unlist()
times_seq_m1 <- lapply(format_seq_m1, function(x) {
  time <- Sys.time()
  tmp <- binary_two_param_sampler(x, params1)
  time2 <- Sys.time()
  as.numeric(difftime(time2, time, units = "secs"))
}) %>% unlist()
times_conc_m2 <- lapply(format_conc_m2, function(x) {
  time <- Sys.time()
  tmp <- binary_two_param_reg_sampler(x, params2)
  time2 <- Sys.time()
  as.numeric(difftime(time2, time, units = "secs"))
}) %>% unlist()
times_seq_m2 <- lapply(format_seq_m2, function(x) {
  time <- Sys.time()
  tmp <- binary_two_param_reg_sampler(x, params2)
  time2 <- Sys.time()
  as.numeric(difftime(time2, time, units = "secs"))
}) %>% unlist()
times_conc_m3 <- lapply(format_conc_m3, function(x) {
  time <- Sys.time()
  tmp <- binary_single_param_sampler(x, params3)
  time2 <- Sys.time()
  as.numeric(difftime(time2, time, units = "secs"))
}) %>% unlist()
times_seq_m3 <- lapply(format_seq_m3, function(x) {
  time <- Sys.time()
  tmp <- binary_single_param_sampler(x, params3)
  time2 <- Sys.time()
  as.numeric(difftime(time2, time, units = "secs"))
}) %>% unlist()
