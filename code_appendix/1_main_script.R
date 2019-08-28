## A script to create all the figures in the paper and generate all results
## Depends on code in the helpers/ folder and 
## stored .RData objects (saved from other helper scripts) in the written_results/ folder
## written by: Andee Kaplan

## Reproducibility ------------------------------------------
set.seed(1022)

## Graph concliques plot, Figure 2----
source("helpers/get_concliques.R")

conc <- edge_concliques(6)
g <- make_full_graph(6)

conc_color <- data.frame(edge = unlist(conc), conc = rep(1:length(conc), each = length(conc[[1]]))) %>%
  arrange(edge) %>%
  select(conc) %>%
  as.matrix()

colors <- c("#252525", "#525252", "#737373", "969696", "#bdbdbd", "#d9d9d9", "#f0f0f0")

lty <- 1:6

plot(g, 
     vertex.color = "white", vertex.frame.color="grey30", vertex.size=40, vertex.label.color = "grey30",
     edge.color = colors[conc_color], edge.lty = lty[conc_color], edge.width=2)


## ----4nn Timing example---------------------------------------------------------
## results saved from helpers/4nn_timing.R
load("written_results/4nn_timing.RData")

summary_times <- timings %>% 
  filter(N == 75) %>% 
  gather(gibbs, time, conclique, sequential) %>% 
  group_by(gibbs, n.iter) %>% 
  summarise(mean_time = mean(time)) %>%
  spread(gibbs, mean_time)


## 4nn Timing plot, Figure 3 --------
timings %>%
  gather(gibbs, time, conclique, sequential) %>% ungroup() %>%
  mutate(n.iter_f = factor(paste("M =", n.iter), levels=unique(paste("M =", timings$n.iter)))) %>%
  group_by(N, n.iter, gibbs) %>%
  mutate(mean_time = mean(time)) %>%
  ggplot() +
  geom_jitter(aes(N, log(time), shape = gibbs), alpha = .1, size = 2) +
  geom_line(aes(N, log(mean_time), lty = gibbs), size = 1) +
  facet_wrap(~n.iter_f, nrow = 1) +
  xlab("m") +
  ylab("Log Time (seconds)") +
  scale_shape_discrete("Gibbs sampler", labels=c("Conclique", "Single-site")) +
  scale_linetype_discrete("Gibbs sampler", labels=c("Conclique", "Single-site")) +
  theme(legend.position = "bottom") +
  theme(panel.spacing = unit(2, "lines"))


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
                 seq = list(m0 = rep(NA, times), m1 = rep(NA, times), m2 = rep(NA, times), m3 = rep(NA, times)),
                 sw = list(m0 = rep(NA, times), m1 = rep(NA, times), m2 = rep(NA, times)), m3 = rep(NA, times))

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

## Compute SW Algorithmic Efficiency ---------------------------------
# parameters formatted for sw
par0 <- list(eta = c(.5, .5), kappa = .5)
par1 <- list(eta = c(.2, .7), kappa = .5)
reg <- cbind(1, as.vector(col(grid))) %*% matrix(c(-10, .5), ncol = 1)
par2 <- list(eta = c(.2, .7), kappa = exp(reg)/(1 + exp(reg)))
# make the ising model with no external field and close to the critical point
par3 <- list(eta = c(.88, .88), kappa = uniroot(function(k, eta) log(k/(1-k)) - 4*eta*(k - 1/2), interval = c(0, 1), eta = .88, tol = 10e-10)$root)


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

## Compute SW Cimputational Efficiency ---------------------------------
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


## Graph example timings ----
## results saved from helpers/graph_timing.R
load("written_results/graph_timing.RData")

summary_times2 <- timings %>% 
  filter(n == 100) %>% 
  gather(gibbs, time, conclique, sequential) %>% 
  group_by(gibbs) %>% 
  summarise(mean_time = mean(time)) %>%
  spread(gibbs, mean_time)

summary_times200 <- timings %>% 
  filter(n == 200) %>% 
  summarise(mean_time = mean(conclique))

timings %>%
  gather(gibbs, time, conclique, sequential) %>% ungroup() %>%
  group_by(n, gibbs) %>%
  mutate(mean_time = mean(time)) %>%
  ggplot() +
  geom_point(aes(n, time, shape = gibbs), alpha = .2, size = 2, position = position_jitter(w = 0.5, h = 0)) +
  geom_line(aes(n, mean_time, lty = gibbs), size = 1) +
  xlab("V") +
  ylab("Time (seconds)") +
  scale_shape_discrete("Gibbs sampler", labels=c("Conclique", "Single-site")) +
  scale_linetype_discrete("Gibbs sampler", labels=c("Conclique", "Single-site"))

