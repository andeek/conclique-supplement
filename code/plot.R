## A script to create all the figures in the paper
## Depends on code in the helpers/ folder and 
## stored .RData objects (saved from other helper scripts) in the written_results/ folder
## be sure packages in 0_libraries.R are installed prior to running
## written by: Andee Kaplan

## load libraries -----
library(ggplot2)
library(dplyr)
library(tidyr)
library(rootSolve)
library(conclique)
library(igraph)


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

## Graph conclique and single-site graph timings ----
## results saved from helpers/graph_timing.R
load("written_results/graph_timing.RData")

## Graph Timing plot, Figure 4 --------
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

