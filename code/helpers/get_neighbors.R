## Set of functions to find the neighbors in a network based on co-incidence
## written by: Andee Kaplan

# libraries -----
library(igraph)
library(dplyr)
library(purrrlyr)

# functions -----
find_neighbors <- function(row, edges) {
  edges[edges[,1] %in% row[,1:2] | edges[,2] %in% row[,1:2], "id"]
}

# n = the number of nodes
edge_neighbors <- function(n) {
  # reference
  g0 <- make_full_graph(n)
  e0 <- data.frame(as_edgelist(g0))
  e0$id <- seq_len(nrow(e0))
  
  e0 %>%
    by_row(find_neighbors, edges = e0, .collate = "cols") %>%
    select(-X1, -X2) %>%
    as.matrix() -> neighbors
  
  list(neighbors)
}






