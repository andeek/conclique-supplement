## Set of functions to find the concliques in a network based on co-incidence
## written by: Andee Kaplan

# libraries -----
library(igraph)
library(dplyr)

# functions -----
edge_concliques_even <- function(n) {
  # label nodes on the circle
  nodes <- data.frame(actual = 1:n)
  v0 <- sample(n, size = 1)
  nodes$label <- NA
  nodes[nodes$actual != v0, "label"] <- 0:(n-2)
  
  Q <- n - 1
  conc <- list(Q)
  for(j in seq_len(Q)) {
    inc <- data.frame(tip = (j - 1 + 1:((n - 2)/2)) %% (n - 1), tail = (j - 1 - 1:((n - 2)/2)) %% (n - 1))
    inc <- rbind(inc, c(NA, j - 1))
    
    inc %>% left_join(nodes, by = c("tip" = "label")) %>% 
      rename(X1 = actual) %>%
      left_join(nodes, by = c("tail" = "label")) %>%
      rename(X2 = actual) %>%
      select(starts_with("X")) -> conc[[j]]
  }
  return(conc)
}
edge_lookup <- function(edge_list, n) {
  # reference
  g0 <- make_full_graph(n)
  e0 <- data.frame(as_edgelist(g0))
  e0$id <- seq_len(nrow(e0))
  
  # list of numeric vectors (edge position)
  conc <- lapply(edge_list, function(co) {
    e0 %>% 
      right_join(co, by = c("X1", "X2")) %>%
      filter(!is.na(id)) %>%
      .$id -> e
    
    e0 %>%
      right_join(co, by = c("X1" = "X2", "X2" = "X1")) %>%
      filter(!is.na(id)) %>%
      .$id %>%
      c(e)
  })
  
  return(conc)
}

# edge concliques -----
# n = the number of nodes
edge_concliques <- function(n) {
  if(n %% 2 == 0) {
    # edge concliques
    edge_list <- edge_concliques_even(n)
  } else {
    # edge concliques w/ fake node
    edge_list <- edge_concliques_even(n + 1)
    
    # remove fake node
    edge_list <- lapply(edge_list, function(conc) {
      conc[rowSums(conc == (n + 1)) == 0,]
    })
  }
  edge_lookup(edge_list, n)
}




