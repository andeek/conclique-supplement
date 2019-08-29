## A script to create the concliques using DSatur and perform the analysis of Section 5.2
## Depends on data in the ../data/ folder
## written by: Daniel J. Nordman

### DSatur first###################
getNColors <- function(x) {
  if (is(x, "matrix")) {
    adj_mat <- x
    diag(adj_mat) <- FALSE
  } else if (is(x, "SpatialPolygons")) {
    adj_mat <- getAM(x)
  } else {
    stop("x must be an adjacency matrix or a SpatialPolygons* object.")
  }
  nColors <- length(unique(dsatur(adj_mat)))
  return(nColors)
}

## Get neighboring verteces
getNeighbors <- function(adj_mat, node_index) {
  nb <- which(adj_mat[node_index,])
  nb <- nb[!(nb==node_index)]
  return(nb)
}

## Count occurrences of color in given nodes
getAmountColor <- function(node_indexes, color_number, coloring) {
  node_colors <- coloring[node_indexes]
  return(sum(node_colors==color_number))
}

## Greedy DSATUR graph coloring algorithm
# Reference: D.Brelaz (1979) - New Methods to color the vertices of a graph. Communications of the ACM: 22(4).
# Ported from Python implementation by Andrei Novikov (pyclustering@yandex.ru)
# Under GNU Public license
dsatur <- function(x, coloring=NULL) {
  
  if (is.null(coloring)) {  # Set up vertex coloring from scratch
    color_counter = 1
    adj_mat <- x
    diag(adj_mat) <- FALSE
    degrees = list()
    saturation_degrees = rep(0, nrow(adj_mat))
    
    coloring = rep(0, nrow(adj_mat))
    uncolored_vertices = 1:nrow(adj_mat)
    
    index_maximum_degree = 0
    maximum_degree = 0
    for (index_node in 1:nrow(adj_mat)) {
      # Fill degree of nodes in the input graph
      degrees[[length(degrees)+1]] <- c(sum(adj_mat[index_node,]), index_node)
      
      # And find node with maximal degree at the same time.
      if ((degrees[[index_node]])[1] > maximum_degree) {
        maximum_degree <- (degrees[[index_node]])[1]
        index_maximum_degree <- index_node
      }
    }
    
    # Update saturation
    neighbors = getNeighbors(adj_mat, index_maximum_degree)
    for (index_neighbor in neighbors){
      saturation_degrees[index_neighbor] <- saturation_degrees[index_neighbor] + 1
    }
    
    # Coloring the first node
    coloring[index_maximum_degree] = color_counter
    uncolored_vertices <- uncolored_vertices[-index_maximum_degree]
    
  } else {  # Set up vertex coloring given input coloring
    color_counter = max(coloring)
    adj_mat <- x
    diag(adj_mat) <- FALSE
    degrees = list()
    saturation_degrees = rep(0, nrow(adj_mat))
    
    uncolored_vertices = 1:nrow(adj_mat)
    uncolored_vertices <- uncolored_vertices[coloring==0]
    
    # Fill degree of nodes in the input graph and update saturation
    for (index_node in 1:nrow(adj_mat)) {
      # Set degree
      degrees[[length(degrees)+1]] <- c(sum(adj_mat[index_node,]), index_node)
      # Set saturation
      index_neighbors <- getNeighbors(adj_mat, index_node)
      index_saturation <- 0
      for (number_color in 1:color_counter) {
        if (getAmountColor(index_neighbors, number_color, coloring) > 0) {
          index_saturation <- index_saturation + 1
        }
      }
      saturation_degrees[index_node] <- index_saturation
    }
  }
  
  
  # Color the remaining verteces
  while(length(uncolored_vertices) > 0) {
    # Get maximum saturation degree
    maximum_satur_degree = -1
    for (index in uncolored_vertices) {
      if (saturation_degrees[index] > maximum_satur_degree) {
        maximum_satur_degree = saturation_degrees[index]
      }
    }
    
    # Get list of indexes with maximum saturation degree
    indexes_maximum_satur_degree <- c()
    for (index in uncolored_vertices) {
      if (saturation_degrees[index] == maximum_satur_degree) {
        indexes_maximum_satur_degree <- c(indexes_maximum_satur_degree, index)
      }
    }
    
    coloring_index = indexes_maximum_satur_degree[1]
    if (length(indexes_maximum_satur_degree) > 1) {  # There are more then one node with maximum saturation
      # Find node with maximum degree
      maximum_degree = -1
      for (index in indexes_maximum_satur_degree) {
        degree <- (degrees[[index]])[1]
        node_index <- (degrees[[index]])[2]
        if (degree > maximum_degree) {
          coloring_index = node_index
          maximum_degree = degree
        }
      }
    }
    
    # Coloring
    node_index_neighbors = getNeighbors(adj_mat, coloring_index)
    for (number_color in 1:(color_counter)) {
      if (getAmountColor(node_index_neighbors, number_color, coloring) == 0) {
        coloring[coloring_index] = number_color
        break;
      }
    }
    
    # If it has not been colored then
    if (coloring[coloring_index] == 0) {
      color_counter <- color_counter + 1  # Add new color
      coloring[coloring_index] = color_counter
    }
    
    # Remove node from uncolored set
    uncolored_vertices <- uncolored_vertices[!(uncolored_vertices==coloring_index)]
    
    # Update degree of saturation
    for (index_neighbor in node_index_neighbors) {
      subneighbors = getNeighbors(adj_mat, index_neighbor)
      if (getAmountColor(subneighbors, coloring[coloring_index], coloring) == 1) {
        saturation_degrees[index_neighbor] <- saturation_degrees[index_neighbor] + 1
      }
    }
  }
  # Return final coloring
  return(coloring)
}

###################################   
####  Alliance Network

tmp<-read.table("../data/edgedat.txt",header=TRUE)
tmp1<-scan("../data/nbhds.txt")

### matrix J stores neighbor indices
J<-matrix(tmp1,ncol=554,byrow=T)


####  use Dsatur to find concliques per year 1946,1947,...,2007

##  get years for observations in J with neighbors
J.year<-tmp$year[J[,1]]

## input J, J.year, and then value for "yr"  (1946,1947,...,2007) 
## into 'dsatur.year'
### dsatur.year returns conclique number and observation index
## 
dsatur.year<-function(yr,J,J.year){
  
  t<-seq(1,nrow(J),1)
  t<-t[J.year==yr]
  
  J1<-J[t,]
  m<-nrow(J1)
  s<-J1[,1]
  s1<-seq(1,m,1)
  indm<-matrix(FALSE,nrow=m,ncol=m)
  
  for(i in 1:m){
    a<-J1[i,]
    a<-a[-1]
    a<-a[a>0]
    m1<-length(a)
    for(j in 1:m1){
      b<-a[j]
      d<-s1[s==b]
      indm[i,j]<-TRUE
      indm[j,i]<-TRUE
    }
  }
  
  dt<-dsatur(indm)
  list(concl=dt,obsindex=s)
  
}


###### collect concliques over years

concl<-as.null(0)
obsindex<-as.null(0)
for(i in 1946:2007){
  a<-dsatur.year(i,J,J.year)
  concl<-c(concl,a[[1]])
  obsindex<-c(obsindex,a[[2]])
}

###### 
###### note number of concliques is 2
concl

##### the number of observations in conclique 1 is 62
obsindex[concl==1]


###### any observation not appearing in obsindex has no neighbors 
######  & may be placed in conclique 2 


###### re-order the matrix J so that first 62 rows correspond to neighbors 
#####  of conclique 1


#### first,  J has 554 columns &
#### order the rows of J by # of neighbors

jt<-rep(0,nrow(J))
for(i in 1:nrow(J)){
  tt<-J[i,]
  tt<-tt[tt>0]
  jt[i]<-length(tt)
}
k<-order(jt,decreasing=TRUE)

for(i in 1:554){
  a<-J[,i]
  a<-a[k]
  J[,i]<-a}


##### now,
####  only the first 62 rows of J indicate more than 1 neighbor
####  each row in the first 62 rows corresponds to an observation in conclique 1

### indices match below
obsindex[concl==1]
J[1:62,1]


#### concordance "1-1" function (% of neighboring edges which are (1,1)) 
#### which is based on neighborhoods J....
#### For the 1st 62 rows of J, need to check which neighbors have an edge "1" 
####  if the first component of the row has edge "1" 
####   
concord<-function(data,J=J){
  s<-0
  for(i in 1:62){
    b<-J[i,]
    t<-b
    t<-t[t>0]
    t<-data[t]
    b<-t[1]
    t<-t[-1]
    if(b==1){
      s<-s+sum(b*t)}
  }
  
  # % of (1,1) edges occurring among neighbor pairs
  return(s/(nrow(J)-62))
}

### observed (1,1) % in network data
concord(tmp$defense,J=J) 

######################################
##########
######### Conclique-based Gibbs next
#### tmp is the original data table  above



### compute a design matrix XX just once

XX<-matrix(0,ncol=5,nrow=45513)
XX[,1]<-1
XX[,2]<-tmp$mil_ratio
XX[,3]<-log(tmp$tot_trade+1)
XX[,4]<-tmp$joint_dem

###  parameter estimates from Chyzh-Kaiser (2019) 

beta<-c(0.094, -2.363, 0.015,0.884,0.016) 

### using parameter estimates & design matrix,
### compute the logit form of conditional probabilities
### up to the one term that depends on sums of neighboring observations
### all else in the logit form follows from design matrix XX & parameters beta  

r<- t(XX)*beta
r<- t(r)
r<- apply(r,1,sum)
r1<-r-beta[5]*XX[,5]
reb1<-r1
re<- exp(r1)
reb<- re/(re+1)


for(j in 1:62){
  b<-J[j,1]
  t<-J[j,]
  t<-t[t>0]
  t<-t[-1]
  reb1[b]<-reb1[b] -beta[5]*sum(reb[t])
}

b1<-J[63:nrow(J),1]
b2<-J[63:nrow(J),2]
reb1[b1]<- reb1[b1]- beta[5]*reb[b2]



### now just have to add "reb1" to "beta[5]*'sum of neighbors'"
##  to obtain the logit conditional prob of '1' per observation



#### conclique-based sampler 
###  performs one full Gibbs iteration/returns network sample
#### incon are intial values, reb1 are root form of logit conditional probabilities of 1
##### where one needs to add beta[5]*"sum of neighbors" to reb to get the logit conditional probability
#####  J is the (ordered) matrix above to indicate neighbors, beta is the parameter vector of length 5   
#
consim<-function(incon,reb1=reb1,J=J,beta=beta){
  h1<-incon
  h<-rep(0,62)
  ## first 62 obs in J/conclique 1 updates from conclique 2
  for(j in 1:62){
    b<-J[j,1]
    t<-J[j,]
    t<-t[t>0]
    t<-t[-1]
    h[j]<-reb1[b]+sum(h1[t])*beta[5]
  }
  ## cond prob of 1
  h<-exp(h)/(exp(h)+1)
  ## compare to uniform to get binary values
  u<-runif(62)
  b<-J[1:62,1]
  ht<-h*0
  ht[u<h]<-1
  h1[b]<-ht
  
  
  ## remaining conclique 2 updates from updated conclique 1
  b1<-J[63:nrow(J),1]
  b2<-J[63:nrow(J),2]
  h<-reb1[b1]+beta[5]*h1[b2]
  ## cond prob of 1
  h<-exp(h)/(exp(h)+1)
  ## compare to uniform to get binary values
  u<-runif(length(h))
  ht<-h*0
  ht[u<h]<-1
  h1[b1]<-ht
  return(h1)
}


## intialize sampler
mm<-consim(rep(0,45513),reb1=reb1,J=J,beta=beta)

## burn-in 10000 full iterations
for(i in 1:10000){
  mm<-consim(mm,reb1=reb1,J=J,beta=beta)}


### sample 2500 network statistics with thinning by 100 
### compute statistic/concordance % on each sample
res<-rep(0,2500)

for(i in 1:2500){
  for(j in 1:100){
    mm<-consim(mm,reb1=reb1,J=J,beta=beta)}
  res[i]<-concord(mm,J=J) 
  print(i)
}

### MCMC p-value
# observed statistic/concordance %
a<-concord(tmp$defense,J=J)
## p-value
length(res[res>a])/length(res)
