library(igraph)
library(rstudioapi)
library(ggplot2)
library(dplyr)
library(ergm)
library(intergraph)
library(coda)
library(sbm)


getwd()
setwd(dirname(getActiveDocumentContext()$path))
set.seed(1)
g=read.graph("mouse_visual.cortex_2.graphml", format = "graphml")
g

table(V(g)$type1)
table(V(g)$type2)

length(E(g))
n=length(V(g))
n
#conteggio archi, conteggio nodi di un tipo ecc..

# undirected
un.g = as.undirected(g, mode = "collapse")


###GRAFICI 

#SENZA NODI COLORATI
#TIPO 1 
colors1 <- ifelse(V(g)$type1=="Dendritic fragment", "deepskyblue", "gold")
colors1[V(g)$type1=="Cell body in EM volume"] <- "limegreen"

# plotta il grafo con i nodi colorati
plot(g, vertex.size = 5, vertex.color = colors1, vertex.label.cex=0.7, vertex.shape = "circle", edge.arrow.size = 0.3, 
     edge.width = 1, main="Type 1")
colors1=c("deepskyblue", "gold", "limegreen")
legend("topleft", legend= c("Dendritic fragment", "Characterized pyramidal neuron", "Cell body in EM volume"), col = colors1, pch = 16, cex = 0.7)

get.edge.attribute(g)

#TIPO 2
# specifica i colori per i diversi tipi di nodi
colors <- ifelse(V(g)$type2=="NA", "gold", "deepskyblue")
colors[V(g)$type2=="Postsynaptic excitatory target"] <- "limegreen"

# plotta il grafo con i nodi colorati
plot(g, vertex.size = 5, vertex.color = colors, vertex.label.cex=0.8, vertex.shape = "circle", edge.arrow.size = 0.4, 
     edge.width = 1, main="Type 2")
colors1=c("gold", "deepskyblue", "limegreen")
legend("topleft", legend= c("NA", "Excitatory", "Inhibitory"), col = colors1, pch = 16, cex = 0.8)

#CREDO CHE I NODI NA MANDINO SOLAMENTE, MA QUANDO RICEVONO LO FANNO SOLAMENTE DA UN ALTRO NA



#### NETWORK STATISTICS ####


# Q1: How connected is the network?
# ----------------------------------- 
# DENSITY
graph.density(g)
#la probabilità di osservare una relazione tra due nodi del network scelti casualmente è 0.005656886, 0.5%.
#Abbastanza sparso.


# Q2: is there a tendency for the nodes in a directed network to return relations?
# ---------------------------------------------------------------------------------
# RECIPROCITY
# DYADS
reciprocity(g)
#non c'è reciprocità nel network - cliques di dimensione 2, star-2
#un nodo (neurone) o è eccitatorio o inibitorio o na (?)

# Q3: is there a tendency of the nodes in the network to cluster together?
# --------------------------------------------------------------------------
# TRANSITIVITY 
# are friends of friends friends too?
# TRIADS

# Let us create a vector of labels with all possible relations
census.labels = c('empty',
                  'A->B, C',
                  'A<->B, C',
                  'A<-B->C',
                  'A->B<-C',
                  'A->B->C',
                  'A<->B<-C',
                  'A<->B->C',
                  'A->B<-C, A->C',
                  'A<-B<-C, A->C',
                  'A<->B<->C',
                  'A<-B->C, A<->C',
                  'A->B<-C, A<->C',
                  'A->B->C, A<->C',
                  'A->B<->C, A<->C',
                  'A<->B<->C, A<->C')

g.tri = triad.census(g)
data.frame(census.labels, g.tri)

#TRANSITIVITY COEFFICIENT

#transitivity
tr.g = transitivity(g)
tr.g
#0.0046875
#è la probabilità condizionata di osservare una relazione tra due nodi che hanno un amico in comune 
#0.4% 

#normalizing the transitivity index  
#let us compare it with the density
dens.g = graph.density(g)
dens.g

# compute the log-odds ratio
odd.dens.g = dens.g/(1-dens.g)
odd.tr.g = tr.g/(1-tr.g)

#odds ratio
odd.tr.g/odd.dens.g
#0.827829 
#La probabilità di osservare un legame tra nodi che condividono un amico in comune è circa il 20% più bassa di quella
#prevista dal caso.


# Q4: are highly connected nodes similar to each other?
# ------------------------------------------------------
# ASSORTATIVE MIXING
#TIPO 1
assortativity(g, as.numeric(V(g)$type1))
cell.nodes=V(g)[which(V(g)$type1=="Cell body in EM volume")]
pyra.nodes=V(g)[which(V(g)$type1=="Characterized pyramidal neuron")]
dend.nodes=V(g)[which(V(g)$type1=="Dendritic fragment")]

V(g)[cell.nodes]$type1 = 0
V(g)[pyra.nodes]$type1 = 1
V(g)[dend.nodes]$type1 = 2
V(g)$type1

g_assortativity1 = assortativity(g, as.numeric(V(g)$type1))
g_assortativity1
#0.009884543 assortative mixing è più probabile osservare una relazione tra nodi che hanno attributi simili
#valore molto vicino allo zero, type1 influisce poco ???
#ha senso farlo col grafo indiretto ???
#si può confrontare con la densità ???


#TIPO 2
assortativity(g, V(g)$type2)
NA.nodes=V(g)[which(V(g)$type2=="NA")]
excit.nodes=V(g)[which(V(g)$type2=="Postsynaptic excitatory target")]
inhib.nodes=V(g)[which(V(g)$type2=="Postsynaptic inhibitory target")]

V(g)[NA.nodes]$type2 = 0
V(g)[excit.nodes]$type2 = 1
V(g)[inhib.nodes]$type2 = 2
V(g)$type2


g_assortativity2 = assortativity(g, as.numeric(V(g)$type2))
g_assortativity2
#Nan, questo perchè per calcolare l'assortatività al denominatore abbiamo il numero di 
#connessioni tra il nodo i e j ma noi tra nodi della stessa categoria type2 non abbiamo connessioni
#CHIEDERE


#### NODAL STATISTICS ####

# 1) DEGREE CENTRALITY
# ---------------------
# A node is central if it is connected to many other nodes
# zeta_i^d = sum_j y_ij
# degree centrality
degree(g)
degree(g, normalized = T) #degree(g)/(n-1)
# how many nodes do point to a given node? 
in.degree = (degree(g, mode="in"))

in.degree_sort = sort(degree(g, mode="in"), decreasing=TRUE)
table(in.degree_sort)
#10 nodi ricevono 0; 162 nodi ricevono 1; 19 nodi ricevono 2; 2 nodi ricevono 3; 2 nodi ricevono 4 

V(g)[which.max(degree(g, mode="in"))]
#il nodo più centrale è il num. 83
V(g)[which.min(degree(g, mode="in"))]



degree(g, v= "83", mode= "in")
#1 nodo che riceve 4

# how many nodes are pointed by a given node?
out.degree=degree(g, mode = "out")

out.degree_sort=sort(degree(g, mode = "out"), decreasing=TRUE)
table(out.degree_sort)

V(g)[which.max(degree(g, mode="out"))]

#Interpretazione dei nodi più popolari ? 

#ISTOGRAMMA
df.degree=data.frame(cbind(type=c(rep("In degree", n), rep("Out degree", n)), value=c(degree(g, mode="in"), degree(g, mode = "out"))))

df.degree %>%
  ggplot( aes(x=as.numeric(value), fill=type)) +
  geom_histogram( color="#e9ecef", alpha=0.6, stat = "count", binwidth = 30, position="dodge") +
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  scale_x_continuous(breaks = 0:35) +
  labs(fill="") +
  xlab("Degree") +
  ylab("Frequency")
#we can see a two main trend: a lot of nodes sent 0 edges and a lot of nodes receive just 1 edge
#so we have few nodes that sent a lot of edges and a lot of nodes that just receive edges from these few nodes
#10 of these 14 nodes that sent a lot of edges don't receive any edges 

#PLOT BY DEGREE
#colors 
#in degree
in.ord = order(in.degree, decreasing = T)
V(g)$color = "orange"
V(g)$color[in.ord[1:3]] = "lightblue"
plot(g, vertex.size = in.degree*5,  edge.arrow.size = 0.5, vertex.label.cex=0.8, vertex.shape = "circle", edge.arrow.size = 0.4, 
     edge.width = 1)

#out degree
out.ord = order(out.degree, decreasing = TRUE)
V(g)$color = "orange"
V(g)$color[out.ord[1:3]] = "lightblue"
plot(g, vertex.size = out.degree,  edge.arrow.size = 0.5, vertex.label.cex=0.8, vertex.shape = "circle", edge.arrow.size = 0.4, 
     edge.width = 1)

# 2) CLOSENESS CENTRALITY
# ------------------------
# A node is central if it is "close" to many other nodes
# z_i^c = 1/sum_j d_ij
# closeness centrality
allin.close = centr_clo(g, mode = "all")
allin.close$centralization
# how easily can a node be reached from other nodes? 
in.close=closeness(g, mode = "in")
centr_clo(g, mode = "in")$centralization
# how easily can a node reach the others? 
out.close=closeness(g, mode = "out")
centr_clo(g, mode = "out")$centralization
#Un nodo è centrale se è vicino ad altri nodi, geodesic distance - diametro. 
#Quanto veloce riesco a raggiungere da un nodo tutti gli altri 

in.close_sort = sort(closeness(g, mode = "in"))
in.close_sort 

out.close_sort = sort(closeness(g, mode = "out"))
out.close_sort


#in closeness
in.close[which(is.nan(in.close))]=0
in.ord = order(in.close, decreasing = T)
V(g)$color = "orange"
V(g)$color[in.ord[1:3]] = "lightblue"
plot(g, vertex.size = in.close*10,  edge.arrow.size = 0.5, vertex.label.cex=0.8, vertex.shape = "circle", edge.arrow.size = 0.4, 
     edge.width = 1, main="In-closeness")
#out closeness
out.close[which(is.nan(out.close))]=0
plot(g, vertex.size = out.close*100,  edge.arrow.size = 0.5, vertex.label.cex=0.8, vertex.shape = "circle", edge.arrow.size = 0.4, 
     edge.width = 1, main="Out-closeness")
#undirected closeness
un.close=closeness(un.g)
un.close[which(is.nan(un.close))]=0
plot(g, vertex.size = un.close, edge.arrow.size = 0.5, vertex.label.cex=0.8, vertex.shape = "circle", edge.arrow.size = 0.4, 
     edge.width = 1)


# 3) BETWEENNESS CENTRALITY
# --------------------------
# a node is central if it is located between many nodes
g.bet=betweenness(g, directed = F)
plot(g, vertex.size = g.bet*0.005, edge.arrow.size = 0.5, vertex.label.cex=0.8, vertex.shape = "circle", edge.arrow.size = 0.4, 
     edge.width = 1, main="Betweenness centrality")
# 4) EIGENVECTOR CENTRALITY
# --------------------------
# a node is central if it is connected to other central nodes
g.eigen=eigen_centrality(g, scale = F)$vector
plot(g, vertex.size = g.eigen*30,  edge.arrow.size = 0.5,
     vertex.label.cex = 1.5)




################ MODELS #################

#### SIMPLE RANDOM GRAPH ####
#### MLE approach ####


# 2. Is the observed network coherent with the family of Binomial random graph models G(n, p)?
# a naive approach based on the best case scenario
# ------------------------------------------------------------------------------------
#assunzione di indipendenza tra i nodi, la probabilità di una relazione è la stessa per ogni coppia di nodi 
Y = get.adjacency(g, sparse = F)
diag(Y) = NA
# number of nodes
n = nrow(Y)

rho.obs = mean(Y, na.rm = T)
rho.obs
#0.005656886
#densità osservata nei dati
#la probabilità di osservare una relazione tra due nodi scelti a caso è del 0.6%
CId.obs = centr_degree(g, mode = "in", loops = F)$centralization # centralization
CId.obs  
#0.01503879 - centralization 
C.obs = transitivity(g) # clustering coefficient
C.obs
#0.0046875

# maximum likelihood estimate of p
p.MLE = mean(Y, na.rm = T)

B = 1000
rho.sim = CId.sim = C.sim = c()
for(b in 1:B){
  tmp = rbinom(n^2,1,p.MLE)  
  Y.sim = matrix(tmp, n,n)
  diag(Y.sim) = NA
  rho.sim[b] = mean(Y.sim, na.rm = TRUE)
  g.sim = graph_from_adjacency_matrix(Y.sim)
  CId.sim[b] = centr_degree(g.sim, mode = "in", loops = F)$centralization
  C.sim[b] = transitivity(g.sim)
}


# Graphical comparison
# density
# of course it's in the middle, we use density for simulation, più è vicina al centro meno evidenza ho contro h0, ovvero che
# using density to understand choerence is not so good
par(mfrow = c(1,3))
low = pmin(min(rho.sim), rho.obs) 
up = pmax(max(rho.sim), rho.obs) 
hist(rho.sim, col = "lightgray", main = "Null distribution", xlim = c(low, up), xlab="Density")
abline(v = rho.obs, col = "red", lwd=2) 

# centralization
low = pmin(min(CId.sim), CId.obs)
up = pmax(max(CId.sim), CId.obs) 
hist(CId.sim, col = "lightgray", main = "Null distribution", xlim = c(low, up), xlab="Centralization")
abline(v = CId.obs, col = "red", lwd=2)

# clustering coefficient
low = pmin(min(C.sim), C.obs) 
up = pmax(max(C.sim), C.obs) 
hist(C.sim, col = "lightgray", main = "Null distribution", xlim = c(low, up), xlab="Clustering coefficient")
abline(v = C.obs, col = "red", lwd=2)


# compute an approximate p-value
mean(rho.sim >= rho.obs)
#0.51
#al livello del 5% non rigettiamo h0
mean(CId.sim < CId.obs) 
#0.138 
#al livello del 5% non rigettiamo h0
mean(C.sim < C.obs)
#0.201 non rigettiamo h0



#### Conditional uniform distribution ####
#creiamo la matrice con il numero di connessioni fissato uguale a quello del network osservato, quindi la degree del grafo è sempre uguale a 0.005656886, cioè alla densità osservata
B = 1000
m =  sum(Y, na.rm = TRUE)
rho.sim = CId.sim = C.sim = c()
for(b in 1:B){ #devo costruire la matrice di adiacenza composta da il numero di 1 ripetuto m volte, cioè le connessioni, mentre il resto della matrice sarà formato dal numero 0 ripetuto per tutte le possibili connessioni meno quelle effettive, quindi n*(n-1) - m
  Y.sim = matrix(, n, n) #creo la matrice vuota
  ones = rep(1, m) #creo le connessioni
  zeros = rep(0, n*(n-1) - m) #creo le non connessioni
  all = c(ones, zeros)
  Y.sim[col(Y.sim) != row(Y.sim)] = sample(all, n*(n-1))
  g.sim = graph_from_adjacency_matrix(Y.sim)
  rho.sim[b] = mean(Y.sim, na.rm = TRUE)
  CId.sim[b] = centr_degree(g.sim, mode = "in", loops = F)$centralization
  C.sim[b] = transitivity(g.sim)
}


# Graphical comparison
par(mfrow = c(1,3))


hist(rho.sim, col = "lightgray", main = "Null distribution", xlab = "Density")
abline(v = rho.obs, col = "red", lwd=2)

hist(CId.sim, col = "lightgray", main = "Null distribution", xlab = "Centralization")
abline(v = CId.obs, col = "red", lwd=2)

hist(C.sim, col = "lightgray", main = "Null distribution", xlab = "Clustering coefficient")
abline(v = C.obs, col = "red", lwd=2)


# compute an approximate p-value
mean(CId.sim < CId.obs) #rigettiamo h0
mean(C.sim < C.obs) #non rigettiamo

#results are pretty similar with the one computed with MLE approach


#### IN-DEGREE, OUT-DEGREE AND RECIPROCITY ####

din = degree(g, mode = "in")
dout = degree(g, mode = "out")

#### MLE approach ####
p = graph.density(g) #MLE
sdIn.sim = sdOut.sim = recip.sim =  c() #I look at the sd to summarize statistics features
for(b in 1:1000){
  tmp = rbinom(n^2,1,p)  
  Y.sim = matrix(tmp, n,n); diag(Y.sim) = NA
  g.sim = graph_from_adjacency_matrix(Y.sim)
  sdIn.sim[b] = sd(degree(g.sim, mode = "in"))
  sdOut.sim[b] = sd(degree(g.sim, mode = "out"))
  recip.sim[b] = reciprocity(g.sim)
}
# sdIn, sdOut and recip are kind of realization of network statistics under the Null hypothesis
# graphical comparison
par(mfrow = c(1,3))
hist(sdIn.sim, breaks = 20, col = "lightgray", main = "", 
     prob = TRUE, xlab = "sd(in-degree)", xlim = c(0.5 ,1.5))
abline(v = sd(din), col = "red", lwd = 2)


hist(sdOut.sim, breaks = 20, col = "lightgray", main = "", 
     prob = TRUE, xlab = "sd(out-degree)", xlim = c(0.5,1.5))
abline(v = sd(dout), col = "red", lwd = 2)

hist(recip.sim, breaks = 20, col = "lightgray", main = "", 
     prob = TRUE, xlab = "reciprocity")
abline(v = reciprocity(g), col = "red", lwd = 2)

#is appropriate only in reciprocity terms
mean(sdIn.sim < sd(din)) #reject
mean(sdOut.sim > sd(dout)) #reject
mean(recip.sim > reciprocity(g)) #cannot reject


#### CONDITIONAL UNIFORM DISTRIBUION MODEL ####

m = ecount(g)
B = 1000
sdIn.sim = sdOut.sim = recip.sim = c()
for(b in 1:B){
  ones = rep(1, m)
  zeros = rep(0, n*(n-1) - m)
  all = c(ones, zeros)
  tmp = sample(all, n*(n-1)) #sampling from conditional uniform distribution, it's like let tmp be random
  Y.sim = matrix(tmp, n,n)
  diag(Y.sim) = NA
  g.sim = graph_from_adjacency_matrix(Y.sim)
  sdIn.sim[b] = sd(degree(g.sim, mode = "in"))
  sdOut.sim[b] = sd(degree(g.sim, mode = "out"))
  recip.sim[b] = reciprocity(g.sim)
}

# graphical comparison
par(mfrow = c(1,3))
hist(sdIn.sim, breaks = 20, col = "lightgray", main = "", 
     prob = TRUE, xlab = "sd(in-degree)", xlim = c(0.5,1.5))
abline(v = sd(din), col = "red", lwd = 2)


hist(sdOut.sim, breaks = 20, col = "lightgray", main = "", 
     prob = TRUE, xlab = "sd(out-degree)", xlim = c(0.5,1.5))
abline(v = sd(dout), col = "red", lwd = 2)


hist(recip.sim, breaks = 20, col = "lightgray", main = "", 
     prob = TRUE, xlab = "reciprocity",xlim = c(0,0.1) )
abline(v = reciprocity(g), col = "red", lwd = 2)

#AIC e BIC 

#### NON HOMOGENEOUS SIMPLE RANDOM GRAPH MODEL ####

# response variable 
# remember how R memorize the elements of the matrix
Y = as.matrix(get.adjacency(g), sparse = F)
diag(Y) = NA
y = c(Y)

# covariates -- sender and receiver effects
rowIdx = row(Y) #this produces a matrix having 1 in the first row, 2 in the second row, ..., n in the n-th row 
colIdx = col(Y) #this produces a matrix having 1 in the first column, 2 in the second column, ..., n in the n-th column
rowIdx[1:4, 1:4]
colIdx[1:4, 1:4]

rowidx = c(rowIdx)
colidx = c(colIdx)

# estimate the parameters
# these are essentially mu, a_i and b_j all grouped together
mod = glm(y ~ factor(rowidx) + factor(colidx), family = "binomial")
summary(mod)
#non c'è nessun coefficiente significativo, non va bene 
pij = mod$fitted.values

# null distribution
sdIn.sim = sdOut.sim = recip.sim = c()
rho.sim = CId.sim = C.sim = c()

for(b in 1:1000){
  tmp = rbinom(n*(n-1),1,pij) # indipendent Bernoulli random variables, so reciprocity it's kinda zero
  Y.sim = matrix(, n,n); Y.sim[row(Y.sim) != col(Y.sim)] = tmp
  g.sim = graph_from_adjacency_matrix(Y.sim)
  
  rho.sim[b] = mean(Y.sim, na.rm = TRUE)
  CId.sim[b] = centr_degree(g.sim, mode = "in", loops = F)$centralization
  C.sim[b] = transitivity(g.sim)
  
  sdIn.sim[b] = sd(degree(g.sim, mode = "in"))
  sdOut.sim[b] = sd(degree(g.sim, mode = "out"))
  recip.sim[b] = reciprocity(g.sim)
}

# graphical comparison

par(mfrow = c(1,3))


hist(rho.sim, col = "lightgray", main = "Null distribution", xlab = "Density")
abline(v = rho.obs, col = "red", lwd=2)

hist(CId.sim, col = "lightgray", main = "Null distribution", xlab = "Centralization")
abline(v = CId.obs, col = "red", lwd=2)

hist(C.sim, col = "lightgray", main = "Null distribution", xlab = "Clustering coefficient")
abline(v = C.obs, col = "red", lwd=2)





par(mfrow = c(1,3))
hist(sdIn.sim, breaks = 20, col = "lightgray", main = "", 
     prob = TRUE, xlab = "sd(in-degree)")
abline(v = sd(din), col = "red", lwd = 2)

hist(sdOut.sim, breaks = 20, col = "lightgray", main = "", 
     prob = TRUE, xlab = "sd(out-degree)")
abline(v = sd(dout), col = "red", lwd = 2)

hist(recip.sim, breaks = 20, col = "lightgray", main = "", 
     prob = TRUE, xlab = "reciprocity")
abline(v = reciprocity(g), col = "red", lwd = 2)

#ora riusciamo a cogliere anche la OUT-DEGREE
mean(sdIn.sim < sd(din))#reject
mean(sdOut.sim < sd(dout)) # cannot reject
mean(recip.sim < reciprocity(g)) #cannot reject



#### WE WANT THE WHOLE FAMILY REPRESENTING THE OBSERVED NETWORK, NOT JUST A MEMBER ####
# simulate from the conditional uniform distribution
# simulating function -- alternating rectangles
aRect.fnc = function(Y, k){
  
  # Y = adjacency matrix
  # k = n. of steps in the alternating rectangles algorithm 
  # I'm drawing four random nodes to extract rows and columns
  # to see if we are connected or not keeping the number of links always constant
  
  #this algo produces networks really similar to the one observed
  #it's like the naive approach that uses the MLE
  
  Y1 = matrix(c(0,1,1,0), 2, 2)
  Y2 = 1 - Y1
  
  n = nrow(Y)
  
  for(s in 1:k){
    # draw 4 distinct indexes
    # two rows and two columns
    ij = sample(1:n,4,replace = F)
    
    # select the corresponding sub-matrix
    rows = ij[1:2]
    cols = ij[3:4]
    Yij = Y[rows, cols]
    
    # perturbation
    if(all(Yij == Y1)) Yij = Y2 else if(all(Yij == Y2))  Yij = Y1 #if(TRUE che Yij == Y1) allora Yij = Y2 else if(TRUE che Yij == Y2) allora Yij = Y
    # se la matrice selezionata non ? uguale n? a Y1 n? a Y2 allora rimarr? uguale al suo valore iniziale
    # altrimenti nella sua posizione inseriamo la nuova matrice, cio? Y1 o Y2
    Y[rows, cols] = Yij
  }
  
  return(Y)
}


sdIn.sim = sdOut.sim = recip.sim = c()
rho.sim = CId.sim = C.sim = c()

for(b in 1:1000){
  Y.sim = aRect.fnc(Y, 100)
  # print number of perturbed elements in Y.sim 
  cat(sum(Y != Y.sim, na.rm = T), "*", sep="") # number of elements that are switched
  g.sim = graph_from_adjacency_matrix(Y.sim)
  
  
  rho.sim[b] = mean(Y.sim, na.rm = TRUE)
  CId.sim[b] = centr_degree(g.sim, mode = "in", loops = F)$centralization
  C.sim[b] = transitivity(g.sim)
  
  
  sdIn.sim[b] = sd(degree(g.sim, mode = "in"))
  sdOut.sim[b] = sd(degree(g.sim, mode = "out"))
  recip.sim[b] = reciprocity(g.sim)
  
}



# graphical comparison
par(mfrow = c(1,3))


hist(rho.sim, col = "lightgray", main = "Null distribution", xlab = "Density")
abline(v = rho.obs, col = "red", lwd=2)

hist(CId.sim, col = "lightgray", main = "Null distribution", xlab = "Centralization")
abline(v = CId.obs, col = "red", lwd=2)

hist(C.sim, col = "lightgray", main = "Null distribution", xlab = "Clustering coefficient")
abline(v = C.obs, col = "red", lwd=2)



par(mfrow = c(1,3))
hist(sdIn.sim) # le marginali non cambiano -> sdIn.sim ? uguale a sd(din)
abline(v = sd(din), col = "red", lty = 2)
sdIn.sim
#0.5234142 degree.in 

hist(sdOut.sim) # le marginali non cambiano -> sdOut.sim ? uguale a sd(dout)
abline(v = sd(dout), col = "red", lty = 2)
sdOut.sim
#5.05493

hist(recip.sim)
abline(v = reciprocity(g), col = "red", lty = 2)

# p-value
mean(sdIn.sim > sd(din))
mean(sdOut.sim > sd(dout))
mean(recip.sim > reciprocity(g))
#reject all


# let us also look at the transitivity of the network
cl.sim = c()
for(b in 1:1000){
  Y.sim = aRect.fnc(Y, 100)
  # print number of perturbed elements in Y.sim 
  cat(sum(Y != Y.sim, na.rm = T), "*", sep="")
  g.sim = graph_from_adjacency_matrix(Y.sim)
  cl.sim[b] = transitivity(g.sim)
}


# graphical comparison
par(mfrow = c(1,1))
hist(cl.sim)
abline(v = transitivity(g), col = "red", lty = 2)

# p-value
mean(cl.sim < transitivity(g))
#reject


#### EXPONENTIAL RANDOM GRAPH MODEL ####
net = network(Y, directed = T)

#aggiungiamo gli attributi degli archi 
attr=vertex.attributes(g)
attr$type2[which(attr$type2=="NA")] <- "Other"
net %v% "type1" = attr$type1
net %v% "type2" = attr$type2
net
net %e% "Weight" = edge_attr(g)$weight
net

#************************************
# NULL model: BRG o SRG
#************************************
mod0 = ergm(net ~ edges)
mod0
#edges -5.16921, significativo 
#lower number of ties that we expected by chance
# look in more depth
summary(mod0)
# indeed...
# the odds of observing a relation between two randomly 
# selected nodes is about 99.5% lower than that of not observing it
mcmc.diagnostics(mod0)

# -------------------------------------------
# non-homogeneous BRG
# -------------------------------------------
# which sufficient statistics for theta?
# n. of edges
# in-degrees -> receiver effects
# out-degrees -> sender effects
mod1 = ergm(net ~ edges + sender + receiver)
summary(mod1)
mcmc.diagnostics(mod1)

# which model is more appropriate? 
BIC(mod0, mod1)
AIC(mod0, mod1)
#visto che molti coeff. del modello 1 sono non significativi seguiamo le indicazioni del BIC e prendiamo mod0
#anche se BIC è più conservativo 


#stiamo inserendo il reciprocity parameter - bidirezionali. noi non ne abbiamo! non ha senso.
#inoltre abbiamo tanti nodi che non mandandano nessun arco o viceversa non ricevono
#la specificazione con sender, receiver e mutual non è buona  

mod2 = ergm(net ~ edges + mutual)
summary(mod2) #torna! perchè non abbiamo reciprocità, infatti il coeff di mutual è infinito

mcmc.diagnostics(mod2)

pdf("C:/Users/David/Desktop/Unifi/SECONDO ANNO/Secondo semestre/Statistical analysis of network data/Progetto/mod2_diagn.pdf")
mcmc.diagnostics(mod2)
dev.off()

#let include the nodal attributes 
#nodefactor(attr) : se i nodi che hanno lo stesso attributo partecipano di più nel network 
#nodematch(attr) : se nodi con lo stesso attributo sono più connessi tra loro 

#type1 e type2
mod3 = ergm(net ~ edges + nodefactor("type2") + nodefactor("type1") + 
              nodematch("type2") + nodematch("type1"), control = control.ergm(seed = 1)) 
summary(mod3)
#l'unico coefficiente significativo è quello di nodematch(type1)
mcmc.diagnostics(mod3) 


#type1
mod4 = ergm(net ~ edges + nodefactor("type1")+ 
              nodematch("type1"), control = control.ergm(seed = 1)) 
summary(mod4)
#convergence
mcmc.diagnostics(mod4)

#type2
mod5 = ergm(net ~ edges + nodefactor("type2")+ 
              nodematch("type2"), control = control.ergm(seed = 1)) 
summary(mod5) #nessun coefficiente è significativo 
mcmc.diagnostics(mod5)
#only main effects
mod6 = ergm(net ~ edges + nodefactor("type1")+ nodefactor("type2"), control = control.ergm(seed = 1)) 
summary(mod6)
mcmc.diagnostics(mod6)
#only homophily
mod7 = ergm(net ~ edges + nodematch("type1")+ nodematch("type2"), control = control.ergm(seed = 1)) 
summary(mod7)

#mix
mod8 = ergm(net ~ edges + nodefactor("type1") + nodematch("type1")+ nodematch("type2"), control = control.ergm(seed = 1)) 
summary(mod8)

# cannot use BIC and AIC for model selection because of some infinite sender coefficients
BIC(mod0, mod1, mod3, mod4, mod5, mod6, mod7, mod8)
AIC(mod0, mod1, mod3, mod4, mod5, mod6, mod7, mod8)


#### MARKOV GRAPH MODEL ####
# Let us add to the model the triangle term, the in- and the out-stars of order 2
# (indicating the tendency to form clusters in the network)

mod9 = ergm(net ~ edges + mutual + nodefactor("type1") + istar(2) + ostar(2) + triangle, 
            control = control.ergm(seed = 1))
summary(mod9)

# let us try to remove triangles
mod10 = ergm(net ~ edges + mutual + nodefactor("type1") + istar(2) + ostar(2), 
             control = control.ergm(seed = 1))

#non girano! neanche alla prof, normale 

# let us try to solve the issue by considering the alternating k-star term
# a standard choice for decay is 1, but model selection can be used!
#in questo modello abbiamo come statistiche sufficienti:
#main effects: nodefactor()
#triangles
#gwidegree: number of in stars 
#gwodegree: number of out stars 

 #in +  main effects 
mod11 = ergm(net ~ edges + triangle + nodefactor("type1") + nodefactor("type2") + gwidegree(decay = 1, fixed = TRUE), 
             control = control.ergm(seed = 1))
summary(mod11)
par(mfrow=c(3,1))
par(ask = TRUE)
mcmc.diagnostics(mod11)

pdf("C:/Users/David/Desktop/Unifi/SECONDO ANNO/Secondo semestre/Statistical analysis of network data/Progetto/mod11_diagn.pdf")
mcmc.diagnostics(mod11)
dev.off()





#out + main effects 
mod12 = ergm(net ~ edges + nodefactor("type1")+ nodefactor("type2") + gwodegree(decay = 1, fixed = TRUE), 
             control = control.ergm(seed = 1))
summary(mod12) #non converge 

# ----------------------------------------
# let us consider the social circuit model 
# ----------------------------------------
# alternating k-triangles --> gwesp(decay = 0, fixed = FALSE) 
# geometrically weighted edge-wise shared partners
# the corresponding parameter expresses the tendency for tied nodes 
# to have multiple shared partners

# alternating k-2-paths --> gwdsp(decay = 0, fixed = FALSE)
# geometrically weighted dyad-wise shared partners
# the corresponding parameter expresses the tendency for dyads 
# (whether tied or not) to have multiple shared partners

mod13 = ergm(net ~ edges + nodefactor("type1")+ nodefactor("type2") + 
               gwesp(decay = 1, fixed = T) + 
               gwdsp(decay = 1, fixed = T), control = control.ergm(seed=1))
summary(mod13)
mcmc.diagnostics(mod13)



BIC(mod0, mod1, mod3, mod4, mod5, mod6, mod7, mod8, mod11, mod13)
AIC(mod0, mod1, mod3, mod4, mod5, mod6, mod7, mod8, mod11,  mod13)
#modello 11 è il preferito 

# let us remove the gwesp term
mod15 = ergm(net ~ edges + nodefactor("type1")+ nodefactor("type2") + 
               gwdsp(decay = 1,fixed = T), control = control.ergm(seed=1))
summary(mod15) 
mcmc.diagnostics(mod15)

BIC(mod0, mod1, mod3, mod4, mod5, mod6, mod7, mod8, mod11, mod13, mod15)
AIC(mod0, mod1, mod3, mod4, mod5, mod6, mod7, mod8, mod11, mod13, mod15)
#model 11

# let us evaluate the goodness of fit
# ----------------------------------
# simulate from the model
sim = simulate(mod0,nsim = 100, verbose = TRUE, seed = 1, 
               control=control.simulate.ergm(MCMC.burnin = 1000))
sim
?asIgraph
rho.sim = CId.sim = recip.sim = 0
fnc = function(xx){
  ig = asIgraph(xx)
  Y.sim = get.adjacency(ig, sparse = F)
  diag(Y.sim) = NA
  tr = transitivity(ig)
  ideg = sd(degree(ig, mode = "in"))
  odeg = sd(degree(ig, mode = "out"))
  
  rho.sim= mean(Y.sim, na.rm = TRUE)
  CId.sim = centr_degree(ig, mode = "in", loops = F)$centralization
  recip.sim = reciprocity(ig)
  
  return(c(tr, ideg, odeg, rho.sim, CId.sim, recip.sim))
}

null.distr = matrix(,100,6)
for(b in 1:100){
  null.distr[b,]  = fnc(sim[[b]])
}
dev.new()
par(mfrow = c(1,3))
hist(unlist(null.distr[,1]), xlab = "transitivity"); abline(v = transitivity(g), col = "red")
hist(unlist(null.distr[,2]), xlab = "in-degree"); abline(v = sd(degree(g, mode = "in")), col = "red")
hist(unlist(null.distr[,3]), xlab = "out-degree"); abline(v = sd(degree(g, mode = "out")), col = "red")

par(mfrow = c(1,3))
hist(unlist(null.distr[,4]), xlab = "density"); abline(v = rho.obs, col = "red")
hist(unlist(null.distr[,5]), xlab = "clustering coefficient"); abline(v = CId.obs, col = "red")
hist(unlist(null.distr[,6]), xlab = "reciprocity"); abline(v = reciprocity(g), col = "red")

mean(null.distr[,5]>CId.obs)
#cattura bene la transitività, l'in-degree ma non l'out degree perchè non è presente nel modello 11



AIC(mod0, mod1, mod2, mod11)
BIC(mod0, mod1, mod2, mod11)


###### SBM ######
#un.g = as.undirected(g, mode = "collapse")
#Yun = get.adjacency(un.g,sparse = F)
#diag(Yun) = NA

Y = get.adjacency(g,sparse = F)
diag(Y) = NA

plotMyMatrix(Y, dimLabels = list(row = 'neuron', col = 'neuron'), plotOptions = list(legend = TRUE))

sbm1 = estimateSimpleSBM(Y, "bernoulli", dimLabels = 'tree', 
                         estimOptions = list(verbosity = 1))

sbm1


# selected number of blocks
sbm1$nbBlocks
# prior block probabilities
sbm1$blockProp
# connectivity parameters
round(sbm1$connectParam$mean,3)

plot(sbm1, type = "data",  plotOptions = list(legend = TRUE))
# nodes are ordered wrt to the block they belong to and blocks are highlighted
plot(sbm1, type = "expected", plotOptions = list(legend = TRUE))
# fitted connection probabilities
plot(sbm1, type = "meso", plotOptions = list(legend = TRUE))
# fitted connection probabilities

# info on all estimated model is given in 
sbm1$storedModels
table(sbm1$memberships, V(g)$type1)
table(sbm1$memberships, V(g)$type2)



# let us consider poisson distributed variables Y_ij
sbm2 = estimateSimpleSBM(Y, "poisson", dimLabels = 'neuron', 
                         estimOptions = list(verbosity = 1))
sbm2

# selected number of blocks
sbm2$nbBlocks
# prior block probabilities
sbm2$blockProp
# connectivity parameters
round(sbm2$connectParam$mean,3)

# Let us graphically represent the data matrix 
# reordering rows and cols according to the estimated block in the SBM

plot(sbm2, type = "data", dimLabels = list(row = 'neuron', col= 'neuron'), plotOptions = list(legend = TRUE))

# or the average number of connections between trees
plot(sbm2, type = "expected", dimLabels = list(row = 'neuron', col= 'neuron'))

# or 
plot(sbm2, type = "meso", dimLabels = list(row = 'neuron', col= 'neuron'))

sbm2$storedModels
table(sbm2$memberships, V(g)$type1)
table(sbm2$memberships, V(g)$type2)
