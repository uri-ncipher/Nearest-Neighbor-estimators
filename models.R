library(lme4)
library(plyr)
library(dplyr)
library(igraph)
library(numDeriv)
library(gtools)
rm(list = ls())

source("https://github.com/uri-ncipher/Nearest-Neighbor-estimators/blob/main/functions.R")
#read in the network net0
nodes=read.csv("https://github.com/uri-ncipher/Nearest-Neighbor-estimators/blob/main/nodes.csv")
edges=read.csv("https://github.com/uri-ncipher/Nearest-Neighbor-estimators/blob/main/edges.csv")

net0=graph_from_data_frame(d=edges, vertices = nodes, directed=F)

#number of sample size, number of network component
n=length(V(net0))
m=components(net0)$no

data=data.frame(id=1:n, na=unlist(lapply(1:n, num_neighbors, net=net0)),
                component=components(net0)$membership)

#assign treatment, outcome, and baseline covariates to data
data$treatment=nodes$treatment
data$outcome=nodes$outcome
data$var1=nodes$var1
data$var2=nodes$var2
data$na_a= unlist(lapply(1:n, trt_neighbors, net=net0))
data$notna_a=data$na-data$na_a

base_covariate=c("var1", "var2")
#averaging baseline covariates
data$avg_var1=unlist(lapply(1:n, avg_neighbors, net=net0, variable="var1"))
data$avg_var2=unlist(lapply(1:n, avg_neighbors, net=net0, variable="var2"))

avg_covariate=c("avg_var1", "avg_var2")
##################################################################
############# IPW1 model #########################################
##################################################################

start_time=Sys.time()
IPW_1_model(data, base_covariate, c(0.25, 0.5, 0.75))
end_time=Sys.time()
end_time-start_time

##################################################################
############# IPW2 model #########################################
##################################################################

formula_1=paste("cbind(na_a, notna_a)", "~", "treatment", "+", 
                paste(base_covariate, collapse = "+"), "+", 
                paste(avg_covariate, collapse = "+"))
formula_2=paste("treatment", "~", paste(base_covariate, collapse = "+"))
M1=glm(formula_1, family = binomial(link = "logit"), data=data)
M2=glm(formula_2, family = binomial(link = "logit"), data=data)

start_time=Sys.time()
IPW_2_model(data, M1, M2, c(0.25, 0.5, 0.75))
end_time=Sys.time()
end_time-start_time
