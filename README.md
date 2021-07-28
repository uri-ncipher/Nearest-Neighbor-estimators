---
Title: "IPW models"
author: "TingFang Lee"
---

Load packages and function script...

```{r message=FALSE}
library(lme4)
library(plyr)
library(dplyr)
library(igraph)
library(numDeriv)
library(gtools)
rm(list = ls())

source("https://github.com/uri-ncipher/Nearest-Neighbor-estimators/blob/main/functions.R?raw=TRUE")
```


## Data Preparation

Read in the synthetic nodes and edges file for creating the simulated network.
```{r}
nodes=read.csv("https://github.com/uri-ncipher/Nearest-Neighbor-estimators/blob/main/nodes.csv?raw=TRUE")
edges=read.csv("https://github.com/uri-ncipher/Nearest-Neighbor-estimators/blob/main/edges.csv?raw=TRUE")
net0=graph_from_data_frame(d=edges, vertices = nodes, directed=F)
```

Network visuliation 
```{r}
plot(net0, vertex.size=1, vertex.label=NA, vertex.color="cadetblue3", vertex.frame.color="cadetblue3")
```

Make the data for modeling
```{r}
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
```

## Modeling

Using IPW1 to evaluate the average potential outcome and causal effects under allocation strategies $\alpha$. The model will output the point estimation and the estimated variance of average potention outcomes $\widehat{Y}(1, \alpha), \widehat{Y}(0, \alpha), \widehat{Y}(\alpha)$.

```{r}
alpha=c(0.25, 0.5, 0.75)
IPW_1_model(data, base_covariate, alpha)
```


###IPW 2
Using IPW2 to evaluate the average potential outcome and causal effects under allocation strategies $\alpha$. The model will output the point estimation and the estimated variance of average potention outcomes $\widehat{Y}(1, \alpha), \widehat{Y}(0, \alpha), \widehat{Y}(\alpha)$.

```{r}
formula_1=paste("cbind(na_a, notna_a)", "~", "treatment", "+", 
                paste(base_covariate, collapse = "+"), "+", 
                paste(avg_covariate, collapse = "+"))
formula_2=paste("treatment", "~", paste(base_covariate, collapse = "+"))
M1=glm(formula_1, family = binomial(link = "logit"), data=data)
M2=glm(formula_2, family = binomial(link = "logit"), data=data)

IPW_2_model(data, M1, M2, alpha)

```


