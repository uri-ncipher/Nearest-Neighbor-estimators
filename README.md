---
Title: "IPW models"
author: "TingFang Lee, URI Avenir Team"
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

Read in the synthetic nodes and edges data sets for creating the simulated network.
```{r}
nodes=read.csv("https://github.com/uri-ncipher/Nearest-Neighbor-estimators/blob/main/nodes.csv?raw=TRUE")
edges=read.csv("https://github.com/uri-ncipher/Nearest-Neighbor-estimators/blob/main/edges.csv?raw=TRUE")
net0=graph_from_data_frame(d=edges, vertices = nodes, directed=F)  # create network based on the nodes and edges data.
```

### Network visualization 
```{r}
# plot the network generated above
plot(net0, vertex.size=1, vertex.label=NA, vertex.color="cadetblue3", vertex.frame.color="cadetblue3")
```

### Make the data for modeling
```{r}
# Save the number of individuals (nodes) and number of components
n=length(V(net0))       # number of individuals
m=components(net0)$no   # number of components

# generate a dataset with three columns, the first column is participant id (node id), 
# the second column is number of nearest neighbors for each participant (is also node degree), 
# the third column indicate which component the participant is in.
data=data.frame(id=1:n, na=unlist(lapply(1:n, num_neighbors, net=net0)),
                component=components(net0)$membership) 

# assign treatment, outcome, and baseline covariates to data
data$treatment=nodes$treatment
data$outcome=nodes$outcome
data$var1=nodes$var1
data$var2=nodes$var2
data$na_a= unlist(lapply(1:n, trt_neighbors, net=net0))
data$notna_a=data$na-data$na_a

base_covariate=c("var1", "var2")

# averaging baseline covariates (average of nearest neighbors' baseline covariate values)
data$avg_var1=unlist(lapply(1:n, avg_neighbors, net=net0, variable="var1"))
data$avg_var2=unlist(lapply(1:n, avg_neighbors, net=net0, variable="var2"))

avg_covariate=c("avg_var1", "avg_var2")
```

## Modeling

In this section, we will apply both IPW1 and IPW2 estimators for the average potential outcomes and corresponding causal effects. An individual's exposure/treatment and allocation strategy $\alpha$ jointly determine the average potential outcome (allocation strategy represents the probability that individuals in nearest neighbors set receiving the exposure/treatment in the counterfactual scenario, so it is a numeric value between 0 and 1). For each estimator, the R programs will output the **point estimation and the estimated variance of average potential outcomes** $\boldsymbol{\widehat{Y}(1, \alpha), \widehat{Y}(0, \alpha), \widehat{Y}(\alpha)}$ under three different allocation strategies $\alpha$ (i.e., 0.25, 0.50, and 0.75), respectively. Also, the R program will output the **point estimation and estimated variance of four causal effects: Direct** (DE), **Indirect** (IE), **Total** (TE), and **Overall effect** (OE).  

Equations for calculating the four causal effects are:
  
* ${\widehat{DE} = \widehat{Y}(1, \alpha) - \widehat{Y}(0, \alpha)}$
* ${\widehat{IE} = \widehat{Y}(0, \alpha_0) - \widehat{Y}(0, \alpha_1)}$
* ${\widehat{TE} = \widehat{Y}(1, \alpha_0) - \widehat{Y}(0, \alpha_1)}$
* ${\widehat{OE} = \widehat{Y}(\alpha_0) - \widehat{Y}(\alpha_1)}$,
  
where $\alpha_0$, $\alpha_1$ both represents allocation strategies, and $\alpha_0 \neq  \alpha_1$ (Note: The associated [paper](https://arxiv.org/abs/2108.04865) compared the estimated average potential outcomes with $\alpha_1$ to $\alpha_0$ for the indirect, total, and overall effects. The code here compared the estimated average potential outcomes with $\alpha_0$ to $\alpha_1$ for the indirect, total, and overall effects). Users can set any values between 0 and 1 for $\alpha_0$ and  $\alpha_1$.

### IPW1
Using IPW1 to estimate the average potential outcomes and causal effects under allocation strategies $\alpha$.


```{r}
alpha=c(0.25, 0.5, 0.75)  # set allocation strategies
IPW_1_model(data, base_covariate, alpha)
```
In the output above, 
  
* The **section [[1]]** displays the point estimates (type = "point estimate" in the output) and estimated variances (type = "variance" in the output) for the average potential outcomes under three different allocation strategies (i.e., $\alpha = 0.25, 0.50, 0.75$) and each individual exposure ($a = 0, a=1$, and margin). "Margin" is the average potential outcomes for a particular allocation strategy, regardless of individual exposure status.
* The **section [[2]]** displays the point estimates (type = "Direct", "Indirect", "Total", "Overall" in the output) and estimated variances (type = "Var DE", "Var IE", "Var TE", "Var OE" in the output) of four causal effects: direct, indirect, total and overall effects, corresponding to particular allocation strategies $(\alpha_0$ and $\alpha_1)$. The estimated values are shown in "estimation" column in the output above.


### IPW2
Using IPW2 to estimate the average potential outcomes and causal effects under allocation strategies $\alpha$. 

```{r}
formula_1=paste("cbind(na_a, notna_a)", "~", "treatment", "+", 
                paste(base_covariate, collapse = "+"), "+", 
                paste(avg_covariate, collapse = "+"))
formula_2=paste("treatment", "~", paste(base_covariate, collapse = "+"))
M1=glm(formula_1, family = binomial(link = "logit"), data=data)
M2=glm(formula_2, family = binomial(link = "logit"), data=data)

IPW_2_model(data, M1, M2, alpha)

```

In the output above, 
  
* The **section [[1]]** displays the point estimates (type = "point estimate" in the output) and estimated variances (type = "variance" in the output) for the average potential outcomes under three different allocation strategies (i.e., $\alpha = 0.25, 0.50, 0.75$) and each individual exposure ($a = 0, a=1$, and margin). "Margin" is the average potential outcomes for a particular allocation strategy, regardless of individual exposure status.
* The **section [[2]]** displays the point estimates (type = "Direct", "Indirect", "Total", "Overall" in the output) and estimated variances (type = "Var DE", "Var IE", "Var TE", "Var OE" in the output) of four causal effects: direct, indirect, total and overall effects, corresponding to particular allocation strategies $(\alpha_0$ and $\alpha_1)$. The estimated values are shown in "estimation" column in the output above.

