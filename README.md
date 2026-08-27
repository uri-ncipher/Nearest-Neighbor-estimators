---
Title: "IPW models"
author: "TingFang Lee, URI Avenir Team"
---

Load packages and function script...

```{r message=FALSE}
library(lme4)
library(plyr)
library(igraph)
library(numDeriv)
library(dplyr)

source("https://github.com/uri-ncipher/Nearest-Neighbor-estimators/blob/main/IPW.R?raw=TRUE")
```


## Data Preparation

Read in the synthetic data and network structure
```{r}
sim.data <- read.csv("https://github.com/uri-ncipher/Nearest-Neighbor-estimators/blob/main/sim.data.csv?raw=TRUE")
net0 <- readRDS(file = url("https://github.com/uri-ncipher/Nearest-Neighbor-estimators/blob/main/net0.100.rds?raw=TRUE", open = "rb"))
```

### Make the data for modeling
The synthetic dataset `sim.data` and network `net0` include observations with missing outcomes. These data can be analyzed using the IPCW estimators, which account for outcome missingness.

Because the IPW estimators require complete outcome data, we first calculate each individual's number of neighbors and number of treated neighbors in the original network. We then remove individuals with missing outcomes from the network, along with any resulting isolated nodes, to create the complete-case network `net1`. The number of neighbors and number of treated neighbors are then recalculated based on this reduced network for use in the IPW analysis.

Finally, for both the original and complete-case networks, we calculate the average covariate value among each individual's neighbors (`avg_L`), which is used in the subsequent IPCW2 and IPW2.
```{r}
V(net0)$name <- as.character(V(net0))
remove.ids <- as.character(sim.data$id[sim.data$censored.cor == 1])
net1 <- delete_vertices(net0, remove.ids)
net1 <- delete_vertices(net1, V(net1)[degree(net1) == 0])

keep.ids <- as.numeric(V(net1)$name)

sim.data1 <- sim.data[match(keep.ids, sim.data$id), ]

sim.data1$n.neighbor <- degree(net1)

sim.data1$n.trt.neighbor <- sapply(sim.data1$id, function(i) {
  nb <- neighbors(net1, as.character(i))
  nb.ids <- as.numeric(V(net1)[nb]$name)
  sum(sim.data1$A[match(nb.ids, sim.data1$id)])
})

nn.sim.data <- data.frame()
for (i in sim.data[["id"]]) {
  nb <- neighborhood(net0, order = 1, nodes = as.character(i))[[1]]
  ids <- as.numeric(V(net0)[nb]$name)
  ids <- setdiff(ids, i)
  temp.data <- sim.data[sim.data[["id"]] %in% ids, ]
  temp.data$ego <- i
  nn.sim.data <- rbind(nn.sim.data, temp.data)
}

nn.sim.data <- nn.sim.data %>%
  group_by(ego) %>%
  dplyr::summarise(avg_L = mean(L), .groups = "drop")
sim.data <- merge(sim.data, nn.sim.data, by.x = "id", by.y = "ego")

nn.sim.data.1 <- data.frame()
for (i in sim.data1[["id"]]) {
  nb <- neighborhood(net1, order = 1, nodes = as.character(i))[[1]]
  ids <- as.numeric(V(net1)[nb]$name)
  ids <- setdiff(ids, i)
  temp.data <- sim.data1[sim.data1[["id"]] %in% ids, ]
  temp.data$ego <- i
  nn.sim.data.1 <- rbind(nn.sim.data.1, temp.data)
}

nn.sim.data.1 <- nn.sim.data.1 %>%
  group_by(ego) %>%
  dplyr::summarise(avg_L = mean(L), .groups = "drop")
sim.data1 <- merge(sim.data1, nn.sim.data.1, by.x = "id", by.y = "ego")
```

## Modeling

In this section, we apply the IPW1, IPW2, IPCW1, and IPCW2 estimators to estimate average potential outcomes and the corresponding causal effects. The average potential outcome is defined by an individual's exposure/treatment status and the allocation strategy (\alpha), where (\alpha) represents the probability that individuals in the nearest-neighbor set receive the exposure/treatment under the counterfactual allocation strategy. Thus, (\alpha) takes a value between 0 and 1.

For each estimator, the R functions return the **point estimates and estimated variances of the average potential outcomes**, (\boldsymbol{\widehat{Y}(1,\alpha), \widehat{Y}(0,\alpha), \widehat{Y}(\alpha)}), under three allocation strategies, (\alpha = 0.25, 0.50,) and (0.75). The functions also return the **point estimates and estimated variances of four causal effects**: the direct effect (DE), indirect effect (IE), total effect (TE), and overall effect (OE).

The IPW1 and IPW2 estimators are applied to complete outcome data, whereas the IPCW1 and IPCW2 estimators accommodate missing outcomes through inverse probability of censoring weighting.

Equations for calculating the four causal effects are:
  
* ${\widehat{DE} = \widehat{Y}(1, \alpha) - \widehat{Y}(0, \alpha)}$
* ${\widehat{IE} = \widehat{Y}(0, \alpha_0) - \widehat{Y}(0, \alpha_1)}$
* ${\widehat{TE} = \widehat{Y}(1, \alpha_0) - \widehat{Y}(0, \alpha_1)}$
* ${\widehat{OE} = \widehat{Y}(\alpha_0) - \widehat{Y}(\alpha_1)}$,
  
where $\alpha_0$, $\alpha_1$ both represents allocation strategies, and $\alpha_0 \neq  \alpha_1$ (Note: The associated [paper](https://arxiv.org/abs/2108.04865) compared the estimated average potential outcomes with $\alpha_1$ to $\alpha_0$ for the indirect, total, and overall effects. The code here compared the estimated average potential outcomes with $\alpha_0$ to $\alpha_1$ for the indirect, total, and overall effects). Users can set any values between 0 and 1 for $\alpha_0$ and  $\alpha_1$. If given $\alpha$ is a vector, then the causal contrasts will produce pairwise comparisons in the sequential order of the values in $\alpha$. For example, if $\alpha = c(0.75, 0.25, 0.50)$, the contrasts for spillover effects will compare 
  
* (i) $\alpha_0 = 0.75$ versus $\alpha_1 = 0.25$; 
* (ii) $\alpha_0 = 0.75$ versus $\alpha_1 = 0.50$;
* (iii) $\alpha_0 = 0.25$ versus $\alpha_1 = 0.50$.

## Function Arguments

The following arguments are used across the IPW and IPCW point and variance estimation functions.

- **`data`**: A data frame containing the individual-level analysis data.

- **`net`**: An `igraph` network object representing the observed network structure. Vertex names should correspond to individual IDs in `data`.

- **`alphas`**: A numeric vector of hypothetical treatment allocation probabilities used to define the target treatment allocation strategies. Values should be between 0 and 1.

- **`id.var`**: A character string specifying the name of the variable that uniquely identifies each individual.

- **`trt.var`**: A character string specifying the name of the binary individual treatment indicator variable, typically coded as 0 = untreated and 1 = treated.

- **`n.neighbor`**: A character string specifying the name of the variable containing the total number of network neighbors for each individual.

- **`n.trt.neighbor`**: A character string specifying the name of the variable containing the number of treated network neighbors for each individual.

- **`propensity.covars`**: A character vector specifying the names of individual-level covariates included in the treatment propensity model.

- **`propensity.nn.covars`**: A character vector specifying the names of covariates included in the model for treatment among an individual's network neighbors. This argument is used in the type 2 estimators.

- **`censor.indicator`**: A character string specifying the name of the binary censoring indicator variable, typically coded as 1 = censored and 0 = observed.

- **`censor.covars`**: A character vector specifying the names of covariates included in the censoring model.

- **`censor.correlation`**: A logical value indicating whether correlation in censoring within network components should be accounted for.

- **`component.var`**: A character string specifying the name of the variable identifying network components. This variable is used to account for clustering within connected network components, including in variance estimation and correlated censoring models.

- **`obs.y`**: A character string specifying the name of the variable containing the observed outcome.

```{r}
point.est.IPCW1 = Y_IPCW1(data = sim.data, net = net0, alphas = c(0.25, 0.5, 0.75), id.var = "id", trt.var = "A", n.neighbor = "n.neighbor", n.trt.neighbor = "n.trt.neighbor", propensity.covars = c("L"), censor.indicator = "censored.cor", censor.covars = c("L"), censor.correlation = TRUE, component.var = "component", obs.y = "obs.y.cor")
var.est.IPCW1 = Var_Y_IPCW1(alphas = c(0.25, 0.5, 0.75), net = net0, data = sim.data, id.var = "id", obs.y= "obs.y.cor", n.neighbor = "n.neighbor", n.trt.neighbor = "n.trt.neighbor", trt.var = "A", propensity.covars = c("L"), censor.indicator = "censored.cor", censor.covars = c("L"), censor.correlation = TRUE, component.var = "component")

point.est.IPCW2 = Y_IPCW2(data = sim.data, net = net0, alphas = c(0.25, 0.5, 0.75), id.var = "id", trt.var = "A", n.neighbor = "n.neighbor", n.trt.neighbor = "n.trt.neighbor", propensity.covars = c("L"), censor.covars = c("L"), censor.correlation = TRUE, component.var = "component", obs.y = "obs.y.cor")
var.est.IPCW2 = Var_Y_IPCW2(alphas = c(0.25, 0.5, 0.75), net = net0, data = sim.data, id.var = "id", obs.y= "obs.y.cor", n.neighbor = "n.neighbor", n.trt.neighbor = "n.trt.neighbor", trt.var = "A", propensity.covars = c("L"), propensity.nn.covars = c("avg_L"), censor.indicator = "censored.cor", censor.covars = c("L"), censor.correlation = TRUE, component.var = "component")

point.est.IPW1 = Y_IPW1(data = sim.data1, net = net1, alphas = c(0.25, 0.5, 0.75), id.var = "id", trt.var = "A", n.neighbor = "n.neighbor", n.trt.neighbor = "n.trt.neighbor", propensity.covars = c("L"), obs.y = "obs.y.cor")
var.est.IPW1 = Var_Y_IPW1(alphas = c(0.25, 0.5, 0.75), net = net1, data = sim.data1, id.var = "id", obs.y= "obs.y.cor", n.neighbor = "n.neighbor", n.trt.neighbor = "n.trt.neighbor", trt.var = "A", propensity.covars = c("L"), component.var = "component")

point.est.IPW2 = Y_IPW2(data = sim.data1, net = net1, alphas = c(0.25, 0.5, 0.75), id.var = "id", trt.var = "A", n.neighbor = "n.neighbor", n.trt.neighbor = "n.trt.neighbor", propensity.nn.covars = c("avg_L"), propensity.covars = c("L"), obs.y = "obs.y.cor")
var.est.IPW2 = Var_Y_IPW2(alphas = c(0.25, 0.5, 0.75), net = net1, data = sim.data1, id.var = "id", obs.y= "obs.y.cor", n.neighbor = "n.neighbor", n.trt.neighbor = "n.trt.neighbor", trt.var = "A", propensity.covars = c("L"), propensity.nn.covars = c("avg_L"), component.var = "component")
```
## References
1. Lee T, Buchanan AL, Katenka NV, Forastiere L, Halloran ME, Friedman SR, Nikolopoulos G. Estimating causal effects of HIV prevention interventions with interference in network-based studies among people who inject drugs. Annals of Applied Statistics. 2023;17(3):2165–2191. doi:10.1214/22-AOAS1713.
2. Lee T, Buchanan AL, Katenka N, Forastiere L, Halloran ME, Nikolopoulos G. Assessing spillover effects: Handling missing outcomes in network-based studies. Statistical Methods in Medical Research. 2025;34(12):2284–2301. doi:10.1177/09622802251382586.
