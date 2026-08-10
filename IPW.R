##############################################################
##### auxiliary functions ####################################
##############################################################

propensity.1 <- function(data, net, id.var, trt.var, propensity.covars) {
  
  nn.data <- data.frame()
  for (i in data[[id.var]]) {
    nb <- neighborhood(net, order = 1, nodes = as.character(i))[[1]]
    ids <- as.numeric(V(net)[nb]$name)
    temp.data <- data[data[[id.var]] %in% ids, ]
    temp.data$ego <- i
    nn.data <- rbind(nn.data, temp.data)
  }
  
  form <- as.formula(paste(trt.var, "~", paste(propensity.covars, collapse = " + "), "+ (1 | ego)"))
  propensity.model <- lme4::glmer(form, data = nn.data, family = binomial(link = "logit"))
  theta <- c(lme4::fixef(propensity.model), as.data.frame(lme4::VarCorr(propensity.model))$sdcor[1])
  
  sapply(data[[id.var]], function(i) {
    nb <- neighborhood(net, order = 1, nodes = as.character(i))[[1]]
    ids <- as.numeric(V(net)[nb]$name)
    temp.data <- data[data[[id.var]] %in% ids, ]
    X <- cbind(1, temp.data[, propensity.covars])
    integrand <- function(b) {
      prob = 1
      for (each.row in 1:nrow(temp.data)){
        prob.new <- plogis(as.vector(X[each.row, ] %*% as.matrix(theta[-length(theta)])) + b)
        if (temp.data[[trt.var]][each.row] == 1){
          prob <- prob*prob.new
        } else {prob <- prob*(1-prob.new)}
      }
      return(prob * dnorm(b, mean = 0, sd = theta[length(theta)]))
    }
    integrate(integrand, -Inf, Inf)$value
  })
}

propensity.2 <- function(data, trt.var, n.trt.neighbor, n.neighbor,
                         propensity.covars, propensity.nn.covars) {
  data$not.trt.neighbor <- data[[n.neighbor]] - data[[n.trt.neighbor]] 
  form.ind <- as.formula(paste(trt.var, "~", paste(propensity.covars, collapse = " + ")))
  form.nn <- as.formula(paste("cbind(", n.trt.neighbor, ", not.trt.neighbor)~", 
                              trt.var, "+", paste(propensity.covars, collapse = " + "), 
                              "+", paste(propensity.nn.covars, collapse = " + ")))
  propensity.model.ind <- glm(form.ind, data = data, family = binomial(link = "logit"))
  propensity.model.nn <- glm(form.nn, data = data, family = binomial(link = "logit"))
  
  p.ind = predict(propensity.model.ind, newdata = data, type = "response")
  p.nn = predict(propensity.model.nn, newdata = data, type = "response")
  return(dbinom(data[[trt.var]], size = 1, prob = p.ind)*(p.nn^data[[n.trt.neighbor]])*((1-p.nn)^data$not.trt.neighbor))
}

censor <- function(data, censor.indicator, censor.covars,
                   censor.correlation = FALSE, component.var = NULL) {
  if (censor.correlation) {
    form <- as.formula(paste(censor.indicator, "~",paste(censor.covars, collapse = " + "), 
                             "+ (1 |", component.var, ")", sep = ""))
    censor.model <- lme4::glmer(form, data = data, family = binomial(link = "logit"))
  } else {
    form <- as.formula(paste(censor.indicator, "~", paste(censor.covars, collapse = " + ")))
    censor.model <- glm(form, data = data, family = binomial(link = "logit"))
  }
  1 - predict(censor.model, newdata = data, type = "response")
}

pi.neighbor <- function(data, n.neighbor, n.trt.neighbor, alpha) {
  alpha^data[[n.trt.neighbor]] *
    (1 - alpha)^(data[[n.neighbor]] - data[[n.trt.neighbor]])
}

##################################################################################
################## IPW1 functions #################################################
##################################################################################

Y_IPW1 <- function(data, net, alphas, id.var, trt.var, n.neighbor, n.trt.neighbor,
                   propensity.covars, obs.y) {
  #p <- pi.neighbor(data = data, n.neighbor = n.neighbor, n.trt.neighbor = n.trt.neighbor,
  #                 alpha = alpha)
  p <- lapply(alphas, function(a)
    pi.neighbor(data, n.neighbor, n.trt.neighbor, alpha = a))
  
  propensity.score <- propensity.1(data = data, net = net, id.var = id.var, 
                                   trt.var = trt.var, propensity.covars = propensity.covars)
  Y_IPW1_1 <- sapply(p, function(pp) {mean(data[[obs.y]] * data[[trt.var]] * pp / propensity.score)})
  Y_IPW1_0 <- sapply(p, function(pp) {mean(data[[obs.y]] * (1 - data[[trt.var]]) * pp / propensity.score)})
  Y_IPW1_margin <- mapply(function(pp, a) {mean(data[[obs.y]] * dbinom(data[[trt.var]], 1, a) * pp / propensity.score)}, p, alphas)
  avg.outcome <- cbind(Y_IPW1_1, Y_IPW1_0, Y_IPW1_margin)
  colnames(avg.outcome) <- c("treat = 1", "treat = 0", "margin")
  rownames(avg.outcome) <- paste0("alpha = ", alphas)
  pairs <- combn(1:length(alphas), 2)
  DE <- data.frame()
  for (k in 1:length(alphas)){
    alpha1 <- alphas[k]
    alpha2 <- alphas[k]
    estimate <- avg.outcome[k, "treat = 1"] - avg.outcome[k, "treat = 0"]
    DE <- rbind(DE, c(alpha1, alpha2, estimate))
  }
  names(DE) <- c("alpha1", "alpha2", "estimate")
  IE <- data.frame()
  for (k in 1:ncol(pairs)){
    alpha1 <- alphas[pairs[1, k]]
    alpha2 <- alphas[pairs[2, k]]
    estimate <- avg.outcome[pairs[2, k], "treat = 0"] - avg.outcome[pairs[1, k], "treat = 0"]
    IE <- rbind(IE, c(alpha1, alpha2, estimate))
  }
  names(IE) <- c("alpha1", "alpha2", "estimate")
  TE <- data.frame()
  for (k in 1:ncol(pairs)){
    alpha1 <- alphas[pairs[1, k]]
    alpha2 <- alphas[pairs[2, k]]
    estimate <- avg.outcome[pairs[2, k], "treat = 1"] - avg.outcome[pairs[1, k], "treat = 0"]
    TE <- rbind(TE, c(alpha1, alpha2, estimate))
  }
  names(TE) <- c("alpha1", "alpha2", "estimate")
  DE$effect <- "DIRECT"
  IE$effect <- "INDIRECT"
  TE$effect <- "TOTAL"
  contrast <- rbind(DE, IE, TE)
  contrast <- contrast[, c(4, 1, 2, 3)]
  return(list(avg.outcome, contrast))
}

Var_Y_IPW1 <- function(net, data, alphas, id.var, obs.y, n.neighbor, n.trt.neighbor, 
                       trt.var, propensity.covars, component.var){
  n <- nrow(data)
  m <- length(unique(data[[component.var]]))
  
  p <- lapply(alphas, function(a)
    pi.neighbor(data, n.neighbor, n.trt.neighbor, alpha = a))
  propensity.score <- propensity.1(data = data, net = net, id.var = id.var, 
                                   trt.var = trt.var, propensity.covars = propensity.covars)
  
  nn.data <- data.frame()
  for (i in data[[id.var]]) {
    nb <- neighborhood(net, order = 1, nodes = as.character(i))[[1]]
    ids <- as.numeric(V(net)[nb]$name)
    temp.data <- data[data[[id.var]] %in% ids, ]
    temp.data$ego <- i
    nn.data <- rbind(nn.data, temp.data)
  }
  form <- as.formula(paste(trt.var, "~", paste(propensity.covars, collapse = " + "), "+ (1 | ego)"))
  propensity.model <- lme4::glmer(form, data = nn.data, family = binomial(link = "logit"))
  theta <- c(lme4::fixef(propensity.model), as.data.frame(lme4::VarCorr(propensity.model))$sdcor[1])
  
  
  propen_U_11 <- matrix(0, nrow = length(theta), ncol = length(theta))
  propen_v_11 <- matrix(NA, nrow = n, ncol = length(theta))
  grad_propensity_inv <- matrix(NA, nrow = n, ncol = length(theta))
  for (i in data[[id.var]]){
    
    #####################################################################
    ## individual propensity and censoring functions given individual i##
    #####################################################################
    propensity_ind <- function(theta) {
      nb <- neighborhood(net, order = 1, nodes = as.character(i))[[1]]
      ids <- as.numeric(V(net)[nb]$name)
      temp.data <- data[data[[id.var]] %in% ids, ]
      X <- model.matrix(as.formula(paste("~", paste(propensity.covars, collapse = " + "))),
                        data = temp.data)
      beta <- theta[-length(theta)]
      sigma.b <- theta[length(theta)]
      X <- cbind(1, temp.data[, propensity.covars])
      integrand <- function(b) {
        prob = 1
        for (each.row in 1:nrow(temp.data)){
          prob.new <- plogis(as.vector(X[each.row, ] %*% as.matrix(theta[-length(theta)])) + b)
          if (temp.data[[trt.var]][each.row] == 1){
            prob <- prob*prob.new
          } else {prob <- prob*(1-prob.new)}
        }
        return(prob * dnorm(b, mean = 0, sd = theta[length(theta)]))
      }
      integrate(integrand, -Inf, Inf)$value
    }
    propensity_ind_inv <- function(theta) {1 / propensity_ind(theta)}
    log_propensity_ind <- function(theta){return(log(propensity_ind(theta)))}
    
    #####################################################################
    ## Construct matrices ###############################################
    #####################################################################
    
    propen_U_11 <- propen_U_11-(1/n)*hessian(log_propensity_ind, x = theta)
    propen_v_11[data[[id.var]] == i, ] <- grad(log_propensity_ind, x = theta)
    grad_propensity_inv[data[[id.var]] == i, ] <- grad(propensity_ind_inv, x = theta)
  }
  U_11 <- propen_U_11
  
  propen_V_11 <- matrix(0, nrow = length(theta), ncol = length(theta))
  for (j in unique(data[[component.var]])){
    idx <- which(data[[component.var]] == j)
    propen_sum <- colSums(propen_v_11[idx, , drop = FALSE])
    propen_V_11 <- propen_V_11 + propen_sum %*% t(propen_sum)
  }
  
  V_11=(m/n^2)*propen_V_11
  
  Var_avg.outcome <- matrix(NA, nrow = length(alphas), ncol = 3)
  for (k in 1:length(alphas)){
    total.par <- length(theta)+1
    U_21.1 <- -colMeans((data[[obs.y]] * data[[trt.var]] * p[[k]]) * grad_propensity_inv)
    B.1 <- data[[obs.y]] * data[[trt.var]] * p[[k]] / propensity.score
    Y.1 <- mean(B.1)
    U_21.0 <- -colMeans((data[[obs.y]]*(1-data[[trt.var]])*p[[k]])*grad_propensity_inv)
    B.0 <- data[[obs.y]]*(1-data[[trt.var]])*p[[k]]/propensity.score
    Y.0 <- mean(B.0)
    p.mar <- dbinom(data[[trt.var]], 1, alphas[[k]])* p[[k]]
    U_21 <- -colMeans((data[[obs.y]]*p.mar)*grad_propensity_inv)
    B <- data[[obs.y]] * dbinom(data[[trt.var]], 1, alphas[[k]]) * p[[k]] / propensity.score
    Y <- mean(B)
    
    v_22.1 <- 0
    v_21.1 <- rep(0, length(theta))
    v_22.0 <- 0
    v_21.0 <- rep(0, length(theta))
    v_22 <- 0
    v_21 <- rep(0, length(theta))
    for (j in unique(data[[component.var]])) {
      idx <- which(data[[component.var]] == j)
      weight <- m / n
      r.1 <- weight * sum(B.1[idx]) - Y.1
      r.0 <- weight * sum(B.0[idx]) - Y.0
      r   <- weight * sum(B[idx])   - Y
      v_21.1 <- v_21.1 + weight * colSums(propen_v_11[idx, , drop = FALSE]) * r.1
      v_22.1 <- v_22.1 + r.1^2
      
      v_21.0 <- v_21.0 + weight * colSums(propen_v_11[idx, , drop = FALSE]) * r.0
      v_22.0 <- v_22.0 + r.0^2
      
      v_21 <- v_21 +  weight * colSums(propen_v_11[idx, , drop = FALSE]) * r
      v_22 <- v_22 + r^2
    }
    U <- cbind(U_11, rep(0, length(theta)))
    U.1 <- rbind(U, c(U_21.1, 1))
    U.0 <- rbind(U, c(U_21.0, 1))
    V.1 <- cbind(V_11, v_21.1/m)
    V.1 <- rbind(V.1, c(v_21.1, v_22.1)/m)
    V.0 <- cbind(V_11, v_21.0/m)
    V.0 <- rbind(V.0, c(v_21.0, v_22.0)/m)
    U_inv.1 <- solve(U.1)
    U_inv.0 <- solve(U.0)
    M.1 <- U_inv.1%*%V.1%*%t(U_inv.1)
    M.0 <- U_inv.0%*%V.0%*%t(U_inv.0)
    U <- rbind(U, c(U_21, 1))
    V <- cbind(V_11, v_21/m)
    V <- rbind(V, c(v_21, v_22)/m)
    U_inv <- solve(U)
    M <- U_inv%*%V%*%t(U_inv)
    var <- c(M.1[total.par, total.par] / m, M.0[total.par, total.par] / m, M[total.par, total.par] / m)
    Var_avg.outcome[k, ] <- var
  }
  colnames(Var_avg.outcome) <- c("treat = 1", "treat = 0", "margin")
  rownames(Var_avg.outcome) <- paste0("alpha = ", alphas)
  
  
  Var_DE <- data.frame()
  for (k in 1:length(alphas)){
    U_21 <- -colMeans((data[[obs.y]]*(data[[trt.var]])*p[[k]])*grad_propensity_inv - 
                        (data[[obs.y]]*(1-data[[trt.var]])*p[[k]])*grad_propensity_inv)
    B <- data[[obs.y]]*data[[trt.var]]*p[[k]]/propensity.score - 
      data[[obs.y]]*(1-data[[trt.var]])*p[[k]]/propensity.score
    Y <- mean(B)
    
    v_22 <- 0
    v_21 <- rep(0, length(theta))
    for (j in unique(data[[component.var]])) {
      idx <- which(data[[component.var]] == j)
      weight <- m / n
      r   <- weight * sum(B[idx]) - Y
      v_21 <- v_21 + weight * colSums(propen_v_11[idx, , drop = FALSE]) * r
      v_22 <- v_22 + r^2
    }
    U <- cbind(U_11, rep(0, length(theta)))
    U <- rbind(U, c(U_21, 1))
    V <- cbind(V_11, v_21/m)
    V <- rbind(V, c(v_21, v_22)/m)
    U_inv <- solve(U)
    M <- U_inv%*%V%*%t(U_inv)
    total.par <- length(theta)+1
    Var_DE <- rbind(Var_DE, c(alphas[k], alphas[k], M[total.par, total.par] / m))
  }
  names(Var_DE) <- c("alpha1", "alpha2", "Variance")
  
  
  pairs <- combn(1:length(alphas), 2)
  
  Var_IE <- data.frame()
  for (k in 1:length(alphas)){
    U_21 <- -colMeans((data[[obs.y]]*(1-data[[trt.var]])*p[[pairs[2, k]]])*grad_propensity_inv - 
                        (data[[obs.y]]*(1-data[[trt.var]])*p[[pairs[1, k]]])*grad_propensity_inv)
    B <- data[[obs.y]]*(1-data[[trt.var]])*p[[pairs[2, k]]]/propensity.score - 
      data[[obs.y]]*(1-data[[trt.var]])*p[[pairs[1, k]]]/propensity.score
    Y <- mean(B)
    
    v_22 <- 0
    v_21 <- rep(0, length(theta))
    for (j in unique(data[[component.var]])) {
      idx <- which(data[[component.var]] == j)
      weight <- m / n
      r   <- weight * sum(B[idx]) - Y
      v_21 <- v_21 + weight * colSums(propen_v_11[idx, , drop = FALSE]) * r
      v_22 <- v_22 + r^2
    }
    U <- cbind(U_11, rep(0, length(theta)))
    U <- rbind(U, c(U_21, 1))
    V <- cbind(V_11, v_21/m)
    V <- rbind(V, c(v_21, v_22)/m)
    U_inv <- solve(U)
    M <- U_inv%*%V%*%t(U_inv)
    total.par <- length(theta)+1
    Var_IE <- rbind(Var_IE, c(alphas[pairs[1, k]], alphas[pairs[2, k]], M[total.par, total.par] / m))
  }
  names(Var_IE) <- c("alpha1", "alpha2", "Variance")
  
  Var_TE <- data.frame()
  for (k in 1:length(alphas)){
    U_21 <- -colMeans((data[[obs.y]]*(data[[trt.var]])*p[[pairs[2, k]]])*grad_propensity_inv - 
                        (data[[obs.y]]*(1-data[[trt.var]])*p[[pairs[1, k]]])*grad_propensity_inv)
    B <- data[[obs.y]]*(data[[trt.var]])*p[[pairs[2, k]]]/propensity.score - 
      data[[obs.y]]*(1-data[[trt.var]])*p[[pairs[1, k]]]/propensity.score
    Y <- mean(B)
    
    v_22 <- 0
    v_21 <- rep(0, length(theta))
    for (j in unique(data[[component.var]])) {
      idx <- which(data[[component.var]] == j)
      weight <- m / n
      r   <- weight * sum(B[idx]) - Y
      v_21 <- v_21 + weight * colSums(propen_v_11[idx, , drop = FALSE]) * r
      v_22 <- v_22 + r^2
    }
    U <- cbind(U_11, rep(0, length(theta)))
    U <- rbind(U, c(U_21, 1))
    V <- cbind(V_11, v_21/m)
    V <- rbind(V, c(v_21, v_22)/m)
    U_inv <- solve(U)
    M <- U_inv%*%V%*%t(U_inv)
    total.par <- length(theta)+1
    Var_TE <- rbind(Var_TE, c(alphas[pairs[1, k]], alphas[pairs[2, k]], M[total.par, total.par] / m))
  }
  names(Var_TE) <- c("alpha1", "alpha2", "Variance")
  
  Var_DE$effect <- "DIRECT"
  Var_IE$effect <- "INDIRECT"
  Var_TE$effect <- "TOTAL"
  
  Var_contrast <- rbind(Var_DE, Var_IE, Var_TE)
  Var_contrast <- Var_contrast[, c(4, 1:3)]
  return(list(Var_avg.outcome, Var_contrast))
  
}

Y_IPCW1 <- function(data, net, alphas, id.var, trt.var, n.neighbor, n.trt.neighbor,
                    propensity.covars, censor.indicator, censor.covars,
                    censor.correlation = FALSE, component.var, obs.y) {
  p <- lapply(alphas, function(a)
    pi.neighbor(data, n.neighbor, n.trt.neighbor, alpha = a))
  propensity.score <- propensity.1(data = data, net = net, id.var = id.var, 
                                   trt.var = trt.var, propensity.covars = propensity.covars)
  censor.score <- censor(data = data, censor.indicator = censor.indicator,
                         censor.covars = censor.covars, censor.correlation = censor.correlation,
                         component.var = component.var)
  denominator = propensity.score * censor.score
  
  Y_IPW1_1 <- sapply(p, function(pp) {mean(data[[obs.y]] * data[[trt.var]] * pp / denominator)})
  Y_IPW1_0 <- sapply(p, function(pp) {mean(data[[obs.y]] * (1 - data[[trt.var]]) * pp / denominator)})
  Y_IPW1_margin <- mapply(function(pp, a) {mean(data[[obs.y]] * dbinom(data[[trt.var]], 1, a) * pp / denominator)}, p, alphas)
  avg.PO <- cbind(Y_IPW1_1, Y_IPW1_0, Y_IPW1_margin)
  colnames(avg.PO) <- c("treat = 1", "treat = 0", "margin")
  rownames(avg.PO) <- paste0("alpha = ", alphas)
  pairs <- combn(1:length(alphas), 2)
  DE <- data.frame()
  for (k in 1:length(alphas)){
    alpha1 <- alphas[k]
    alpha2 <- alphas[k]
    estimate <- avg.PO[k, "treat = 1"] - avg.PO[k, "treat = 0"]
    DE <- rbind(DE, c(alpha1, alpha2, estimate))
  }
  names(DE) <- c("alpha1", "alpha2", "estimate")
  IE <- data.frame()
  for (k in 1:ncol(pairs)){
    alpha1 <- alphas[pairs[1, k]]
    alpha2 <- alphas[pairs[2, k]]
    estimate <- avg.PO[pairs[2, k], "treat = 0"] - avg.PO[pairs[1, k], "treat = 0"]
    IE <- rbind(IE, c(alpha1, alpha2, estimate))
  }
  names(IE) <- c("alpha1", "alpha2", "estimate")
  TE <- data.frame()
  for (k in 1:ncol(pairs)){
    alpha1 <- alphas[pairs[1, k]]
    alpha2 <- alphas[pairs[2, k]]
    estimate <- avg.PO[pairs[2, k], "treat = 1"] - avg.PO[pairs[1, k], "treat = 0"]
    TE <- rbind(TE, c(alpha1, alpha2, estimate))
  }
  names(TE) <- c("alpha1", "alpha2", "estimate")
  DE$effect <- "DIRECT"
  IE$effect <- "INDIRECT"
  TE$effect <- "TOTAL"
  contrast <- rbind(DE, IE, TE)
  contrast <- contrast[, c(4, 1, 2, 3)]
  return(list(avg.PO, contrast))
}

Var_Y_IPCW1 <- function(alphas, net, data, id.var, obs.y, n.neighbor, n.trt.neighbor, 
                        trt.var, propensity.covars, censor.indicator, censor.covars, 
                        censor.correlation = FALSE, component.var){
  n <- nrow(data)
  m <- length(unique(data[[component.var]]))
  
  p <- lapply(alphas, function(a)
    pi.neighbor(data, n.neighbor, n.trt.neighbor, alpha = a))
  propensity.score <- propensity.1(data = data, net = net, id.var = id.var, 
                                   trt.var = trt.var, propensity.covars = propensity.covars)
  censor.score <- censor(data = data, censor.indicator = censor.indicator,
                         censor.covars = censor.covars, censor.correlation = censor.correlation,
                         component.var = component.var)
  nn.data <- data.frame()
  for (i in data[[id.var]]) {
    ids <- as.integer(neighborhood(net, order = 1, nodes = as.character(i))[[1]])
    temp.data <- data[data[[id.var]] %in% ids, ]
    temp.data$ego <- i
    nn.data <- rbind(nn.data, temp.data)
  }
  form <- as.formula(paste(trt.var, "~", paste(propensity.covars, collapse = " + "), "+ (1 | ego)"))
  propensity.model <- lme4::glmer(form, data = nn.data, family = binomial(link = "logit"))
  theta <- c(lme4::fixef(propensity.model), as.data.frame(lme4::VarCorr(propensity.model))$sdcor[1])
  
  if (censor.correlation) {
    form <- as.formula(paste(censor.indicator, "~",paste(censor.covars, collapse = " + "), 
                             "+ (1 |", component.var, ")"))
    censor.model <- lme4::glmer(form, data = data, family = binomial(link = "logit"))
    eta <- c(lme4::fixef(censor.model), as.data.frame(lme4::VarCorr(censor.model))$sdcor[1])
  } else {
    form <- as.formula(paste(censor.indicator, "~", paste(censor.covars, collapse = " + ")))
    censor.model <- glm(form, data = data, family = binomial(link = "logit"))
    eta <- coef(censor.model)
  }
  
  propen_U_11 <- matrix(0, nrow = length(theta), ncol = length(theta))
  censor_U_11 <- matrix(0, nrow = length(eta), ncol = length(eta))
  propen_v_11 <- matrix(NA, nrow = n, ncol = length(theta))
  censor_v_11 <- matrix(NA, nrow = n, ncol = length(eta))
  grad_propensity_inv <- matrix(NA, nrow = n, ncol = length(theta))
  grad_censor_inv <- matrix(NA, nrow = n, ncol = length(eta))
  for (i in data[[id.var]]){
    
    #####################################################################
    ## individual propensity and censoring functions given individual i##
    #####################################################################
    propensity_ind <- function(theta) {
      ids <- as.integer(neighborhood(net, order = 1, nodes = as.character(i))[[1]])
      temp.data <- data[data[[id.var]] %in% ids, ]
      X <- model.matrix(as.formula(paste("~", paste(propensity.covars, collapse = " + "))),
                        data = temp.data)
      beta <- theta[-length(theta)]
      sigma.b <- theta[length(theta)]
      X <- cbind(1, temp.data[, propensity.covars])
      integrand <- function(b) {
        prob = 1
        for (each.row in 1:nrow(temp.data)){
          prob.new <- plogis(as.vector(X[each.row, ] %*% as.matrix(theta[-length(theta)])) + b)
          if (temp.data[[trt.var]][each.row] == 1){
            prob <- prob*prob.new
          } else {prob <- prob*(1-prob.new)}
        }
        return(prob * dnorm(b, mean = 0, sd = theta[length(theta)]))
      }
      integrate(integrand, -Inf, Inf)$value
    }
    propensity_ind_inv <- function(theta) {1 / propensity_ind(theta)}
    log_propensity_ind <- function(theta){return(log(propensity_ind(theta)))}
    censor_ind <- function(eta) {
      subdata <- data[data[[id.var]] == i, , drop = FALSE]
      X <- model.matrix(as.formula(paste("~", paste(censor.covars, collapse = " + "))), data = subdata)
      if (censor.correlation) {
        beta <- eta[-length(eta)]
        sigma.b <- eta[length(eta)]
        prob <- plogis(as.vector(X %*% beta + sigma.b))
      } else {
        prob <- plogis(as.vector(X %*% eta))
      }
      prod(dbinom(subdata[[censor.indicator]], 1, prob))
    }
    censor_ind_inv <- function(eta){return(1/censor_ind(eta))}
    log_censor_ind <- function(eta){return(log(censor_ind(eta)))}
    
    #####################################################################
    ## Construct matrices ###############################################
    #####################################################################
    
    propen_U_11 <- propen_U_11-(1/n)*hessian(log_propensity_ind, x = theta)
    censor_U_11 <- censor_U_11-(1/n)*hessian(log_censor_ind, x = eta)
    propen_v_11[data[[id.var]] == i, ] <- grad(log_propensity_ind, x = theta)
    censor_v_11[data[[id.var]] == i, ] <- grad(log_censor_ind, x = eta)
    grad_propensity_inv[data[[id.var]] == i, ] <- grad(propensity_ind_inv, x = theta)
    grad_censor_inv[data[[id.var]] == i, ] <- grad(censor_ind_inv, x = eta)
  }
  U_11=rbind(cbind(propen_U_11, matrix(0, nrow = length(theta), ncol = length(eta))), 
             cbind(matrix(0, nrow = length(eta), ncol = length(theta)), censor_U_11))
  
  propen_V_11=matrix(0, nrow = length(theta), ncol = length(theta))
  censor_V_11=matrix(0, nrow = length(eta), ncol = length(eta))
  for (j in unique(data[[component.var]])){
    idx <- which(data[[component.var]] == j)
    propen_sum <- colSums(propen_v_11[idx, , drop = FALSE])
    censor_sum <- colSums(censor_v_11[idx, , drop = FALSE])
    
    propen_V_11 <- propen_V_11 + propen_sum %*% t(propen_sum)
    censor_V_11 <- censor_V_11 + censor_sum %*% t(censor_sum)
  }
  
  V_11=(m/n^2)*rbind(cbind(propen_V_11, matrix(0, nrow = length(theta), ncol = length(eta))),
                     cbind(matrix(0, nrow = length(eta), ncol = length(theta)), censor_V_11))
  
  Var_avg.outcome <- matrix(NA, ncol = 3, nrow = length(alphas))
  for (k in 1:length(alphas)){
    U_21.1 <- c(-colMeans((data[[obs.y]]*data[[trt.var]]*p[[k]])*grad_propensity_inv), 
                -colMeans((data[[obs.y]]*data[[trt.var]]*p[[k]])*grad_censor_inv))
    B.1 <- data[[obs.y]]*data[[trt.var]]*p[[k]]/(propensity.score*censor.score)
    Y.1 <- mean(B.1)
    
    U_21.0 <- c(-colMeans((data[[obs.y]]*(1-data[[trt.var]])*p[[k]])*grad_propensity_inv), 
                -colMeans((data[[obs.y]]*(1-data[[trt.var]])*p[[k]])*grad_censor_inv))
    B.0 <- data[[obs.y]]*(1-data[[trt.var]])*p[[k]]/(propensity.score*censor.score)
    Y.0 <- mean(B.0)
    
    B <- data[[obs.y]]*dbinom(data[[trt.var]], 1, alphas[[k]])*p[[k]]/(propensity.score*censor.score)
    Y <- mean(B)
    
    v_22.1 <- 0
    v_21.1 <- rep(0, length(theta)+length(eta))
    v_22.0 <- 0
    v_21.0 <- rep(0, length(theta)+length(eta))
    v_22 <- 0
    v_21 <- rep(0, length(theta)+length(eta))
    for (j in unique(data[[component.var]])) {
      idx <- which(data[[component.var]] == j)
      weight <- m / n
      r.1 <- weight * sum(B.1[idx]) - Y.1
      r.0 <- weight * sum(B.0[idx]) - Y.0
      r   <- weight * sum(B[idx])   - Y
      v_21.1 <- v_21.1 + c(
        weight * colSums(propen_v_11[idx, , drop = FALSE]) * r.1,
        weight * colSums(censor_v_11[idx, , drop = FALSE]) * r.1
      )
      v_22.1 <- v_22.1 + r.1^2
      
      v_21.0 <- v_21.0 + c(
        weight * colSums(propen_v_11[idx, , drop = FALSE]) * r.0,
        weight * colSums(censor_v_11[idx, , drop = FALSE]) * r.0
      )
      v_22.0 <- v_22.0 + r.0^2
      
      v_21 <- v_21 + c(
        weight * colSums(propen_v_11[idx, , drop = FALSE]) * r,
        weight * colSums(censor_v_11[idx, , drop = FALSE]) * r
      )
      v_22 <- v_22 + r^2
    }
    
    U=cbind(U_11, rep(0, length(theta)+length(eta)))
    U.1=rbind(U, c(U_21.1, 1))
    U.0=rbind(U, c(U_21.0, 1))
    V.1=cbind(V_11, v_21.1/m)
    V.1=rbind(V.1, c(v_21.1, v_22.1)/m)
    V.0=cbind(V_11, v_21.0/m)
    V.0=rbind(V.0, c(v_21.0, v_22.0)/m)
    U_inv.1=solve(U.1)
    U_inv.0=solve(U.0)
    M.1=U_inv.1%*%V.1%*%t(U_inv.1)
    M.0=U_inv.0%*%V.0%*%t(U_inv.0)
    total.par = length(theta)+length(eta)+1
    
    p.mar <- dbinom(data[[trt.var]], 1, alphas[k])* p[[k]]
    U <- cbind(U_11, rep(0, length(theta)+length(eta)))
    U <- rbind(U, c(-colMeans((data[[obs.y]]*p.mar)*grad_propensity_inv), 
                    -colMeans((data[[obs.y]]*p.mar)*grad_censor_inv), 1))
    V=cbind(V_11, v_21/m)
    V=rbind(V, c(v_21, v_22)/m)
    U_inv=solve(U)
    M=U_inv%*%V%*%t(U_inv)
    Var_avg.outcome[k, ] = c(M.1[total.par, total.par] / m, M.0[total.par, total.par] / m, 
                             M[total.par, total.par] / m)
  }
  colnames(Var_avg.outcome) <- c("treat = 1", "treat = 0", "margin")
  rownames(Var_avg.outcome) <- paste("alpha = ", alphas, sep = "")
  
  
  Var_DE <- data.frame()
  for (k in 1:length(alphas)){
    U_21 <- c(-colMeans((data[[obs.y]]*data[[trt.var]]*p[[k]])*grad_propensity_inv), 
              -colMeans((data[[obs.y]]*(1-data[[trt.var]])*p[[k]])*grad_censor_inv))
    B <- data[[obs.y]]*data[[trt.var]]*p[[k]]/(propensity.score*censor.score) - 
      data[[obs.y]]*(1-data[[trt.var]])*p[[k]]/(propensity.score*censor.score)
    Y <- mean(B)
    
    v_22 <- 0
    v_21 <- rep(0, length(theta)+length(eta))
    for (j in unique(data[[component.var]])) {
      idx <- which(data[[component.var]] == j)
      weight <- m / n
      r   <- weight * sum(B[idx])   - Y
      v_21 <- v_21 + c(
        weight * colSums(propen_v_11[idx, , drop = FALSE]) * r,
        weight * colSums(censor_v_11[idx, , drop = FALSE]) * r
      )
      v_22 <- v_22 + r^2
    }
    
    U <- cbind(U_11, rep(0, length(theta)+length(eta)))
    U <- rbind(U, c(U_21, 1))
    V <- cbind(V_11, v_21/m)
    V <- rbind(V, c(v_21, v_22)/m)
    U_inv <- solve(U)
    M <- U_inv%*%V%*%t(U_inv)
    total.par <- length(theta)+length(eta)+1
    Var_DE <- rbind(Var_DE, c(alphas[k], alphas[k], M[total.par, total.par] / m))
  }
  pairs <- combn(1:length(alphas), 2)
  
  Var_IE <- data.frame()
  for (k in 1:length(alphas)){
    U_21 <- c(-colMeans((data[[obs.y]]*(1-data[[trt.var]])*p[[pairs[2, k]]])*grad_propensity_inv), 
              -colMeans((data[[obs.y]]*(1-data[[trt.var]])*p[[pairs[1, k]]])*grad_censor_inv))
    B <- data[[obs.y]]*(1-data[[trt.var]])*p[[pairs[2, k]]]/(propensity.score*censor.score) -
      data[[obs.y]]*(1-data[[trt.var]])*p[[pairs[1, k]]]/(propensity.score*censor.score)
    Y <- mean(B)
    
    v_22 <- 0
    v_21 <- rep(0, length(theta)+length(eta))
    for (j in unique(data[[component.var]])) {
      idx <- which(data[[component.var]] == j)
      weight <- m / n
      r   <- weight * sum(B[idx])   - Y
      v_21 <- v_21 + c(
        weight * colSums(propen_v_11[idx, , drop = FALSE]) * r,
        weight * colSums(censor_v_11[idx, , drop = FALSE]) * r
      )
      v_22 <- v_22 + r^2
    }
    
    U <- cbind(U_11, rep(0, length(theta)+length(eta)))
    U <- rbind(U, c(U_21, 1))
    V <- cbind(V_11, v_21/m)
    V <- rbind(V, c(v_21, v_22)/m)
    U_inv <- solve(U)
    M <- U_inv%*%V%*%t(U_inv)
    total.par <- length(theta)+length(eta)+1
    Var_IE <- rbind(Var_IE, c(alphas[pairs[1, k]], alphas[pairs[2, k]], M[total.par, total.par] / m))
  }
  
  Var_TE <- data.frame()
  for (k in 1:length(alphas)){
    U_21 <- c(-colMeans((data[[obs.y]]*(data[[trt.var]])*p[[pairs[2, k]]])*grad_propensity_inv), 
              -colMeans((data[[obs.y]]*(1-data[[trt.var]])*p[[pairs[1, k]]])*grad_censor_inv))
    B <- data[[obs.y]]*(data[[trt.var]])*p[[pairs[2, k]]]/(propensity.score*censor.score) -
      data[[obs.y]]*(1-data[[trt.var]])*p[[pairs[1, k]]]/(propensity.score*censor.score)
    Y <- mean(B)
    
    v_22 <- 0
    v_21 <- rep(0, length(theta)+length(eta))
    for (j in unique(data[[component.var]])) {
      idx <- which(data[[component.var]] == j)
      weight <- m / n
      r   <- weight * sum(B[idx])   - Y
      v_21 <- v_21 + c(
        weight * colSums(propen_v_11[idx, , drop = FALSE]) * r,
        weight * colSums(censor_v_11[idx, , drop = FALSE]) * r
      )
      v_22 <- v_22 + r^2
    }
    
    U <- cbind(U_11, rep(0, length(theta)+length(eta)))
    U <- rbind(U, c(U_21, 1))
    V <- cbind(V_11, v_21/m)
    V <- rbind(V, c(v_21, v_22)/m)
    U_inv <- solve(U)
    M <- U_inv%*%V%*%t(U_inv)
    total.par <- length(theta)+length(eta)+1
    Var_TE <- rbind(Var_TE, c(alphas[pairs[1, k]], alphas[pairs[2, k]], M[total.par, total.par] / m))
  }
  
  names(Var_DE) <- c("alpha1", "alpha2", "variance")
  names(Var_IE) <- c("alpha1", "alpha2", "variance")
  names(Var_TE) <- c("alpha1", "alpha2", "variance")
  Var_DE$effect <- "DIRECT"
  Var_IE$effect <- "INDIRECT"
  Var_TE$effect <- "TOTAL"
  
  Var_contrast <- rbind(Var_DE, Var_IE, Var_TE)
  Var_contrast <- Var_contrast[, c(4, 1:3)]
  return(list(Var_avg.outcome, Var_contrast))
  
}

##################################################################################
################## IPW2 functions ################################################
##################################################################################

Y_IPW2 <- function(data, net, alphas, id.var, trt.var, n.neighbor, n.trt.neighbor,
                   propensity.covars, propensity.nn.covars, component.var = NULL, obs.y) {
  p <- lapply(alphas, function(a)
    pi.neighbor(data, n.neighbor, n.trt.neighbor, alpha = a))
  propensity.score <- propensity.2(data = data, trt.var = trt.var, n.trt.neighbor = n.trt.neighbor, 
                                   n.neighbor = n.neighbor, propensity.covars = propensity.covars, 
                                   propensity.nn.covars = propensity.nn.covars)
  
  Y_IPW1_1 <- sapply(p, function(pp) {mean(data[[obs.y]] * data[[trt.var]] * pp / propensity.score)})
  Y_IPW1_0 <- sapply(p, function(pp) {mean(data[[obs.y]] * (1 - data[[trt.var]]) * pp / propensity.score)})
  Y_IPW1_margin <- mapply(function(pp, a) {mean(data[[obs.y]] * dbinom(data[[trt.var]], 1, a) * pp / propensity.score)}, p, alphas)
  avg.outcome <- cbind(Y_IPW1_1, Y_IPW1_0, Y_IPW1_margin)
  colnames(avg.outcome) <- c("treat = 1", "treat = 0", "margin")
  rownames(avg.outcome) <- paste0("alpha = ", alphas)
  pairs <- combn(1:length(alphas), 2)
  DE <- data.frame()
  for (k in 1:length(alphas)){
    alpha1 <- alphas[k]
    alpha2 <- alphas[k]
    estimate <- avg.outcome[k, "treat = 1"] - avg.outcome[k, "treat = 0"]
    DE <- rbind(DE, c(alpha1, alpha2, estimate))
  }
  names(DE) <- c("alpha1", "alpha2", "estimate")
  IE <- data.frame()
  for (k in 1:ncol(pairs)){
    alpha1 <- alphas[pairs[1, k]]
    alpha2 <- alphas[pairs[2, k]]
    estimate <- avg.outcome[pairs[2, k], "treat = 0"] - avg.outcome[pairs[1, k], "treat = 0"]
    IE <- rbind(IE, c(alpha1, alpha2, estimate))
  }
  names(IE) <- c("alpha1", "alpha2", "estimate")
  TE <- data.frame()
  for (k in 1:ncol(pairs)){
    alpha1 <- alphas[pairs[1, k]]
    alpha2 <- alphas[pairs[2, k]]
    estimate <- avg.outcome[pairs[2, k], "treat = 1"] - avg.outcome[pairs[1, k], "treat = 0"]
    TE <- rbind(TE, c(alpha1, alpha2, estimate))
  }
  names(TE) <- c("alpha1", "alpha2", "estimate")
  DE$effect <- "DIRECT"
  IE$effect <- "INDIRECT"
  TE$effect <- "TOTAL"
  contrast <- rbind(DE, IE, TE)
  contrast <- contrast[, c(4, 1, 2, 3)]
  return(list(avg.outcome, contrast))
}

Var_Y_IPW2 <- function(alphas, net, data, id.var, obs.y, n.neighbor, n.trt.neighbor, 
                       trt.var, propensity.covars, propensity.nn.covars,
                       component.var){
  n <- nrow(data)
  m <- length(unique(data[[component.var]]))
  
  p <- lapply(alphas, function(a)
    pi.neighbor(data, n.neighbor, n.trt.neighbor, alpha = a))
  propensity.score <- propensity.2(data = data, trt.var = trt.var, n.trt.neighbor = n.trt.neighbor, 
                                   n.neighbor = n.neighbor, propensity.covars = propensity.covars, 
                                   propensity.nn.covars = propensity.nn.covars)
  
  data$not.trt.neighbor <- data[[n.neighbor]] - data[[n.trt.neighbor]] 
  form.ind <- as.formula(paste(trt.var, "~", paste(propensity.covars, collapse = " + ")))
  form.nn <- as.formula(paste("cbind(", n.trt.neighbor, ", not.trt.neighbor)~", 
                              trt.var, "+", paste(propensity.covars, collapse = " + "), 
                              "+", paste(propensity.nn.covars, collapse = " + ")))
  propensity.model.ind <- glm(form.ind, data = data, family = binomial(link = "logit"))
  propensity.model.nn <- glm(form.nn, data = data, family = binomial(link = "logit"))
  
  theta.ind <- coef(propensity.model.ind)
  theta.nn <- coef(propensity.model.nn)
  theta <- c(theta.ind, theta.nn)
  
  propen_U_11 <- matrix(0, nrow = length(theta), ncol = length(theta))
  propen_v_11 <- matrix(NA, nrow = n, ncol = length(theta))
  grad_propensity_inv <- matrix(NA, nrow = n, ncol = length(theta))
  for (i in data[[id.var]]){
    
    #####################################################################
    ## individual propensity and censoring functions given individual i##
    #####################################################################
    propensity.2_ind <- function(theta) {
      temp.data <- data[data[[id.var]] == i, ]
      theta1 <- theta[1:length(theta.ind)]
      theta2 <- theta[(length(theta.ind)+1):length(theta)]
      X1 <- model.matrix(as.formula(paste("~", paste(propensity.covars, collapse = " + "))),
                         data = temp.data)
      X2 <- model.matrix(as.formula(paste("~", trt.var, "+",
                                          paste(propensity.covars, collapse = " + "), "+",
                                          paste(propensity.nn.covars, collapse = " + "))), 
                         data = temp.data)
      p.ind <- as.numeric(plogis(X1 %*% theta1))
      p.nn <- as.numeric(plogis(X2 %*% theta2))
      
      return(dbinom(temp.data[[trt.var]], size = 1, prob = p.ind)*(p.nn^temp.data[[n.trt.neighbor]])*((1-p.nn)^(temp.data$not.trt.neighbor)))
    }
    propensity.2_ind_inv <- function(theta) {1 / propensity.2_ind(theta)}
    log_propensity.2_ind <- function(theta){return(log(propensity.2_ind(theta)))}
    
    #####################################################################
    ## Construct matrices ###############################################
    #####################################################################
    
    propen_U_11 <- propen_U_11-(1/n)*hessian(log_propensity.2_ind, x = theta)
    propen_v_11[data[[id.var]] == i, ] <- grad(log_propensity.2_ind, x = theta)
    grad_propensity_inv[data[[id.var]] == i, ] <- grad(propensity.2_ind_inv, x = theta)
  }
  U_11 = propen_U_11
  
  propen_V_11 <- matrix(0, nrow = length(theta), ncol = length(theta))
  for (j in unique(data[[component.var]])){
    idx <- which(data[[component.var]] == j)
    propen_sum <- colSums(propen_v_11[idx, , drop = FALSE])
    propen_V_11 <- propen_V_11 + propen_sum %*% t(propen_sum)
  }
  
  V_11=(m/n^2)*propen_V_11
  
  Var_avg.outcome <- matrix(NA, nrow = length(alphas), ncol = 3)
  for (k in 1:length(alphas)){
    total.par <- length(theta)+1
    U_21.1 <- -colMeans((data[[obs.y]] * data[[trt.var]] * p[[k]]) * grad_propensity_inv)
    B.1 <- data[[obs.y]] * data[[trt.var]] * p[[k]] / propensity.score
    Y.1 <- mean(B.1)
    U_21.0 <- -colMeans((data[[obs.y]]*(1-data[[trt.var]])*p[[k]])*grad_propensity_inv)
    B.0 <- data[[obs.y]]*(1-data[[trt.var]])*p[[k]]/propensity.score
    Y.0 <- mean(B.0)
    p.mar <- dbinom(data[[trt.var]], 1, alphas[[k]])* p[[k]]
    U_21 <- -colMeans((data[[obs.y]]*p.mar)*grad_propensity_inv)
    B <- data[[obs.y]] * dbinom(data[[trt.var]], 1, alphas[[k]]) * p[[k]] / propensity.score
    Y <- mean(B)
    
    v_22.1 <- 0
    v_21.1 <- rep(0, length(theta))
    v_22.0 <- 0
    v_21.0 <- rep(0, length(theta))
    v_22 <- 0
    v_21 <- rep(0, length(theta))
    for (j in unique(data[[component.var]])) {
      idx <- which(data[[component.var]] == j)
      weight <- m / n
      r.1 <- weight * sum(B.1[idx]) - Y.1
      r.0 <- weight * sum(B.0[idx]) - Y.0
      r   <- weight * sum(B[idx])   - Y
      v_21.1 <- v_21.1 + weight * colSums(propen_v_11[idx, , drop = FALSE]) * r.1
      v_22.1 <- v_22.1 + r.1^2
      
      v_21.0 <- v_21.0 + weight * colSums(propen_v_11[idx, , drop = FALSE]) * r.0
      v_22.0 <- v_22.0 + r.0^2
      
      v_21 <- v_21 +  weight * colSums(propen_v_11[idx, , drop = FALSE]) * r
      v_22 <- v_22 + r^2
    }
    U <- cbind(U_11, rep(0, length(theta)))
    U.1 <- rbind(U, c(U_21.1, 1))
    U.0 <- rbind(U, c(U_21.0, 1))
    V.1 <- cbind(V_11, v_21.1/m)
    V.1 <- rbind(V.1, c(v_21.1, v_22.1)/m)
    V.0 <- cbind(V_11, v_21.0/m)
    V.0 <- rbind(V.0, c(v_21.0, v_22.0)/m)
    U_inv.1 <- solve(U.1)
    U_inv.0 <- solve(U.0)
    M.1 <- U_inv.1%*%V.1%*%t(U_inv.1)
    M.0 <- U_inv.0%*%V.0%*%t(U_inv.0)
    U <- rbind(U, c(U_21, 1))
    V <- cbind(V_11, v_21/m)
    V <- rbind(V, c(v_21, v_22)/m)
    U_inv <- solve(U)
    M <- U_inv%*%V%*%t(U_inv)
    var <- c(M.1[total.par, total.par] / m, M.0[total.par, total.par] / m, M[total.par, total.par] / m)
    Var_avg.outcome[k, ] <- var
  }
  colnames(Var_avg.outcome) <- c("treat = 1", "treat = 0", "margin")
  rownames(Var_avg.outcome) <- paste0("alpha = ", alphas)
  
  
  Var_DE <- data.frame()
  for (k in 1:length(alphas)){
    U_21 <- -colMeans((data[[obs.y]]*(data[[trt.var]])*p[[k]])*grad_propensity_inv - 
                        (data[[obs.y]]*(1-data[[trt.var]])*p[[k]])*grad_propensity_inv)
    B <- data[[obs.y]]*data[[trt.var]]*p[[k]]/propensity.score - 
      data[[obs.y]]*(1-data[[trt.var]])*p[[k]]/propensity.score
    Y <- mean(B)
    
    v_22 <- 0
    v_21 <- rep(0, length(theta))
    for (j in unique(data[[component.var]])) {
      idx <- which(data[[component.var]] == j)
      weight <- m / n
      r   <- weight * sum(B[idx]) - Y
      v_21 <- v_21 + weight * colSums(propen_v_11[idx, , drop = FALSE]) * r
      v_22 <- v_22 + r^2
    }
    U <- cbind(U_11, rep(0, length(theta)))
    U <- rbind(U, c(U_21, 1))
    V <- cbind(V_11, v_21/m)
    V <- rbind(V, c(v_21, v_22)/m)
    U_inv <- solve(U)
    M <- U_inv%*%V%*%t(U_inv)
    total.par <- length(theta)+1
    Var_DE <- rbind(Var_DE, c(alphas[k], alphas[k], M[total.par, total.par] / m))
  }
  names(Var_DE) <- c("alpha1", "alpha2", "Variance")
  
  
  pairs <- combn(1:length(alphas), 2)
  
  Var_IE <- data.frame()
  for (k in 1:length(alphas)){
    U_21 <- -colMeans((data[[obs.y]]*(1-data[[trt.var]])*p[[pairs[2, k]]])*grad_propensity_inv - 
                        (data[[obs.y]]*(1-data[[trt.var]])*p[[pairs[1, k]]])*grad_propensity_inv)
    B <- data[[obs.y]]*(1-data[[trt.var]])*p[[pairs[2, k]]]/propensity.score - 
      data[[obs.y]]*(1-data[[trt.var]])*p[[pairs[1, k]]]/propensity.score
    Y <- mean(B)
    
    v_22 <- 0
    v_21 <- rep(0, length(theta))
    for (j in unique(data[[component.var]])) {
      idx <- which(data[[component.var]] == j)
      weight <- m / n
      r   <- weight * sum(B[idx]) - Y
      v_21 <- v_21 + weight * colSums(propen_v_11[idx, , drop = FALSE]) * r
      v_22 <- v_22 + r^2
    }
    U <- cbind(U_11, rep(0, length(theta)))
    U <- rbind(U, c(U_21, 1))
    V <- cbind(V_11, v_21/m)
    V <- rbind(V, c(v_21, v_22)/m)
    U_inv <- solve(U)
    M <- U_inv%*%V%*%t(U_inv)
    total.par <- length(theta)+1
    Var_IE <- rbind(Var_IE, c(alphas[pairs[1, k]], alphas[pairs[2, k]], M[total.par, total.par] / m))
  }
  names(Var_IE) <- c("alpha1", "alpha2", "Variance")
  
  Var_TE <- data.frame()
  for (k in 1:length(alphas)){
    U_21 <- -colMeans((data[[obs.y]]*(data[[trt.var]])*p[[pairs[2, k]]])*grad_propensity_inv - 
                        (data[[obs.y]]*(1-data[[trt.var]])*p[[pairs[1, k]]])*grad_propensity_inv)
    B <- data[[obs.y]]*(data[[trt.var]])*p[[pairs[2, k]]]/propensity.score - 
      data[[obs.y]]*(1-data[[trt.var]])*p[[pairs[1, k]]]/propensity.score
    Y <- mean(B)
    
    v_22 <- 0
    v_21 <- rep(0, length(theta))
    for (j in unique(data[[component.var]])) {
      idx <- which(data[[component.var]] == j)
      weight <- m / n
      r   <- weight * sum(B[idx]) - Y
      v_21 <- v_21 + weight * colSums(propen_v_11[idx, , drop = FALSE]) * r
      v_22 <- v_22 + r^2
    }
    U <- cbind(U_11, rep(0, length(theta)))
    U <- rbind(U, c(U_21, 1))
    V <- cbind(V_11, v_21/m)
    V <- rbind(V, c(v_21, v_22)/m)
    U_inv <- solve(U)
    M <- U_inv%*%V%*%t(U_inv)
    total.par <- length(theta)+1
    Var_TE <- rbind(Var_TE, c(alphas[pairs[1, k]], alphas[pairs[2, k]], M[total.par, total.par] / m))
  }
  names(Var_TE) <- c("alpha1", "alpha2", "Variance")
  
  Var_DE$effect <- "DIRECT"
  Var_IE$effect <- "INDIRECT"
  Var_TE$effect <- "TOTAL"
  
  Var_contrast <- rbind(Var_DE, Var_IE, Var_TE)
  Var_contrast <- Var_contrast[, c(4, 1:3)]
  return(list(Var_avg.outcome, Var_contrast))
  
}

Y_IPCW2 <- function(data, net, alphas, id.var, trt.var, n.neighbor, n.trt.neighbor,
                    propensity.covars, propensity.nn.covars, censor.indicator, censor.covars,
                    censor.correlation = FALSE, component.var = NULL, obs.y) {
  p <- lapply(alphas, function(a)
    pi.neighbor(data, n.neighbor, n.trt.neighbor, alpha = a))
  propensity.score <- propensity.2(data = data, trt.var = trt.var, n.trt.neighbor = n.trt.neighbor, 
                                   n.neighbor = n.neighbor, propensity.covars = propensity.covars, 
                                   propensity.nn.covars = propensity.nn.covars)
  censor.score <- censor(data = data, censor.indicator = censor.indicator,
                         censor.covars = censor.covars, censor.correlation = censor.correlation,
                         component.var = component.var)
  denominator = propensity.score * censor.score
  Y_IPW1_1 <- sapply(p, function(pp) {mean(data[[obs.y]] * data[[trt.var]] * pp / denominator)})
  Y_IPW1_0 <- sapply(p, function(pp) {mean(data[[obs.y]] * (1 - data[[trt.var]]) * pp / denominator)})
  Y_IPW1_margin <- mapply(function(pp, a) {mean(data[[obs.y]] * dbinom(data[[trt.var]], 1, a) * pp / denominator)}, p, alphas)
  avg.PO <- cbind(Y_IPW1_1, Y_IPW1_0, Y_IPW1_margin)
  colnames(avg.PO) <- c("treat = 1", "treat = 0", "margin")
  rownames(avg.PO) <- paste0("alpha = ", alphas)
  pairs <- combn(1:length(alphas), 2)
  DE <- data.frame()
  for (k in 1:length(alphas)){
    alpha1 <- alphas[k]
    alpha2 <- alphas[k]
    estimate <- avg.PO[k, "treat = 1"] - avg.PO[k, "treat = 0"]
    DE <- rbind(DE, c(alpha1, alpha2, estimate))
  }
  names(DE) <- c("alpha1", "alpha2", "estimate")
  IE <- data.frame()
  for (k in 1:ncol(pairs)){
    alpha1 <- alphas[pairs[1, k]]
    alpha2 <- alphas[pairs[2, k]]
    estimate <- avg.PO[pairs[2, k], "treat = 0"] - avg.PO[pairs[1, k], "treat = 0"]
    IE <- rbind(IE, c(alpha1, alpha2, estimate))
  }
  names(IE) <- c("alpha1", "alpha2", "estimate")
  TE <- data.frame()
  for (k in 1:ncol(pairs)){
    alpha1 <- alphas[pairs[1, k]]
    alpha2 <- alphas[pairs[2, k]]
    estimate <- avg.PO[pairs[2, k], "treat = 1"] - avg.PO[pairs[1, k], "treat = 0"]
    TE <- rbind(TE, c(alpha1, alpha2, estimate))
  }
  names(TE) <- c("alpha1", "alpha2", "estimate")
  DE$effect <- "DIRECT"
  IE$effect <- "INDIRECT"
  TE$effect <- "TOTAL"
  contrast <- rbind(DE, IE, TE)
  contrast <- contrast[, c(4, 1, 2, 3)]
  return(list(avg.PO, contrast))
}

Var_Y_IPCW2 <- function(alphas, net, data, id.var, obs.y, n.neighbor, n.trt.neighbor, 
                        trt.var, propensity.covars, propensity.nn.covars,
                        censor.indicator, censor.covars, 
                        censor.correlation = FALSE, component.var){
  n <- nrow(data)
  m <- length(unique(data[[component.var]]))
  
  p <- lapply(alphas, function(a)
    pi.neighbor(data, n.neighbor, n.trt.neighbor, alpha = a))
  propensity.score <- propensity.2(data = data, trt.var = trt.var, n.trt.neighbor = n.trt.neighbor, 
                                   n.neighbor = n.neighbor, propensity.covars = propensity.covars, 
                                   propensity.nn.covars = propensity.nn.covars)
  censor.score <- censor(data = data, censor.indicator = censor.indicator,
                         censor.covars = censor.covars, censor.correlation = censor.correlation,
                         component.var = component.var)
  
  data$not.trt.neighbor <- data[[n.neighbor]] - data[[n.trt.neighbor]] 
  form.ind <- as.formula(paste(trt.var, "~", paste(propensity.covars, collapse = " + ")))
  form.nn <- as.formula(paste("cbind(", n.trt.neighbor, ", not.trt.neighbor)~", 
                              trt.var, "+", paste(propensity.covars, collapse = " + "), 
                              "+", paste(propensity.nn.covars, collapse = " + ")))
  propensity.model.ind <- glm(form.ind, data = data, family = binomial(link = "logit"))
  propensity.model.nn <- glm(form.nn, data = data, family = binomial(link = "logit"))
  
  theta.ind <- coef(propensity.model.ind)
  theta.nn <- coef(propensity.model.nn)
  theta <- c(theta.ind, theta.nn)
  
  if (censor.correlation) {
    form <- as.formula(paste(censor.indicator, "~",paste(censor.covars, collapse = " + "), 
                             "+ (1 |", component.var, ")"))
    censor.model <- lme4::glmer(form, data = data, family = binomial(link = "logit"))
    eta <- c(lme4::fixef(censor.model), as.data.frame(lme4::VarCorr(censor.model))$sdcor[1])
  } else {
    form <- as.formula(paste(censor.indicator, "~", paste(censor.covars, collapse = " + ")))
    censor.model <- glm(form, data = data, family = binomial(link = "logit"))
    eta <- coef(censor.model)
  }
  
  propen_U_11 <- matrix(0, nrow = length(theta), ncol = length(theta))
  censor_U_11 <- matrix(0, nrow = length(eta), ncol = length(eta))
  propen_v_11 <- matrix(NA, nrow = n, ncol = length(theta))
  censor_v_11 <- matrix(NA, nrow = n, ncol = length(eta))
  grad_propensity_inv <- matrix(NA, nrow = n, ncol = length(theta))
  grad_censor_inv <- matrix(NA, nrow = n, ncol = length(eta))
  for (i in data[[id.var]]){
    
    #####################################################################
    ## individual propensity and censoring functions given individual i##
    #####################################################################
    propensity.2_ind <- function(theta) {
      temp.data <- data[data[[id.var]] == i, ]
      theta1 <- theta[1:length(theta.ind)]
      theta2 <- theta[(length(theta.ind)+1):length(theta)]
      X1 <- model.matrix(as.formula(paste("~", paste(propensity.covars, collapse = " + "))),
                         data = temp.data)
      X2 <- model.matrix(as.formula(paste("~", trt.var, "+",
                                          paste(propensity.covars, collapse = " + "), "+",
                                          paste(propensity.nn.covars, collapse = " + "))), 
                         data = temp.data)
      p.ind <- as.numeric(plogis(X1 %*% theta1))
      p.nn <- as.numeric(plogis(X2 %*% theta2))
      
      return(dbinom(temp.data[[trt.var]], size = 1, prob = p.ind)*(p.nn^temp.data[[n.trt.neighbor]])*((1-p.nn)^(temp.data$not.trt.neighbor)))
    }
    propensity.2_ind_inv <- function(theta) {1 / propensity.2_ind(theta)}
    log_propensity.2_ind <- function(theta){return(log(propensity.2_ind(theta)))}
    censor_ind <- function(eta) {
      subdata <- data[data[[id.var]] == i, , drop = FALSE]
      X <- model.matrix(as.formula(paste("~", paste(censor.covars, collapse = " + "))), data = subdata)
      if (censor.correlation) {
        beta <- eta[-length(eta)]
        sigma.b <- eta[length(eta)]
        prob <- plogis(as.vector(X %*% beta + sigma.b))
      } else {
        prob <- plogis(as.vector(X %*% eta))
      }
      prod(dbinom(subdata[[censor.indicator]], 1, prob))
    }
    censor_ind_inv <- function(eta){return(1/censor_ind(eta))}
    log_censor_ind <- function(eta){return(log(censor_ind(eta)))}
    
    #####################################################################
    ## Construct matrices ###############################################
    #####################################################################
    
    propen_U_11 <- propen_U_11-(1/n)*hessian(log_propensity.2_ind, x = theta)
    censor_U_11 <- censor_U_11-(1/n)*hessian(log_censor_ind, x = eta)
    propen_v_11[data[[id.var]] == i, ] <- grad(log_propensity.2_ind, x = theta)
    censor_v_11[data[[id.var]] == i, ] <- grad(log_censor_ind, x = eta)
    grad_propensity_inv[data[[id.var]] == i, ] <- grad(propensity.2_ind_inv, x = theta)
    grad_censor_inv[data[[id.var]] == i, ] <- grad(censor_ind_inv, x = eta)
  }
  U_11=rbind(cbind(propen_U_11, matrix(0, nrow = length(theta), ncol = length(eta))), 
             cbind(matrix(0, nrow = length(eta), ncol = length(theta)), censor_U_11))
  
  propen_V_11=matrix(0, nrow = length(theta), ncol = length(theta))
  censor_V_11=matrix(0, nrow = length(eta), ncol = length(eta))
  for (j in unique(data[[component.var]])){
    idx <- which(data[[component.var]] == j)
    propen_sum <- colSums(propen_v_11[idx, , drop = FALSE])
    censor_sum <- colSums(censor_v_11[idx, , drop = FALSE])
    
    propen_V_11 <- propen_V_11 + propen_sum %*% t(propen_sum)
    censor_V_11 <- censor_V_11 + censor_sum %*% t(censor_sum)
  }
  
  V_11=(m/n^2)*rbind(cbind(propen_V_11, matrix(0, nrow = length(theta), ncol = length(eta))),
                     cbind(matrix(0, nrow = length(eta), ncol = length(theta)), censor_V_11))
  
  Var_avg.outcome <- matrix(NA, ncol = 3, nrow = length(alphas))
  for (k in 1:length(alphas)){
    U_21.1 <- c(-colMeans((data[[obs.y]]*data[[trt.var]]*p[[k]])*grad_propensity_inv), 
                -colMeans((data[[obs.y]]*data[[trt.var]]*p[[k]])*grad_censor_inv))
    B.1 <- data[[obs.y]]*data[[trt.var]]*p[[k]]/(propensity.score*censor.score)
    Y.1 <- mean(B.1)
    
    U_21.0 <- c(-colMeans((data[[obs.y]]*(1-data[[trt.var]])*p[[k]])*grad_propensity_inv), 
                -colMeans((data[[obs.y]]*(1-data[[trt.var]])*p[[k]])*grad_censor_inv))
    B.0 <- data[[obs.y]]*(1-data[[trt.var]])*p[[k]]/(propensity.score*censor.score)
    Y.0 <- mean(B.0)
    
    B <- data[[obs.y]]*dbinom(data[[trt.var]], 1, alphas[[k]])*p[[k]]/(propensity.score*censor.score)
    Y <- mean(B)
    
    v_22.1 <- 0
    v_21.1 <- rep(0, length(theta)+length(eta))
    v_22.0 <- 0
    v_21.0 <- rep(0, length(theta)+length(eta))
    v_22 <- 0
    v_21 <- rep(0, length(theta)+length(eta))
    for (j in unique(data[[component.var]])) {
      idx <- which(data[[component.var]] == j)
      weight <- m / n
      r.1 <- weight * sum(B.1[idx]) - Y.1
      r.0 <- weight * sum(B.0[idx]) - Y.0
      r   <- weight * sum(B[idx])   - Y
      v_21.1 <- v_21.1 + c(
        weight * colSums(propen_v_11[idx, , drop = FALSE]) * r.1,
        weight * colSums(censor_v_11[idx, , drop = FALSE]) * r.1
      )
      v_22.1 <- v_22.1 + r.1^2
      
      v_21.0 <- v_21.0 + c(
        weight * colSums(propen_v_11[idx, , drop = FALSE]) * r.0,
        weight * colSums(censor_v_11[idx, , drop = FALSE]) * r.0
      )
      v_22.0 <- v_22.0 + r.0^2
      
      v_21 <- v_21 + c(
        weight * colSums(propen_v_11[idx, , drop = FALSE]) * r,
        weight * colSums(censor_v_11[idx, , drop = FALSE]) * r
      )
      v_22 <- v_22 + r^2
    }
    
    U=cbind(U_11, rep(0, length(theta)+length(eta)))
    U.1=rbind(U, c(U_21.1, 1))
    U.0=rbind(U, c(U_21.0, 1))
    V.1=cbind(V_11, v_21.1/m)
    V.1=rbind(V.1, c(v_21.1, v_22.1)/m)
    V.0=cbind(V_11, v_21.0/m)
    V.0=rbind(V.0, c(v_21.0, v_22.0)/m)
    U_inv.1=solve(U.1)
    U_inv.0=solve(U.0)
    M.1=U_inv.1%*%V.1%*%t(U_inv.1)
    M.0=U_inv.0%*%V.0%*%t(U_inv.0)
    total.par = length(theta)+length(eta)+1
    
    p.mar <- dbinom(data[[trt.var]], 1, alphas[k])* p[[k]]
    U <- cbind(U_11, rep(0, length(theta)+length(eta)))
    U <- rbind(U, c(-colMeans((data[[obs.y]]*p.mar)*grad_propensity_inv), 
                    -colMeans((data[[obs.y]]*p.mar)*grad_censor_inv), 1))
    V=cbind(V_11, v_21/m)
    V=rbind(V, c(v_21, v_22)/m)
    U_inv=solve(U)
    M=U_inv%*%V%*%t(U_inv)
    Var_avg.outcome[k, ] = c(M.1[total.par, total.par] / m, M.0[total.par, total.par] / m, 
                             M[total.par, total.par] / m)
  }
  colnames(Var_avg.outcome) <- c("treat = 1", "treat = 0", "margin")
  rownames(Var_avg.outcome) <- paste("alpha = ", alphas, sep = "")
  
  
  Var_DE <- data.frame()
  for (k in 1:length(alphas)){
    U_21 <- c(-colMeans((data[[obs.y]]*data[[trt.var]]*p[[k]])*grad_propensity_inv), 
              -colMeans((data[[obs.y]]*(1-data[[trt.var]])*p[[k]])*grad_censor_inv))
    B <- data[[obs.y]]*data[[trt.var]]*p[[k]]/(propensity.score*censor.score) - 
      data[[obs.y]]*(1-data[[trt.var]])*p[[k]]/(propensity.score*censor.score)
    Y <- mean(B)
    
    v_22 <- 0
    v_21 <- rep(0, length(theta)+length(eta))
    for (j in unique(data[[component.var]])) {
      idx <- which(data[[component.var]] == j)
      weight <- m / n
      r   <- weight * sum(B[idx])   - Y
      v_21 <- v_21 + c(
        weight * colSums(propen_v_11[idx, , drop = FALSE]) * r,
        weight * colSums(censor_v_11[idx, , drop = FALSE]) * r
      )
      v_22 <- v_22 + r^2
    }
    
    U <- cbind(U_11, rep(0, length(theta)+length(eta)))
    U <- rbind(U, c(U_21, 1))
    V <- cbind(V_11, v_21/m)
    V <- rbind(V, c(v_21, v_22)/m)
    U_inv <- solve(U)
    M <- U_inv%*%V%*%t(U_inv)
    total.par <- length(theta)+length(eta)+1
    Var_DE <- rbind(Var_DE, c(alphas[k], alphas[k], M[total.par, total.par] / m))
  }
  pairs <- combn(1:length(alphas), 2)
  
  Var_IE <- data.frame()
  for (k in 1:length(alphas)){
    U_21 <- c(-colMeans((data[[obs.y]]*(1-data[[trt.var]])*p[[pairs[2, k]]])*grad_propensity_inv), 
              -colMeans((data[[obs.y]]*(1-data[[trt.var]])*p[[pairs[1, k]]])*grad_censor_inv))
    B <- data[[obs.y]]*(1-data[[trt.var]])*p[[pairs[2, k]]]/(propensity.score*censor.score) -
      data[[obs.y]]*(1-data[[trt.var]])*p[[pairs[1, k]]]/(propensity.score*censor.score)
    Y <- mean(B)
    
    v_22 <- 0
    v_21 <- rep(0, length(theta)+length(eta))
    for (j in unique(data[[component.var]])) {
      idx <- which(data[[component.var]] == j)
      weight <- m / n
      r   <- weight * sum(B[idx])   - Y
      v_21 <- v_21 + c(
        weight * colSums(propen_v_11[idx, , drop = FALSE]) * r,
        weight * colSums(censor_v_11[idx, , drop = FALSE]) * r
      )
      v_22 <- v_22 + r^2
    }
    
    U <- cbind(U_11, rep(0, length(theta)+length(eta)))
    U <- rbind(U, c(U_21, 1))
    V <- cbind(V_11, v_21/m)
    V <- rbind(V, c(v_21, v_22)/m)
    U_inv <- solve(U)
    M <- U_inv%*%V%*%t(U_inv)
    total.par <- length(theta)+length(eta)+1
    Var_IE <- rbind(Var_IE, c(alphas[pairs[1, k]], alphas[pairs[2, k]], M[total.par, total.par] / m))
  }
  
  Var_TE <- data.frame()
  for (k in 1:length(alphas)){
    U_21 <- c(-colMeans((data[[obs.y]]*(data[[trt.var]])*p[[pairs[2, k]]])*grad_propensity_inv), 
              -colMeans((data[[obs.y]]*(1-data[[trt.var]])*p[[pairs[1, k]]])*grad_censor_inv))
    B <- data[[obs.y]]*(data[[trt.var]])*p[[pairs[2, k]]]/(propensity.score*censor.score) -
      data[[obs.y]]*(1-data[[trt.var]])*p[[pairs[1, k]]]/(propensity.score*censor.score)
    Y <- mean(B)
    
    v_22 <- 0
    v_21 <- rep(0, length(theta)+length(eta))
    for (j in unique(data[[component.var]])) {
      idx <- which(data[[component.var]] == j)
      weight <- m / n
      r   <- weight * sum(B[idx])   - Y
      v_21 <- v_21 + c(
        weight * colSums(propen_v_11[idx, , drop = FALSE]) * r,
        weight * colSums(censor_v_11[idx, , drop = FALSE]) * r
      )
      v_22 <- v_22 + r^2
    }
    
    U <- cbind(U_11, rep(0, length(theta)+length(eta)))
    U <- rbind(U, c(U_21, 1))
    V <- cbind(V_11, v_21/m)
    V <- rbind(V, c(v_21, v_22)/m)
    U_inv <- solve(U)
    M <- U_inv%*%V%*%t(U_inv)
    total.par <- length(theta)+length(eta)+1
    Var_TE <- rbind(Var_TE, c(alphas[pairs[1, k]], alphas[pairs[2, k]], M[total.par, total.par] / m))
  }
  
  names(Var_DE) <- c("alpha1", "alpha2", "variance")
  names(Var_IE) <- c("alpha1", "alpha2", "variance")
  names(Var_TE) <- c("alpha1", "alpha2", "variance")
  Var_DE$effect <- "DIRECT"
  Var_IE$effect <- "INDIRECT"
  Var_TE$effect <- "TOTAL"
  
  Var_contrast <- rbind(Var_DE, Var_IE, Var_TE)
  Var_contrast <- Var_contrast[, c(4, 1:3)]
  return(list(Var_avg.outcome, Var_contrast))
  
}
