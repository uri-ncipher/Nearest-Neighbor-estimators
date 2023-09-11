num_neighbors=function(net, node){
  return(length(neighbors(net, node)))
}
trt_neighbors=function(net, node){
  return(sum(subset(data, data$id %in% neighbors(net, node))$treatment))
}
avg_neighbors=function(net, node, variable){
  return(mean(subset(data, data$id %in% neighbors(net, node))[, variable]))
}

##################################################################################
################## IPW1 functions ################################################
##################################################################################
propensity=function(theta){
  prop = foreach(i = 1:n, .combine = "c") %dopar% {
    subdata=subset(data, data$id %in% neighborhood(net0, order = 1, nodes = i)[[1]])
    integrand=function(b){
      p=1
      for (j in 1:nrow(subdata)){
        pred=c(theta[1], theta[2:(2+length(base_covariate)-1)]) %*% c(1, unlist(subdata[j, base_covariate]))
        p=p*dbinom(subdata$treatment[j], size = 1, prob = plogis(as.vector(pred)+as.vector(b)))
      }
      return(p*dnorm(b, mean = 0, sd=theta[length(theta)]))
    }
    integrate(integrand, -Inf, Inf)$value
  }
  return(prop)
}
pi=function(alpha){
  return(alpha^data$na_a*(1-alpha)^(data$na-data$na_a))
}
Y_IPW=function(alpha, score){
  p=pi(alpha)
  return(c(mean(data$outcome*data$treatment*p/score), 
           mean(data$outcome*(1-data$treatment)*p/score), 
           mean(data$outcome*dbinom(data$treatment, 1, alpha)*p/score)))
}
Y_DE=function(alpha, score){
  p=pi(alpha)
  return(mean(data$outcome*data$treatment*p/score)- 
           mean(data$outcome*(1-data$treatment)*p/score))
}
Y_IE=function(alpha, score){
  p0=pi(alpha[1])
  p1=pi(alpha[2])
  return(mean(data$outcome*(1-data$treatment)*p0/score)-
           mean(data$outcome*(1-data$treatment)*p1/score))
}
Y_TE=function(alpha, score){
  p0=pi(alpha[1])
  p1=pi(alpha[2])
  return(mean(data$outcome*(data$treatment)*p0/score)-
           mean(data$outcome*(1-data$treatment)*p1/score))
}
Y_OE=function(alpha, score){
  p0=dbinom(data$treatment, 1, alpha[1])*pi(alpha[1])
  p1=dbinom(data$treatment, 1, alpha[2])*pi(alpha[2])
  return(mean(data$outcome*p0/score)-
           mean(data$outcome*p1/score))
}

Var_Y=function(a, alpha, U_11, v_11, V_11, grad_propensity_inv, score, theta){
  p=pi(alpha)
  if (a==1){
    U_21=-colMeans((data$outcome*data$treatment*p)*grad_propensity_inv)
    B=data$outcome*data$treatment*p/score
    Y=mean(B)
  } else {
    U_21=-colMeans((data$outcome*(1-data$treatment)*p)*grad_propensity_inv)
    B=data$outcome*(1-data$treatment)*p/score
    Y=mean(B)
  }
  v_22=0
  v_21=rep(0, length(theta))
  for (j in 1:m){
    c=data$id[data$component==j]
    r=(m/n)*sum(B[c])-Y
    v_21=v_21+(m/n)*colSums(v_11[c,])*r
    v_22=v_22+r^2
  }
  U=cbind(U_11, rep(0, length(theta)))
  U=rbind(U, c(U_21, 1))
  V=cbind(V_11, v_21/m)
  V=rbind(V, c(v_21, v_22)/m)
  U_inv=solve(U)
  M=U_inv%*%V%*%t(U_inv)
  return(M[(length(theta)+1), (length(theta)+1)]/m)
}
Var_Y_margin=function(alpha, U_11, v_11, V_11, grad_propensity_inv, score, theta){
  p=pi(alpha)
  v_22=0
  v_21=rep(0, length(theta))
  B=data$outcome*dbinom(data$treatment, 1, alpha)*p/score
  Y=mean(B)
  for (j in 1:m){
    c=data$id[data$component==j]
    r=(m/n)*sum(B[c])-Y
    v_21=v_21+(m/n)*r*colSums(v_11[c,])
    v_22=v_22+r^2
  }
  U=cbind(U_11, rep(0, length(theta)))
  U=rbind(U, c(-colMeans((data$outcome*dbinom(data$treatment, 1, alpha)*p)*grad_propensity_inv), 1))
  V=cbind(V_11, v_21/m)
  V=rbind(V, c(v_21, v_22)/m)
  U_inv=solve(U)
  M=U_inv%*%V%*%t(U_inv)
  return(M[(length(theta)+1), (length(theta)+1)]/m)
}

Var_DE=function(alpha, U_11, v_11, V_11, grad_propensity_inv, score, theta){
  p=pi(alpha)
  U_21=-colMeans((data$outcome*(2*data$treatment-1)*p)*grad_propensity_inv)
  B=data$outcome*(2*data$treatment-1)*p/score
  Y=mean(B)
  v_22=0
  v_21=rep(0, length(theta))
  for (j in 1:m){
    c=data$id[data$component==j]
    r=(m/n)*sum(B[c])-Y
    v_21=v_21+(m/n)*colSums(v_11[c,])*r
    v_22=v_22+r^2
  }
  U=cbind(U_11, rep(0, length(theta)))
  U=rbind(U, c(U_21, 1))
  V=cbind(V_11, v_21/m)
  V=rbind(V, c(v_21, v_22)/m)
  U_inv=solve(U)
  M=U_inv%*%V%*%t(U_inv)
  return(M[(length(theta)+1), (length(theta)+1)]/m)
}
Var_IE=function(alpha, U_11, v_11, V_11, grad_propensity_inv, score, theta){
  p0=pi(alpha[1])
  p1=pi(alpha[2])
  U_21=-colMeans(data$outcome*(1-data$treatment)*(p0-p1)*grad_propensity_inv)
  B=data$outcome*(1-data$treatment)*(p0-p1)/score
  Y=mean(B)
  v_22=0
  v_21=rep(0, length(theta))
  for (j in 1:m){
    c=data$id[data$component==j]
    r=(m/n)*sum(B[c])-Y
    v_21=v_21+(m/n)*colSums(v_11[c,])*r
    v_22=v_22+r^2
  }
  U=cbind(U_11, rep(0, length(theta)))
  U=rbind(U, c(U_21, 1))
  V=cbind(V_11, v_21/m)
  V=rbind(V, c(v_21, v_22)/m)
  U_inv=solve(U)
  M=U_inv%*%V%*%t(U_inv)
  return(M[(length(theta)+1), (length(theta)+1)]/m)
}
Var_TE=function(alpha, U_11, v_11, V_11, grad_propensity_inv, score, theta){
  p0=pi(alpha[1])
  p1=pi(alpha[2])
  U_21=-colMeans(data$outcome*(data$treatment*p0-(1-data$treatment)*p1)*grad_propensity_inv)
  B=data$outcome*(data$treatment*p0-(1-data$treatment)*p1)/score
  Y=mean(B)
  v_22=0
  v_21=rep(0, length(theta))
  for (j in 1:m){
    c=data$id[data$component==j]
    r=(m/n)*sum(B[c])-Y
    v_21=v_21+(m/n)*colSums(v_11[c,])*r
    v_22=v_22+r^2
  }
  U=cbind(U_11, rep(0, length(theta)))
  U=rbind(U, c(U_21, 1))
  V=cbind(V_11, v_21/m)
  V=rbind(V, c(v_21, v_22)/m)
  U_inv=solve(U)
  M=U_inv%*%V%*%t(U_inv)
  return(M[(length(theta)+1), (length(theta)+1)]/m)
}
Var_OE=function(alpha, U_11, v_11, V_11, grad_propensity_inv, score, theta){
  p0=dbinom(data$treatment, 1, alpha[1])*pi(alpha[1])
  p1=dbinom(data$treatment, 1, alpha[2])*pi(alpha[2])
  v_22=0
  v_21=rep(0, length(theta))
  B=data$outcome*(p0-p1)/score
  Y=mean(B)
  for (j in 1:m){
    c=data$id[data$component==j]
    r=(m/n)*sum(B[c])-Y
    v_21=v_21+(m/n)*r*colSums(v_11[c,])
    v_22=v_22+r^2
  }
  U=cbind(U_11, rep(0, length(theta)))
  U=rbind(U, c(-colMeans(data$outcome*(p0-p1)*grad_propensity_inv), 1))
  V=cbind(V_11, v_21/m)
  V=rbind(V, c(v_21, v_22)/m)
  U_inv=solve(U)
  M=U_inv%*%V%*%t(U_inv)
  return(M[(length(theta)+1), (length(theta)+1)]/m)
}

Var=function(alpha, U_11, v_11, V_11, grad_propensity_inv, score, theta){
  return(c(Var_Y(1, alpha, U_11, v_11, V_11, grad_propensity_inv, score, theta), 
           Var_Y(0, alpha, U_11, v_11, V_11, grad_propensity_inv, score, theta), 
           Var_Y_margin(alpha, U_11, v_11, V_11, grad_propensity_inv, score, theta)))
}

IPW_1_model=function(data, base_covariate, alpha){
  nn_data=data.frame()
  for (i in 1:n){
    sub=subset(data, data$id %in% neighborhood(net0, order=1, nodes = i)[[1]])
    sub$group=i
    nn_data=rbind(nn_data, sub)
  }
  formula=paste("treatment~", paste(base_covariate, collapse = "+"), "+", "(1|group)")
  M=glmer(formula, data = nn_data, family = binomial(link= "logit"))
  theta=c(fixef(M), as.data.frame(VarCorr(M))$sdcor)
  score=propensity(c(fixef(M), as.data.frame(VarCorr(M))$sdcor))
  U_11 = foreach(i=1:n, .combine = "+") %dopar% {
    log_propensity_ind=function(theta){
      subdata=subset(data, data$id %in% neighborhood(net0, order = 1, nodes = i)[[1]])
      integrand=function(b){
        p=1
        for (j in 1:nrow(subdata)){
          pred=c(theta[1], theta[2:(2+length(base_covariate)-1)]) %*% c(1, unlist(subdata[j, base_covariate]))
          p=p*dbinom(subdata$treatment[j], size = 1, prob = plogis(as.vector(pred)+as.vector(b)))
        }
        return(p*dnorm(b, mean = 0, sd=theta[length(theta)]))
      }
      return(log(integrate(integrand, -Inf, Inf)$value))
    }
    -(1/n)*hessian(log_propensity_ind, x=theta)
  }
  v_11 = foreach(i=1:n, .combine = "rbind") %dopar% {
    log_propensity_ind=function(theta){
      subdata=subset(data, data$id %in% neighborhood(net0, order = 1, nodes = i)[[1]])
      integrand=function(b){
        p=1
        for (j in 1:nrow(subdata)){
          pred=c(theta[1], theta[2:(2+length(base_covariate)-1)]) %*% c(1, unlist(subdata[j, base_covariate]))
          p=p*dbinom(subdata$treatment[j], size = 1, prob = plogis(as.vector(pred)+as.vector(b)))
        }
        return(p*dnorm(b, mean = 0, sd=theta[length(theta)]))
      }
      return(log(integrate(integrand, -Inf, Inf)$value))
    }
    grad(log_propensity_ind, x=theta)
  }
  grad_propensity_inv = foreach (i = 1:n, .combine = "rbind") %dopar% {
    propensity_ind_inv=function(theta){
      subdata=subset(data, data$id %in% neighborhood(net0, order = 1, nodes = i)[[1]])
      integrand=function(b){
        p=1
        for (j in 1:nrow(subdata)){
          pred=c(theta[1], theta[2:(2+length(base_covariate)-1)]) %*% c(1, unlist(subdata[j, base_covariate]))
          p=p*dbinom(subdata$treatment[j], size = 1, prob = plogis(as.vector(pred)+as.vector(b)))
        }
        return(p*dnorm(b, mean = 0, sd=theta[length(theta)]))
      }
      return(1/integrate(integrand, -Inf, Inf)$value)
    }
    grad(propensity_ind_inv, x=theta)
  }
  V_11 = foreach (j = 1:m, .combine = "+") %do% {
    colSums(v_11[data$id[data$component==j],])%*%t(colSums(v_11[data$id[data$component==j],]))
  }
  V_11=V_11*(m/n^2)
  
  APO=rbind(cbind(ldply(alpha, Y_IPW, score=score), alpha=alpha, 
                  type="point estimate"), 
            cbind(ldply(alpha, Var, U_11, v_11, V_11, 
                        grad_propensity_inv, score=score, theta=theta), 
                  alpha=alpha, type="variance"))
  names(APO)=c("a=1", "a=0", "margin", "alpha", "type")
  contrast=t(combn(alpha,2)) 
  
  CE=rbind(cbind(ldply(alpha, Y_DE, score=score), 
                 alpha0=alpha, alpha1=alpha, type="Direct"), 
           cbind(ldply(alpha, Var_DE, U_11=U_11, v_11=v_11, V_11=V_11, 
                       grad_propensity_inv=grad_propensity_inv, score=score, 
                       theta=theta), 
                 alpha0=alpha, alpha1=alpha, type="Var DE"),
           cbind(adply(contrast, 1, Y_IE, score=score), 
                 alpha0=contrast[, 1], alpha1=contrast[, 2], type="Indirect")[, -1],
           cbind(adply(contrast, 1, Var_IE, U_11=U_11, v_11=v_11, V_11=V_11, 
                       grad_propensity_inv=grad_propensity_inv, score=score, 
                       theta=theta), 
                 alpha0=contrast[, 1], alpha1=contrast[, 2], type="Var IE")[, -1],
           cbind(adply(contrast, 1, Y_TE, score=score), 
                 alpha0=contrast[, 1], alpha1=contrast[, 2], type="Total")[, -1],
           cbind(adply(contrast, 1, Var_TE, U_11=U_11, v_11=v_11, V_11=V_11, 
                       grad_propensity_inv=grad_propensity_inv, score=score, 
                       theta=theta), 
                 alpha0=contrast[, 1], alpha1=contrast[, 2], type="Var TE")[, -1],
           cbind(adply(contrast, 1, Y_OE, score=score), 
                 alpha0=contrast[, 1], alpha1=contrast[, 2], type="Overall")[, -1],
           cbind(adply(contrast, 1, Var_OE, U_11=U_11, v_11=v_11, V_11=V_11, 
                       grad_propensity_inv=grad_propensity_inv, score=score, 
                       theta=theta), 
                 alpha0=contrast[, 1], alpha1=contrast[, 2], type="Var OE")[, -1]
  )
  names(CE)=c("estimation", "alpha0", "alpha1", "type")
  return(list(APO, CE))
}

##################################################################################
################## IPW2 functions ################################################
##################################################################################

propensity_2=function(M1, M2){
  p1=predict(M1, data, type = "response")
  p2=predict(M2, data, type = "response")
  return((p1^data$na_a)*((1-p1)^data$notna_a)*dbinom(data$treatment, size = 1, prob = p2))
}

IPW_2_model=function(data, M1,M2, alpha){
  
  score=propensity_2(M1, M2)
  theta=c(M1$coefficients, M2$coefficients)
  
  U_11 = foreach(i = 1:n, .combine = "+") %dopar% {
    log_propensity_2_ind=function(theta){
      v1=theta[1:length(M1$coefficients)]
      v2=theta[(length(M1$coefficients)+1):length(theta)]
      p1=plogis(v1 %*% unlist(c(1, data[i, c("treatment", base_covariate, avg_covariate)])))
      p2=plogis(v2 %*% unlist(c(1, data[i, base_covariate])))
      #prop=dbinom(data$na_a[i], size = data$na[i], prob = p1)*dbinom(data$treatment[i], size = 1, prob = p2)
      prop=(p1^data$na_a[i])*((1-p1)^data$notna_a[i])*dbinom(data$treatment[i], size = 1, prob = p2)
      return(log(prop))
    }
    -(1/n)*hessian(log_propensity_2_ind, x=theta)
  }
  v_11 = foreach(i = 1:n, .combine = "rbind") %dopar% {
    log_propensity_2_ind=function(theta){
      v1=theta[1:length(M1$coefficients)]
      v2=theta[(length(M1$coefficients)+1):length(theta)]
      p1=plogis(v1 %*% unlist(c(1, data[i, c("treatment", base_covariate, avg_covariate)])))
      p2=plogis(v2 %*% unlist(c(1, data[i, base_covariate])))
      #prop=dbinom(data$na_a[i], size = data$na[i], prob = p1)*dbinom(data$treatment[i], size = 1, prob = p2)
      prop=(p1^data$na_a[i])*((1-p1)^data$notna_a[i])*dbinom(data$treatment[i], size = 1, prob = p2)
      return(log(prop))
    }
    grad(log_propensity_2_ind, x=theta)
  }
  grad_propensity_inv = foreach (i = 1:n, .combine = "rbind") %dopar% {
    propensity_2_ind_inv = function(theta){
      v1=theta[1:length(M1$coefficients)]
      v2=theta[(length(M1$coefficients)+1):length(theta)]
      p1=plogis(v1 %*% unlist(c(1, data[i, c("treatment", base_covariate, avg_covariate)])))
      p2=plogis(v2 %*% unlist(c(1, data[i, base_covariate])))
      #prop=dbinom(data$na_a[i], size = data$na[i], prob = p1)*dbinom(data$treatment[i], size = 1, prob = p2)
      prop=(p1^data$na_a[i])*((1-p1)^data$notna_a[i])*dbinom(data$treatment[i], size = 1, prob = p2)
      return(1/prop)
    }
    grad(propensity_2_ind_inv, x=theta)
  }
  V_11 = foreach(j = 1:m, .combine = "+") %dopar% {
    colSums(v_11[data$id[data$component==j],])%*%t(colSums(v_11[data$id[data$component==j],]))
  }
  V_11=V_11*(m/n^2)
  APO=rbind(cbind(ldply(alpha, Y_IPW, score=score), alpha=alpha, 
                  type="point estimate"), 
            cbind(ldply(alpha, Var, U_11, v_11, V_11, 
                        grad_propensity_inv, score=score, theta=theta), 
                  alpha=alpha, type="variance"))
  names(APO)=c("a=1", "a=0", "margin", "alpha", "type")
  contrast=t(combn(alpha,2)) 
  
  CE=rbind(cbind(ldply(alpha, Y_DE, score=score), 
                 alpha0=alpha, alpha1=alpha, type="Direct"), 
           cbind(ldply(alpha, Var_DE, U_11, v_11, V_11, 
                       grad_propensity_inv, score=score, theta=theta), 
                 alpha0=alpha, alpha1=alpha, type="Var DE"),
           cbind(adply(contrast, 1, Y_IE, score=score), 
                 alpha0=contrast[, 1], alpha1=contrast[, 2], type="Indirect")[, -1],
           cbind(adply(contrast, 1, Var_IE, U_11, v_11, V_11, 
                       grad_propensity_inv, score=score, theta=theta), 
                 alpha0=contrast[, 1], alpha1=contrast[, 2], type="Var IE")[, -1],
           cbind(adply(contrast, 1, Y_TE, score=score), 
                 alpha0=contrast[, 1], alpha1=contrast[, 2], type="Total")[, -1],
           cbind(adply(contrast, 1, Var_TE, U_11, v_11, V_11, 
                       grad_propensity_inv, score=score, theta=theta), 
                 alpha0=contrast[, 1], alpha1=contrast[, 2], type="Var TE")[, -1],
           cbind(adply(contrast, 1, Y_OE, score=score), 
                 alpha0=contrast[, 1], alpha1=contrast[, 2], type="Overall")[, -1],
           cbind(adply(contrast, 1, Var_OE, U_11, v_11, V_11, 
                       grad_propensity_inv, score=score, theta=theta), 
                 alpha0=contrast[, 1], alpha1=contrast[, 2], type="Var OE")[, -1]
  )
  names(CE)=c("estimation", "alpha0", "alpha1", "type")
  return(list(APO, CE))
}
