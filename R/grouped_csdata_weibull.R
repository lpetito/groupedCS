##########
### Generate data from a weibull distribution, then censor uniformly
##########
gen.data.weibull.unif <- function(n, k, shape, scale, quantile = 0.99, alpha = 1, beta = 1){
  require(gtools)
  #Generate grouping
  groups <- permute(rep(1:ceiling(n/k), length.out = n))
  #Generate Ts and Cs
  Ts <- rweibull( n, shape, scale )
  Cs <- qweibull(quantile, shape, scale)*runif(n)
  delta.ind <- as.numeric(Cs > Ts)
  
  #Create dataset
  data = data.frame(cbind(groups, Ts, Cs, delta.ind))
  data <- data[order(data$Cs),]
  
  initial.p <- runif(n)
  initial.p <- initial.p[order(initial.p)]
  data$initial.p <- initial.p
  #Create grouped deltas
  data<- data[order(data$groups, data$Cs),]
  Delta.groups <- sapply(1:ceiling(n/k), FUN = function(X){
    temp <- delta.ind[groups==X]
    ifelse(sum(temp) != 0, 1, 0)
  })
  data$delta.group <- Delta.groups[data$groups]
  #Create misclassified grouped ys
  if (alpha != 1 | beta != 1){
    #Create misclassified individual ys
    data$y.ind <- sapply(1:n, FUN = function(x) {
      ifelse(data$delta.ind[x]==1,
             rbinom(1, 1, alpha), rbinom(1, 1, 1-beta))
    })
    Y.groups <- sapply(1:ceiling(n/k), FUN = function(X){
      temp <- data$delta.group[data$groups==X]
      z <- ifelse(sum(temp) !=0, 1, 0) #Indicator of group 1 or group 0
      #Introduce misclassification: if 
      ifelse(z == 1, #logical statement evaluated
             rbinom(1, 1, prob = alpha), #command evaluated if true
             rbinom(1, 1, prob = 1 - beta)) #command evaluated if false
    })
    data$y.group <- Y.groups[data$groups]
    return(data)
  }
  else{
    data$y.ind <- data$delta.ind
    data$y.group <- data$delta.group
    return(data)
  }
}

#PAVA for individual data
pava.cs <- function(Cs, initial){
  require(isotone)
  gpava(Cs, initial, ties = "secondary")$x
}

#Altered PAVA for misclassified current status data
pava.cs.mc <- function(Cs, initial, alpha, beta){
  require(isotone)
  Gs <- gpava(Cs, initial, ties = "secondary")$x
  Gs[Gs <= (1 - beta)] <- 1 - beta
  Gs[Gs >= alpha] <- alpha
  Fs <- (Gs + beta - 1) / (alpha + beta - 1)
  return(Fs)
}

#Expectation function
expectation <- function(p, delta, alpha=1, beta=1){
  Fs = p
  Ss = 1 - p
  gamma = alpha + beta - 1
  k = length(p)
  if (sum(delta) != 0){
    update <- sapply(1:length(Fs), FUN = function(X){
      alpha * Fs[X] / (alpha - gamma * prod(Ss))
    })
  }
  if (sum(delta) == 0){
    update <- sapply(1:length(Fs), FUN = function(X){
      (1 - alpha) * Fs[X] / (1 - alpha + gamma * prod(Ss))
    })
  }
  if (!(sum(delta) %in% c(0, length(p)))){
    stop("group test result not all 0 or 1")
  }
  return(update)
}

#Hybrid EM-PAV algorithm
hybrid.em.pav <- function(data, initial.p, threshold=0.01, alpha=1, beta=1){
  require(isotone)
  delta <- data$y.group
  if(alpha == beta & alpha == 1) {delta <- data$delta.group}
  diff = 1
  it = 0
  p.0 = initial.p
  m = max(as.numeric(data$groups))
  while(diff > threshold){
    p.star = do.call("c", lapply(1:m, FUN = function(X){
      expectation(p.0[data$groups==X], delta[data$groups==X], alpha, beta)
    }))
    out = gpava(data$Cs, p.star, ties="secondary")$x
    it = it + 1
    diff = sqrt(sum((out - p.0)^2))
    p.0 = out
  }
  p.star = out
  return(list(results = out, num.iterations = it, diff = diff))
}

##########
### Create simulation function
##########

simulation <- function(n, k, shape, scale, quantile, x, alpha=1, beta=1, t=0.01){
  data <- gen.data.weibull.unif(n, k, shape, scale, quantile, alpha, beta)
  desc.ind <- with(data, xtabs(~delta.ind + y.ind))
  desc.group <- with(data, xtabs(~delta.group + y.group)) / k
  Cs = seq(0, x, by=0.05)
  #For grouped data
  group.result <- hybrid.em.pav(data, data$initial.p, t, alpha, beta)
  num.it <- group.result$num.iterations
  diff <- group.result$diff
  temp <- data.frame(Cs = data$Cs, group.result = group.result$result)
  temp <- temp[order(temp$Cs),]
  result <- stepfun(temp$Cs, c(0, temp$group.result), right=F)
  group.res <- result(Cs)
  #For individual data
  data <- data[order(data$Cs),]
  if (alpha == 1 & beta ==1){
    ind.result <- with(data, pava.cs(Cs, delta.ind))
  }
  else {
    ind.result <- with(data, pava.cs.mc(Cs, y.ind, alpha, beta))
  }
  result <- stepfun(data$Cs, c(0, ind.result), right=F)
  ind.res <- result(Cs)
  
  return(list(desc.ind = desc.ind, desc.group = desc.group,
              num.it = num.it,
              ind.result = ind.res,
              group.result = group.res))
}