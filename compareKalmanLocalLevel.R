
rm(list = ls())
gc()
graphics.off()
resultsfolder <- "/home/pierre/Dropbox/py-smc2/results/"
resultsfile <- "/home/pierre/Dropbox/py-smc2/results/temp.RData"
library(ggplot2)
cat("working with file", resultsfile, "...
")
load(file = resultsfile)
trueparameters
dlm <- list("FF" = 1, "GG" = 1, "V" = trueparameters[2], "W" = trueparameters[1],
             "m0" = 0, "C0" = 1)
POSTERIOR <- FALSE
FILTERING <- TRUE

KF <- function(observations, somedlm){
  T <- length(observations)
  m <- rep(0, T + 1); C <- rep(1, T + 1)
  a <- rep(0, T); R <- rep(0, T)
  f <- rep(0, T); Q <- rep(0, T)
  m[1] <- somedlm$m0; C[1] <- somedlm$C0
  for (t in 1:T){
    a[t] <- somedlm$GG * m[t]
    R[t] <- somedlm$GG * C[t] * somedlm$GG + somedlm$W
    f[t] <- somedlm$FF * a[t]
    Q[t] <- somedlm$FF * R[t] * somedlm$FF + somedlm$V
    m[t+1] <- a[t] + R[t] * somedlm$FF * (1 / Q[t]) * (observations[t] - f[t])
    C[t+1] <- R[t] - R[t] * somedlm$FF * (1 / Q[t]) * somedlm$FF * R[t]
  }
  return(list(observations = observations, NextObsMean = f, NextObsVar = Q,
              NextStateMean = a, NextStatevar = R,
              FiltStateMean = m[2:(T+1)], FiltStateVar = C[2:(T+1)]))
}
getLoglikelihood <- function(KFresults){
  IncrLogLike <- log(dnorm(KFresults$observations, 
            mean = KFresults$NextObsMean, 
            sd = sqrt(KFresults$NextObsVar)))
  loglikelihood <- sum(IncrLogLike)
  return(list(IncrLogLike = IncrLogLike, loglikelihood = loglikelihood))
}
KFLL <- function(observations, dlm){
  KFres <- KF(observations, dlm)
  return(getLoglikelihood(KFres)$loglikelihood)
}

###################################
# posterior distribution of the parameters
if (POSTERIOR){
  temptrueloglikelihood <- function(theta){
      somedlm <- dlm
      somedlm["V"] <- theta[2]
      somedlm["W"] <- theta[1]
      return(KFLL(observations, somedlm))
  }
  trueloglikelihood <- function(theta1, theta2){
    theta <- cbind(theta1, theta2)
    return(apply(X=theta, MARGIN=1, FUN=temptrueloglikelihood))
  }
  priorfunction <- function(x){
      shape <- 1.00000 
      scale <- 1.00000
      return(scale**shape / gamma(shape) * x**(- shape - 1) * exp(-scale / x))
  }
  truelogposterior <- function(theta1, theta2){
    prior <- log(priorfunction(theta1)) + log(priorfunction(theta2))
    theta <- cbind(theta1, theta2)
    return(prior + apply(X=theta, MARGIN=1, FUN=temptrueloglikelihood))
  }
  range1 <- range(thetahistory[1,1,])
  range2 <- range(thetahistory[1,2,])
  x = seq(from = range1[1], to = range1[2], length.out = 6)
  y = seq(from = range2[1], to = range2[2], length.out = 6)
  z = outer(X=x,Y=y,FUN=truelogposterior)
  maxlogposterior <- max(z)
  tmpunnormmarginal1 <- function(theta1){
    temp <- function(theta2){
      exp(truelogposterior(theta1, theta2) - maxlogposterior)
    }
    integrate(temp, lower = range2[1], upper = range2[2])$value
  }
  unnormmarginal1 <- function(theta1) sapply(X=theta1, FUN=tmpunnormmarginal1)
  normconstant1 <- integrate(unnormmarginal1, lower = range1[1], upper = range1[2])$value
  normmarginal1 <- function(theta1) unnormmarginal1(theta1) / normconstant1
  tmpunnormmarginal2 <- function(theta2){
    temp <- function(theta1){
      exp(truelogposterior(theta1, theta2) - maxlogposterior)
    }
    integrate(temp, lower = range1[1], upper = range1[2])$value
  }
  unnormmarginal2 <- function(theta2) sapply(X=theta2, FUN=tmpunnormmarginal2)
  normconstant2 <- integrate(unnormmarginal2, lower = range2[1], upper = range2[2])$value
  normmarginal2 <- function(theta2) unnormmarginal2(theta2) / normconstant2
  marginals <- c(normmarginal1, normmarginal2)
  
  indexhistory <- length(savingtimes)
  w <- weighthistory[indexhistory,]
  w <- w / sum(w)
  i <- 1
  g <- qplot(x = thetahistory[indexhistory,i,], weight = w, geom = "blank") + 
    geom_histogram(aes(y = ..density..)) + geom_density(fill = "blue", alpha = 0.5) +
      xlab(expression(rho))
  g <- g + geom_vline(xintercept = trueparameters[i], linetype = 2, size = 1)
  
  g <- g + stat_function(fun = priorfunction, aes(colour = "prior"), linetype = 1, size = 1)
  g <- g + stat_function(fun = marginals[[i]], aes(colour = "posterior"), size = 2)
  g <- g + scale_colour_discrete(name = "")
  print(g)
  i <- 2
  g <- qplot(x = thetahistory[indexhistory,i,], weight = w, geom = "blank") + 
    geom_histogram(aes(y = ..density..)) + geom_density(fill = "blue", alpha = 0.5) +
      xlab(expression(rho))
  g <- g + geom_vline(xintercept = trueparameters[i], linetype = 2, size = 1)
  
  g <- g + stat_function(fun = priorfunction, aes(colour = "prior"), linetype = 1, size = 1)
  g <- g + stat_function(fun = marginals[[i]], aes(colour = "posterior"), size = 2)
  g <- g + scale_colour_discrete(name = "")
  print(g)
}
###################################
# filtering
if (FILTERING){
  truestates <- as.matrix(truestates, ncol = 1)
  kalmanresults <- KF(observations, dlm)
  g <- qplot(x = 1:T, y = filteredstate1, geom = "line", colour = "SMC2 mean")
  g <- g + geom_line(aes(y = kalmanresults$FiltStateMean, colour = "KF mean"), alpha = 1.)
  g <- g + geom_point(aes(y = kalmanresults$FiltStateMean, colour = "KF mean"))
  g <- g + xlab("time") + ylab("hidden states")
  g <- g + scale_colour_discrete(name = "")
  print(g)
}




