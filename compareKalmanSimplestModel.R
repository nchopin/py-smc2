
rm(list = ls())
gc()
graphics.off()
resultsfolder <- "/home/pierre/Dropbox/py-smc2/results/"
resultsfile <- "/home/pierre/Dropbox/py-smc2/results/temp.RData"
library(ggplot2)
cat("working with file", resultsfile, "...
")
load(file = resultsfile)
dlm <- list("FF" = 1, "GG" = 0.8, "V" = 0.25, "W" = 0.25,
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
kalmanresults <- KF(observations, dlm)
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
trueLogLikelihood <- KFLL(observations, dlm)
trueIncrlogLikelihood <- getLoglikelihood(KF(observations, dlm))$IncrLogLike
trueCumLogLikelihood <- cumsum(trueIncrlogLikelihood)
###################################
# posterior distribution of the parameters
if (POSTERIOR){
  temptrueloglikelihood <- function(theta){
     somedlm <- dlm
     somedlm["GG"] <- theta
     return(KFLL(observations, somedlm))
  }
  trueloglikelihood <- function(theta){
     return(sapply(X= theta, FUN= temptrueloglikelihood))
  }
  trueunnormlikelihood <- function(theta) exp(trueloglikelihood(theta))
  normlikelihood <- integrate(f = trueunnormlikelihood, lower = 0, upper = 1)$value
  truelikelihood <- function(theta) trueunnormlikelihood(theta) / normlikelihood
  priorfunction <- function(x){
      return(1)
  }
  trueposterior <- function(x) priorfunction(x) * truelikelihood(x)
  
  indexhistory <- length(savingtimes)
  w <- weighthistory[indexhistory,]
  w <- w / sum(w)
  i <- 1
  g <- qplot(x = thetahistory[indexhistory,i,], weight = w, geom = "blank") + 
    geom_histogram(aes(y = ..density..)) + geom_density(fill = "blue", alpha = 0.5) +
      xlab(expression(rho))
  g <- g + geom_vline(xintercept = trueparameters[i], linetype = 2, size = 1)
  
  g <- g + stat_function(fun = priorfunction, aes(colour = "prior"), linetype = 1, size = 1)
  g <- g + stat_function(fun = trueposterior, aes(colour = "posterior"), size = 2)
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

