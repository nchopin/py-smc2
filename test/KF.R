# rm(list = ls())
### DLM:
# m0, C0 : mean and variance of x_0
# x_t \vert x_{t-1} ~ Normal(G x_{t-1}, W)
# y_t \vert x_t ~ Normal(F x_t, V)
#
# cf http://www.jstatsoft.org/v36/i12/

rm(list = ls())
gc()
graphics.off()
resultsfolder <- "/home/pierre/Dropbox/py-smc2/results/"
resultsfile <- "/tmp/testSMCres.RData"
pdffile <- "testSMCres.pdf"
library(ggplot2)
setwd(resultsfolder)
cat("working with file", resultsfile, "...
")
load(file = resultsfile)

KF <- function(observations, somedlm){
  # the notation comes from Dynamical Linear Models with R
  T <- length(observations)
  m <- rep(0, T + 1)
  C <- rep(1, T + 1)
  a <- rep(0, T)
  R <- rep(0, T)
  f <- rep(0, T)
  Q <- rep(0, T)
  m[1] <- somedlm$m0
  C[1] <- somedlm$C0
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

dlm <- list("FF" = 1, "GG" = 0.8, "V" = 0.25, "W" = 0.25,
             "m0" = 0, "C0" = 1)
kalmanresults <- KF(observations, dlm)
sometime <- 37
truefilteringdistribution <- function(x){
  dnorm(x, mean = kalmanresults$FiltStateMean[sometime], 
            sd = sqrt(kalmanresults$FiltStateVar[sometime]))
}
qplot(x = xhistory[,,sometime + 1], geom = "blank") + 
#   geom_histogram(aes(y = ..density..), binwidth = 0.01) + 
  stat_function(fun=truefilteringdistribution,colour = "red") + 
  geom_density(fill = "blue", alpha = 0.3)

xhat <- apply(X=xhistory[,1,2:(T+1)], MARGIN = 2, FUN=mean)
qplot(x = 1:T, y = xhat, geom = "line") +
 geom_line(aes(y = truestates), colour= "green") + 
 geom_point(aes(y = kalmanresults$FiltStateMean), alpha = 0.5,
            colour = "orange")



