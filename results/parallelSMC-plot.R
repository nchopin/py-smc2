
rm(list = ls())
gc()
graphics.off()
resultsfolder <- "/home/pierre/Dropbox/py-smc2/results/"
resultsfile <- "/tmp/parallelSMCres.RData"
pdffile <- "parallelSMCres.pdf"
library(ggplot2)
setwd(resultsfolder)
cat("working with file", resultsfile, "...
")
load(file = resultsfile)
pdf(file = pdffile, useDingbats = FALSE, title = "parallel SMC results")

dlm <- list("FF" = 1, "GG" = 0.800, "V" = 0.250, "W" = 0.250,
             "m0" = 0, "C0" = 1)


KF <- function(observations, somedlm){
  # the notation comes from the package dlm
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
trueLogLikelihood <- KFLL(observations, dlm)
trueIncrlogLikelihood <- getLoglikelihood(KF(observations, dlm))$IncrLogLike
trueCumLogLikelihood <- cumsum(trueIncrlogLikelihood)

g <- qplot(x = SMCLL, geom = "blank") + geom_histogram(aes(y = ..density..)) + geom_density(alpha = 0.2, fill = "blue") +
    geom_vline(xintercept = trueLogLikelihood, size = 2, colour = "red", linetype = 2) + 
    xlab("log likelihoods")
print(g)

SMCcumloglikelihoods <- apply(allLL, 2, cumsum)
SMCcumlikelihoods <- exp(SMCcumloglikelihoods)
means <- apply(SMCcumlikelihoods, 1, mean)
vars <- apply(SMCcumlikelihoods, 1, var)
normvars <- vars / (exp(trueCumLogLikelihood))^2
T <- length(observations)
g <- qplot(x = 1:T, y = normvars, geom = "line")
g <- g + xlab("time") + ylab("normalised variance")
print(g)
upper <- means + sqrt(vars)
lower <- means - sqrt(vars)
g <- qplot(x = 1:T, y = means, geom = "line") + scale_y_log() + 
    geom_line(aes(y = upper), colour = "green") +
    geom_line(aes(y = lower), colour = "green") +
    geom_point(aes(y = exp(trueCumLogLikelihood)), colour = "red") +
    xlab("time") + ylab("likelihood estimates")
print(g)

dev.off()

