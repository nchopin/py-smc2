# The file requires two RData files to be available. They should be
# runs from SMC^2 and BSMC on the same model, same observations, same horizon T.
rm(list = ls())
gc()
graphics.off()
resultsfile <- 
"/home/pierre/Dropbox/py-smc2/results/BSMC-simplestmodel-synthetic-T250-N10000-h0.100(0).RData"
load(file = resultsfile)
BSMCevidences <- evidences
rm(list = ls()[ls() != "BSMCevidences"])
resultsfile <- 
"/home/pierre/Dropbox/py-smc2/results/SMC2-simplestmodel-synthetic-T250-independent-dynamicNx250-Ntheta1000(0).RData"
load(file = resultsfile)
SMC2evidences <- evidences
rm(list = ls()[!(ls() %in% c("BSMCevidences", "SMC2evidences", "T", "observations"))])

library(ggplot2)
dlm <- list("FF" = 1, "GG" = 0.800, "V" = 0.250, "W" = 0.250,
             "m0" = 0, "C0" = 1)

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

likelihood <- function(y, theta, time){
  changeddlm <- dlm
  changeddlm$GG <- theta
  KFres <- KF(observations[1:time], changeddlm)
  return(dnorm(y, mean = KFres$NextObsMean[time], sd = sqrt(KFres$NextObsVar[time])))
}
trueevidences <- rep(0, T)
for (time in 1:T){
  marginaltheta <- function(theta){
    sapply(X=theta, FUN=function(theta) likelihood(observations[time], theta, time))
  }
  trueevidences[time] <- integrate(f=marginaltheta, lower = 0, upper = 1)$value
}

g <- qplot(x = 1:T, y = SMC2evidences, geom = "line", linetype = "SMC2 evidence",
           colour = "SMC2 evidence")
g <- g + geom_line(aes(y = trueevidences, linetype = "true evidence", colour = "true evidence"))
g <- g + geom_line(aes(y = BSMCevidences, linetype = "BSMC evidence", colour = "BSMC evidence"))
g <- g + xlab("iterations") + ylab("evidence") + scale_colour_discrete(name = "") +
  scale_linetype_discrete(name = "")
# pdf(file = "evidencecomparison.pdf", width = 10, height = 6)
print(g)
# dev.off()

cumlogtrueevidences <- cumsum(log(trueevidences))
cumlogSMC2evidences <- cumsum(log(SMC2evidences)) - cumlogtrueevidences
cumlogBSMCevidences <- cumsum(log(BSMCevidences)) - cumlogtrueevidences
g <- qplot(x = 1:T, y = cumlogSMC2evidences, geom = "line", linetype = "SMC2 evidence",
           colour = "SMC2 evidence")
g <- g + geom_line(aes(y = cumlogBSMCevidences, linetype = "BSMC evidence", colour = "BSMC evidence"))
# pdf(file = "cumlogevidencecomparison.pdf", width = 10, height = 6)
print(g)
# dev.off()
