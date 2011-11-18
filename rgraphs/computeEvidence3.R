# The file requires two RData files to be available. They should be
# runs from SMC^2 and BSMC on the same model, same observations, same horizon T.
rm(list = ls())
gc()
graphics.off()
resultsfolder <- "/home/pierre/Dropbox/py-smc2/results/simplestmodel/synthetic/"
setwd(resultsfolder)
nsimBSMC  <- 6
T <- 1000
BSMCevidences <- matrix(ncol = nsimBSMC, nrow = T)

baseBSMCname <- "BSMC-T1000-N100000-h0.100"
for (counter in 1:nsimBSMC){
  resultsfile <- paste(baseBSMCname, "(", (counter-1), ").RData", sep = "")
  load(file = resultsfile)
  BSMCevidences[,counter] <- evidences
}  
listtokeep <- c("BSMCevidences", "T", "nsimBSMC", "listtokeep", "nsimBSMC")
rm(list = ls()[!(ls() %in% listtokeep)])
nsimSMC2 <- 6
SMC2evidences <- matrix(ncol = nsimSMC2, nrow = T)
baseSMC2name <- "SMC2-T1000-independent-dynamicNx250-Ntheta500"
for (counter in 1:nsimSMC2){
  resultsfile <- paste(baseSMC2name, "(", (counter-1), ").RData", sep = "")
  load(file = resultsfile)
  SMC2evidences[,counter] <- evidences
}
listtokeep <- c(listtokeep, "SMC2evidences", "observations", "nsimSMC2")
rm(list = ls()[!(ls() %in% listtokeep)])

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

cumlikelihood <- function(theta, time){
  changeddlm <- dlm
  changeddlm$GG <- theta
  KFresults <- KF(observations[1:time], changeddlm)
  IncrLogLike <- log(dnorm(KFresults$observations[1:time], 
            mean = KFresults$NextObsMean[1:time], 
            sd = sqrt(KFresults$NextObsVar[1:time])))
  loglikelihood <- sum(IncrLogLike)
  exp(loglikelihood)
}

# alternateevidences <- rep(0, T)
# for (time in 1:T){
#   cat(time/T, "% -")
#   marginaltheta <- function(theta){
#     sapply(X=theta, FUN=function(theta) cumlikelihood(theta, time))
#   }
#   alternateevidences[time] <- integrate(f=marginaltheta, lower = 0, upper = 1)$value
# }; cat("\n")
# save(alternateevidences, file = "truthT1000.RData")
load(file = "truthT1000.RData")
alternateevidences <- alternateevidences[1:T]
qplot(x = 1:T, y = alternateevidences, geom = "line") +
  scale_y_log()
cumlogalternateevidences <- log(alternateevidences)
cumlogSMC2evidences <- apply(log(SMC2evidences), 2, cumsum)
cumlogBSMCevidences <- apply(log(BSMCevidences), 2, cumsum)
cumlogSMC2evidences <- apply(cumlogSMC2evidences, 2, function(x) x - cumlogalternateevidences)
cumlogBSMCevidences <- apply(cumlogBSMCevidences, 2, function(x) x - cumlogalternateevidences)

g <- qplot(x = 1:T, geom = "blank")
for (counter in 1:nsimSMC2){
  g <- g + geom_line(aes_string(y = paste("cumlogSMC2evidences[,", counter, "]", sep = ""), 
                                linetype = '"SM2 evidence"', colour = '"SMC2 evidence"'))
}
for (counter in 1:nsimBSMC){
  g <- g + geom_line(aes_string(y = paste("cumlogBSMCevidences[,", counter, "]", sep = ""), 
                                linetype = '"BSMC evidence"', colour = '"BSMC evidence"'))
}
g <- g + scale_colour_discrete(name = "") + scale_linetype_discrete(name = "")
g <- g + xlim(0, 600)
print(g)

