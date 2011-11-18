
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

print(cumlikelihood(0.8, 25))
alternateevidences <- rep(0, T)
for (time in 1:T){
  marginaltheta <- function(theta){
    sapply(X=theta, FUN=function(theta) cumlikelihood(theta, time))
  }
  alternateevidences[time] <- integrate(f=marginaltheta, lower = 0, upper = 1)$value
}
alternateevidences
cumprod(trueevidences)

qplot(x = 1:T, y = alternateevidences, geom = "line") + geom_line(aes(y = cumprod(trueevidences))) +
  scale_y_log()


cumlogalternateevidences <- log(alternateevidences)
cumlogSMC2evidences <- cumsum(log(SMC2evidences)) - cumlogalternateevidences
cumlogBSMCevidences <- cumsum(log(BSMCevidences)) - cumlogalternateevidences
g <- qplot(x = 1:T, y = cumlogSMC2evidences, geom = "line", linetype = "SMC2 evidence",
           colour = "SMC2 evidence")
g <- g + geom_line(aes(y = cumlogBSMCevidences, linetype = "BSMC evidence", colour = "BSMC evidence"))
pdf(file = "cumlogevidencecomparison.pdf", width = 10, height = 6)
print(g)
dev.off()

