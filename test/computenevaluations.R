Ntheta
nevaluations <- 0
currentnx <- Nxlist[1]
for (t in 1:T){
  if (t %in% increaseindices){
    cat("increase of Nx at time", t, "\n")
    currentnx <- Nxlist[match(t, increaseindices)]
  }
  nevaluations <- nevaluations + currentnx
}

nevaluations * Ntheta / T