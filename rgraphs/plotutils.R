library(ggplot2)
theme_update(
axis.title.x = theme_text(size=20),
axis.title.y = theme_text(size=20, angle = 90),
 axis.text.x = theme_text(size=15),
 axis.text.y = theme_text(size=15),
 strip.text.x = theme_text(size=15))

loadSomeQuantity <- function(basename, nruns, somequantity, somename = NULL, functiontoapply = NULL){
    if (is.null(somename)){
        if (is.null(dim(get(somequantity)))){
            somename <- somequantity
        } else {
            somename <- paste(somequantity, 1:dim(get(somequantity))[2], sep = "")
        }
    }
    library(foreach)
    df <- foreach(run = 1:nruns, .combine = rbind) %do% {
        load(file = paste(basename, run - 1, ").RData", sep = ""))
        if (is.null(functiontoapply)){
            quantityDF <- as.data.frame(cbind(get(somequantity), 1:T, run))
        } else {
            quantityDF <- as.data.frame(cbind(functiontoapply(get(somequantity)), 1:T, run))
        }
        names(quantityDF) <- c(somename, "Time", "Run")
        quantityDF
    }
    df$Run <- as.factor(df$Run)
    df
}

loadSMCThetaParticles <- function(basename, nruns, time){
    library(foreach)
    df <- foreach(run = 1:nruns, .combine = rbind) %do% {
        load(file = paste(basename, run - 1, ").RData", sep = ""))
        nbparameters <- dim(thetahistory)[2]
        indexhistory <- which(time == savingtimes)
        w <- weighthistory[indexhistory,]
        w <- w / sum(w)
        thetas <- as.data.frame(t(thetahistory[indexhistory,,]))
        thetasDF <- cbind(thetas, w, run, time)
        names(thetasDF) <- c(paste("Theta", 1:(nbparameters), sep = ""), "w", "Run", "Time")
        thetasDF
    }
    df$Run <- as.factor(df$Run)
    df$Time <- as.factor(df$Time)
    df
}

loadMCMCresult <- function(basename, nruns, time, burnin){
    library(foreach)
    df <- foreach(run = 1:nruns, .combine = rbind) %do% {
        load(file = paste(basename, run - 1, ").RData", sep = ""))
        thetasDF <- data.frame(t(chain[,(burnin + 1):(dim(chain)[2])]),  w = 1, run, time)
        names(thetasDF) <- c(paste("Theta", 1:(nbparameters), sep = ""), "w", "Run", "Time")
        thetasDF
    }
    df$Run <- as.factor(df$Run)
    df$w <- 1 / (dim(df)[1])
    df$Time <- as.factor(df$Time)
    df
}


