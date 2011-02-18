theme_update(
axis.title.x = theme_text(size=20),
axis.title.y = theme_text(size=20, angle = 90),
 axis.text.x = theme_text(size=15),
 axis.text.y = theme_text(size=15))

filteringDF <- data.frame(observations)
filteringDF$year <- 1976:2010
filteringDF$xfiltered <- filteredfirststate
filteringDF$CIlow <- filteredCIlow
filteringDF$CIhigh <- filteredCIhigh
names(filteringDF) <- c("y1", "y2", "year", "x", "CIlow", "CIhigh")
g <- ggplot(filteringDF, aes(x = year))
g <- g + geom_line(aes(y = x), size = 1, linetype = 2, colour = "black")
g <- g + geom_line(aes(y = CIlow), size = 1, colour = "blue")
g <- g + geom_line(aes(y = CIhigh), size = 1, colour = "blue")
g <- g + geom_point(aes(y = y1), size = 4, colour = "green")
g <- g + geom_point(aes(y = y2), size = 4, colour = "red")
g <- g + ylab("Times (seconds)") + ylim(450, 540)
X11(); print(g)

#pdf(file = pdffile, useDingbats = FALSE, title = "SMC2 results")

smoothingDF <- data.frame(observations)
smoothingDF$year <- 1976:2010
smoothingDF$xsmoothed <- smoothedmeansfirststate35
smoothingDF$CIlow <- smoothedmeansCIlow35
smoothingDF$CIhigh <- smoothedmeansCIhigh35
names(smoothingDF) <- c("y1", "y2", "year", "x", "CIlow", "CIhigh")
g <- ggplot(smoothingDF, aes(x = year))
#g <- g + geom_line(aes(y = x), size = 1, linetype = 2, colour = "black")
#g <- g + geom_line(aes(y = CIlow), size = 1, colour = "black")
#g <- g + geom_line(aes(y = CIhigh), size = 1, colour = "black")
g <- g + geom_point(aes(y = y1), size = 4, colour = "green")
g <- g + geom_point(aes(y = y2), size = 4, colour = "red")
g <- g + ylab("Times (seconds)") + xlab("Year")
g <- g + ylim(480, 530)
X11(); print(g)
pdf(file = "athletics-obs.pdf", useDingbats = FALSE)
print(g)
dev.off()


beatDF <- data.frame(observations)
beatDF$year <- 1976:2010
beatDF$smoothbeat1985 <- smoothedmeansbeat198535
beatDF$filterbeat1985 <- filteredbeat1985
beatDF$smoothbeat1993 <- smoothedmeansbeat199335
beatDF$filterbeat1993 <- filteredbeat1993
names(beatDF) <- c("y1", "y2", "year", "smooth1985", "filter1985", "smooth1993", "filter1993")
g <- ggplot(beatDF, aes(x = year))
g <- g + geom_line(aes(y = smooth1985), colour = "blue")
g <- g + geom_point(aes(y = smooth1985), colour = "blue")
g <- g + ylab("Probability of beating XXX seconds")
g <- g + geom_line(aes(y = filter1985), colour = "darkblue")
g <- g + geom_point(aes(y = filter1985), colour = "darkblue")
g <- g + geom_line(aes(y = smooth1993), colour = "green")
g <- g + geom_point(aes(y = smooth1993), colour = "green")
g <- g + geom_line(aes(y = filter1993), colour = "darkgreen")
g <- g + geom_point(aes(y = filter1993), colour = "darkgreen")
g <- g + geom_line(aes(y = smooth1993 / smooth1985), colour = "orange", linetype = 2)
g <- g + geom_point(aes(y = smooth1993 / smooth1985), colour = "orange")
g <- g + scale_y_log10()
#pdf(file = "probabilityOfBeatingXXXseconds.pdf", useDingbats = FALSE, title = "Athletics records")
#X11(); print(g)
print(g)
#dev.off()


indexhistory <- length(savingtimes)
t <- savingtimes[indexhistory]
w <- weighthistory[indexhistory,]
w <- w / sum(w)
thetas <- as.data.frame(t(thetahistory[indexhistory,,]))
thetasDF <- cbind(thetas, w, smoothedvaluesbeat199335, smoothedvaluesbeat198535)
names(thetasDF) <- c(paste("Theta", 1:(nbparameters), sep = ""), "w", "smooth1993", "smooth1985")
#pdf(file = "probaVersusParameters.pdf", useDingbats = FALSE, title = "SMC2 results")

g <- ggplot(thetasDF, aes(x = thetasDF[[1]], y = smooth1985, weight = w))
g <- g + geom_point(aes(alpha = w, size = w), colour = "blue", shape = 1)
g <- g + geom_point(aes(alpha = w, size = w, y = smooth1993), colour = "red", shape = 2)
g <- g + opts(legend.position = "none")
g <- g + xlab(expression(nu)) + ylab("Probability")
g <- g + scale_y_log()
g <- g + xlim(0, 5)
#pdf(file = "probaVersusNu.pdf", useDingbats = FALSE)
print(g)
ggsave("probaVersusNu.png", width = 7, height = 7, dpi = 150)
#dev.off()


g <- ggplot(thetasDF, aes(x = -thetasDF[[2]], y = smooth1985, weight = w))
g <- g + geom_point(aes(alpha = w, size = w), colour = "blue", shape = 1)
g <- g + geom_point(aes(alpha = w, size = w, y = smooth1993), colour = "red", shape = 2)
g <- g + opts(legend.position = "none")
g <- g + xlab(expression(xi)) + ylab("Probability")
g <- g + xlim(-0.8, 0)
g <- g + scale_y_log()
print(g)
ggsave("probaVersusXi.png", width = 7, height = 7, dpi = 150)

g <- ggplot(thetasDF, aes(x = thetasDF[[3]], y = smooth1985, weight = w))
g <- g + geom_point(aes(alpha = w, size = w), colour = "blue", shape = 1)
g <- g + geom_point(aes(alpha = w, size = w, y = smooth1993), colour = "red", shape = 2)
g <- g + opts(legend.position = "none")
g <- g + xlab(expression(sigma)) + ylab("Probability")
g <- g + xlim(2, 7)
g <- g + scale_y_log()
#pdf(file = "probaVersusSigma.pdf", useDingbats = FALSE)
print(g)
ggsave("probaVersusSigma.png", width = 7, height = 7, dpi = 150)
#dev.off()




