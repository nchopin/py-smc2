
additionalsfilename = strsplit(pdffile, ".pdf")[[1]]
pdf(file = paste(additionalsfilename, "-additionals.pdf", sep = ""))
g <- ggplot(thetasDF, aes(x = thetasDF[[4]], y = thetasDF[[6]], weight = w))
g <- g + geom_point(aes(size = w), alpha = 0.2) + opts(legend.position = "none")
g <- g + geom_density2d()
g <- g + xlab("r") + ylab(expression(theta))
print(g)


g <- ggplot(thetasDF, aes(x = thetasDF[[3]], y = thetasDF[[5]], weight = w))
g <- g + geom_point(aes(size = w), alpha = 0.2) + opts(legend.position = "none")
g <- g + geom_density2d()
g <- g + xlab(expression(log(N[0]))) + ylab("K")
print(g)


g <- ggplot(thetasDF, aes(x = thetasDF[[6]], y = thetasDF[[5]], weight = w))
g <- g + geom_point(aes(size = w), alpha = 0.2) + opts(legend.position = "none")
g <- g + geom_density2d()
g <- g + xlab(expression(theta)) + ylab("K")
print(g)

g <- ggplot(thetasDF, aes(x = thetasDF[[4]], y = thetasDF[[5]], weight = w))
g <- g + geom_point(aes(size = w), alpha = 0.2) + opts(legend.position = "none")
g <- g + geom_density2d()
g <- g + xlab("r") + ylab("K")
print(g)
dev.off()


