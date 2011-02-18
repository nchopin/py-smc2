
rm(list = ls())
gc()
graphics.off()
resultsfolder <- "/home/pierre/Dropbox/pySMC2/results/"
Multiresultsfile <- "/home/pierre/Dropbox/pySMC2/results/SVmultifactor/SP500recentyears/SMC2-T1006-dynamicNx125-Nth1000(0).RData"
library(ggplot2)
setwd(resultsfolder)
theme_update(
axis.title.x = theme_text(size=20),
axis.title.y = theme_text(size=20, angle = 90),
 axis.text.x = theme_text(size=15),
 axis.text.y = theme_text(size=15))

load(file = Multiresultsfile)
Multievidencedataframe <- as.data.frame(cbind(1:length(evidences), evidences, cumsum(log(evidences)))) 
names(Multievidencedataframe) <- c("iterations", "evidence", "cumevidence")

g <- ggplot(data = Multievidencedataframe, aes(x = iterations, y= evidence))
g <- g + geom_line() + xlab("iterations") + ylab("evidence")
g <- g + scale_y_log()
#print(g)
g <- ggplot(data = Multievidencedataframe, aes(x = iterations, y= cumevidence))
g <- g + geom_line() + xlab("iterations") + ylab("cumevidence")
X11(); print(g)

Oneresultsfile <- "/home/pierre/Dropbox/pySMC2/results/SVonefactor/SP500recentyears/SMC2-T1006-dynamicNx100-Nth1000(0).RData"
load(file = Oneresultsfile)
Oneevidencedataframe <- as.data.frame(cbind(1:length(evidences), evidences, cumsum(log(evidences)))) 
names(Oneevidencedataframe) <- c("iterations", "evidence", "cumevidence")

g <- ggplot(data = Oneevidencedataframe, aes(x = iterations, y= evidence))
g <- g + geom_line() + xlab("iterations") + ylab("evidence")
g <- g + scale_y_log()
#print(g)
g <- ggplot(data = Oneevidencedataframe, aes(x = iterations, y= cumevidence))
g <- g + geom_line() + xlab("iterations") + ylab("cumevidence")
g <- g + geom_line(data = Multievidencedataframe, aes(x = iterations, y = cumevidence))
X11(); print(g)

