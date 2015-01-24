## This is an R script to produce the figures that are attached to the qcML format

#options
options(digits=10)
require(graphics)

ids_in<-commandArgs(TRUE)[1]
ms2s_in<-commandArgs(TRUE)[2]
post<-commandArgs(TRUE)[3]
png(post)
##########################
###Mass accuracy time course
##########################
ids <- read.table(ids_in, header=TRUE, sep="\t", na.strings="?", dec=".", strip.white=TRUE)
ms2s <- read.table(ms2s_in, header=TRUE, sep="\t", na.strings="?", dec=".", strip.white=TRUE)
names(ms2s) <- c('RT','MZ')
names(ids)[1] <- 'rt'
knime.in<-merge(x = ms2s, y = ids, by.x = "MZ", by.y = "MZ",  all=TRUE)
identified<-knime.in[!is.na(knime.in$PeptideSequence)&is.na(knime.in$rt), ]
plot(identified$"RT"/60,identified$"DeltaPpm",xlab="RT (min)",ylab="mass error in [ppm]",main="", ylim=c(-11,11), pch='.')
lines(lowess(identified$"RT"/60, identified$"DeltaPpm", f=1/20), col="red")
abline(h=0, col="blue")

######################################
dev.off()
