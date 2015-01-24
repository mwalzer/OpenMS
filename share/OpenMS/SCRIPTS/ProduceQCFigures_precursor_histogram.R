## This is an R script to produce the figures that are attached to the qcML format

#options
options(digits=10)

ms2s_in<-commandArgs(TRUE)[1]
post<-commandArgs(TRUE)[2]
png(post)
##########################
###IDs on rt/mz map vs precursors
##########################

ms2s <- read.table(ms2s_in, header=TRUE, sep="\t", na.strings="?", dec=".", strip.white=TRUE)
names(ms2s) <- c('RT','MZ')
hist(ms2s$"RT",breaks=100, xlab="RT (sec)", main="Precursors histogram")

######################################
dev.off()
