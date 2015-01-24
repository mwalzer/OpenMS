## This is an R script to produce the figures that are attached to the qcML format

#options
options(digits=10)

ms2s_in<-commandArgs(TRUE)[1]
post<-commandArgs(TRUE)[2]
png(post)
##########################
###precursors map as substitute for a MS1 map
##########################

ms2s <- read.table(ms2s_in, header=TRUE, sep="\t", na.strings="?", dec=".", strip.white=TRUE)
names(ms2s) <- c('RT','MZ')
plot(ms2s$RT,ms2s$MZ,pch=16,,xlab="RT (sec)",ylab="m/z",cex=0.3, main="Precursors map")

######################################
dev.off()
