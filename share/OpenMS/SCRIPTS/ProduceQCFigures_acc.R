## This is an R script to produce the figures that are attached to the qcML format

#options
options(digits=10)

ids_in<-commandArgs(TRUE)[1]
post<-commandArgs(TRUE)[2]
png(post)
##########################
###Mass accuracy
##########################
knime.in <- read.table(file=ids_in, header=TRUE, sep="\t", na.strings="?", dec=".", strip.white=TRUE)
names(knime.in) <- c('RT','MZ', 'Score', 'PeptideSequence', 'Charge', 'TheoreticalWeight', 'DeltaPpm', 'Oxidation(M)')
hist(knime.in$"DeltaPpm",xlim=c(-10,10),breaks=seq(min(knime.in$"DeltaPpm")-0.01, max(knime.in$"DeltaPpm")+0.01, 0.01),xlab="ppm", main="")
abline(v=median(knime.in$"DeltaPpm"),col="red", lwd=2)
mtext(paste("median(accuracy)=",round(median(knime.in$"DeltaPpm"),3)," ppm",sep=""))
######################################
dev.off()
