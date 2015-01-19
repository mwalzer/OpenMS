## This is an R script to produce the figures that are attached to the qcML format

#options
options(digits=10)

##########################
###IDs on rt/mz map vs precursors
##########################

ids <- read.table("/tmp/ids.csv", header=TRUE, sep="\t", na.strings="?", dec=".", strip.white=TRUE)
ms2s <- read.table("/tmp/ms2s.csv", header=TRUE, sep="\t", na.strings="?", dec=".", strip.white=TRUE)
names(ms2s) <- c('RT','Precursor')
#names(ids) <- c('RT','MZ', 'Score', 'PeptideSequence', 'Charge', 'TheoreticalWeight', 'DeltaPpm')
knime.in<-merge(x = ms2s, y = ids, by.x = "RT", by.y = "RT",  all=TRUE)
recorded<-knime.in[knime.in$Precursor != "?", ]
plot(recorded$RT/60,recorded$Precursor,pch=16,,xlab="RT (min)",ylab="m/z",cex=0.3)
identified<-knime.in[knime.in$MZ != "?", ]
points(identified$RT/60,identified$MZ,col="red",pch=4,cex=0.3)
legend("topleft",c("recorded spectra","identified spectra"),pch=19,col=c(1,2))

######################################
dev.off()
