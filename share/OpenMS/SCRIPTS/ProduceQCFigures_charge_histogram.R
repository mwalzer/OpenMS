## This is an R script to produce the figures that are attached to the qcML format

#options
options(digits=10)

ids_in<-commandArgs(TRUE)[1]
#ms2s_in<-commandArgs(TRUE)[2]
post<-commandArgs(TRUE)[2] #[3]
png(post)
##########################
###IDs on rt/mz map vs precursors
##########################

ids <- read.table(ids_in, header=TRUE, sep="\t", na.strings="?", dec=".", strip.white=TRUE)
#ms2s <- read.table(ms2s_in, header=TRUE, sep="\t", na.strings="?", dec=".", strip.white=TRUE)
#names(ms2s) <- c('RT','MZ')
#names(ids)[1] <- 'rt'
#knime.in<-merge(x = ms2s, y = ids, by.x = "MZ", by.y = "MZ",  all=TRUE)
#recorded<-knime.in[is.na(knime.in$rt), ]
#identified<-knime.in[!is.na(knime.in$PeptideSequence)&is.na(knime.in$rt), ]
#TODO fetch the ms2s charge as well (from the mgf)
h<-hist(ids$"Charge", xlab="Charge", main="Precursors charge distribution", breaks=seq(0,max(ids$"Charge")), xaxt="n")
axis(side=1,at=h$mids,labels=seq(1,max(ids$"Charge")))
######################################
dev.off()
