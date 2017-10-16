library(readr)
library(tximport)

samplefile<-args[1]
kallisto.biasC.dir<-args[2]
kallisto.nobiasC.dir<-args[3]
salmon.biasC.dir<-args[4]
salmon.nobiasC.dir<-args[5]

samples <- scan(samplefile,what="character",quiet=TRUE)

### kallisto
files<-file.path(kallisto.biasC.dir,samples,"abundance.tsv")
filesa<-file.path(kallisto.nobiasC.dir,samples,"abundance.tsv")
all(file.exists(files))
all(file.exists(filesa))
names(files)<-paste("biasC",samples,sep=".")
names(filesa)<-paste("nobiasC",samples,sep=".")

tx2gene<-read.table("/srv/gevaertlab/data/Hong/RNASeq/RNASeq_pipeline/prepare_reference/tx2gene.txt")
kallisto.txim.biasC <- tximport(files,type="kallisto",tx2gene = tx2gene)
kallisto.txim.nobiasC <- tximport(filesa,type="kallisto",tx2gene = tx2gene)

write.csv(kallisto.txim.biasC,file=paste0(kallisto.biasC.dir,"/kallisto.biasC.txt"),row.names = T)
write.csv(kallisto.txim.nobiasC,file=paste0(kallisto.nobiasC.dir,"/kallisto.nobiasC.txt"),row.names = T)

### salmon
files<-file.path(salmon.biasC.dir,samples,"quant.sf")
filesa<-file.path(salmon.nobiasC.dir,samples,"quant.sf")
all(file.exists(files))
all(file.exists(filesa))
names(files)<-paste("biasC",samples,sep=".")
names(filesa)<-paste("nobiasC",samples,sep=".")

salmon.txim.biasC <- tximport(files,type="salmon",tx2gene = tx2gene)
salmon.txim.nobiasC <- tximport(filesa,type="salmon",tx2gene = tx2gene)

write.csv(salmon.txim.biasC,file=paste0(salmon.biasC.dir,"/salmon/salmon.biasC.txt"),row.names = T)
write.csv(salmon.txim.nobiasC,file=paste0(salmon.nobiasC.dir,"/salmon.nobiasC.txt"),row.names = T)



