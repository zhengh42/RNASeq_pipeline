args = commandArgs(trailingOnly=TRUE)
library(tximport)

workDir=args[1]
projectid=args[2]
organismname=args[3]        
transcriptome=args[4]      
refDir=args[5] 

#workDir="/srv/gevaertlab/data/Hong/RNASeq"
#projectid="202002_juneho_micecellline"
#organismname="musmusculus" 
#transcriptome="gencodevM24"
#refDir="/srv/gevaertlab/reference"

tx2gene<-read.table(paste0(refDir,"/",organismname,"/",transcriptome,".tx2gene"))
samplemeta=read.table(paste0(workDir,"/",projectid,"/","samples.meta"),header = T)
sampleID=samplemeta$sampleID

files <- paste0(workDir,"/",projectid,"/kallisto/","/",sampleID,".",transcriptome,"/abundance.h5")
all(file.exists(files))
names(files)<-sampleID 

kallisto.expr<- tximport(files,type="kallisto",tx2gene = tx2gene,countsFromAbundance = "lengthScaledTPM")

saveRDS(kallisto.expr,file=paste0(workDir,"/",projectid,"/kallisto/tximport.",transcriptome,".rds"))
