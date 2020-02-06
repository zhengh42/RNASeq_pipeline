require(data.table)
require(gridExtra)
require(knitr)
require(ComplexHeatmap)
require(circlize)
require(RColorBrewer)
require(matrixStats)
require(ggplot2)
require(grid)
require(gridExtra)
require(reshape2)
require(circlize)
require(RColorBrewer)
require(extrafont)
require(captioner)

args = commandArgs(trailingOnly=TRUE)

workDir=args[1]
projectid=args[2]
organismname=args[3]        
transcriptome=args[4]      
scriptDir=args[5] 
datatype=args[6] # count or from tximport
tool=args[7] # kallisto

source(paste0(scriptDir,"/RNASeqScripts.R"))

#workDir="/srv/gevaertlab/data/Hong/RNASeq"
#projectid="202002_juneho_micecellline"
#organismname="musmusculus" 
#transcriptome="gencodevM24"
#scriptDir="/home/zhengh42/R/scripts"
#datatype="tximport"
#tool="kallisto"

obj=readRDS(paste0(workDir,"/",projectid,"/",tool,"/tximport.",transcriptome,".rds"))
samplemeta=read.table(paste0(workDir,"/",projectid,"/","samples.meta"),header = T)
coldata=data.frame(condition=samplemeta$condition)
rownames(coldata)=samplemeta$sampleID

if(datatype=="tximport"){
  obj$length[which(obj$length==0)]=1
  DESeq2.obj=DESeq2Analysis(obj,coldata,type=datatype,normmethod = "vst")
}

if(datatype=="count"){
  colnames(obj)[-1]=substr(colnames(obj)[-1],1,16)
  obj.count=as.matrix(obj[,c(sample.case,sample.case2,sample.case3,sample.ctrl)])
  rownames(obj.count)=obj$gene_id
  DESeq2.obj=DESeq2Analysis(obj.count,samplemeta,type=datatype,normmethod = "vst")
}
                   
saveRDS(DESeq2.obj,file=paste0(workDir,"/",projectid,"/",tool,"/deseq2.",transcriptome,".rds"))
                   
                   
