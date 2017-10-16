#######
# Download data from GDC data portal 
#######
# 1. Query in GDC data portal Repository: cases.project.project_id in ["TCGA-BRCA","TCGA-CESC","TCGA-COAD","TCGA-ESCA","TCGA-HNSC","TCGA-KIRC","TCGA-KIRP","TCGA-LAML","TCGA-LGG","TCGA-LIHC","TCGA-LUAD","TCGA-LUSC","TCGA-PAAD","TCGA-PCPG","TCGA-PRAD","TCGA-READ","TCGA-SARC","TCGA-SKCM","TCGA-STAD","TCGA-TGCT","TCGA-THCA","TCGA-THYM","TCGA-UCEC"] and files.data_format in ["BAM"] and files.experimental_strategy in ["RNA-Seq"]
# 2. Download the JSON file, covert to CSV format in http://www.convertcsv.com/json-to-csv.htm, and get sample.annotation.csv


library(dplyr)
library(data.table)
metd <- read.csv("sample.annotation.csv",head=T,sep = ",")
manifest <- read.table("gdc_manifest.2017-06-29T16-16-02.567550.tsv",head=T)

names(metd)[names(metd) == 'file_name'] <- 'filename'

manifest.metd<-dplyr::left_join(manifest,metd,by="filename")
manifest.metd<-manifest.metd[order(manifest.metd$cases.0.project.project_id),]

summary(manifest.metd$size)
manifest.metd<- manifest.metd %>% filter(size>=quantile(manifest.metd$size)[2] & size<=quantile(manifest.metd$size)[4] & cases.0.project.project_id != "TCGA-ESCA" & cases.0.project.project_id != 'TCGA_THCA') %>% droplevels()

tcga.freq<-as.data.frame(table(manifest.metd$cases.0.project.project_id))

index <- NULL
for (i in 1:nrow(tcga.freq)){
  n<-as.numeric(tcga.freq$Freq[i])
  a<-rep(0,n)
  a[sample(n,10)] = 1
  index<-c(index,a)
}

manifest.metd$index<-index
table(manifest.metd %>% filter(index==1) %>% select(cases.0.project.project_id))

tmp<-rep(1:10,21)
manifest.metd.chosen<- manifest.metd %>% filter(index==1) %>% select(c(names(manifest),"cases.0.project.project_id"))
manifest.metd.chosen<- manifest.metd.chosen %>% mutate(sampleID=paste(cases.0.project.project_id,tmp,sep=""))

write.table(manifest.metd.chosen,"gdc_manifest_TCGA210samples.txt",quote = F,row.names = F,sep = "\t")

