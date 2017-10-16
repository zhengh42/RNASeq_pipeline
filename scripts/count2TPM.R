#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
data<- read.table(args[1],head=T)
count2Tpm <- function(df)
{
	counts <- df$count
	effLen <- df$EffectiveLength
	rate <- log(counts) - log(effLen)
	denom <- log(sum(exp(rate)))
	exp(rate - denom + log(1e6))
}

data$TPM <- count2Tpm(data)
write.table(data,file=args[2],sep="\t",quote=F,row.names = F)

