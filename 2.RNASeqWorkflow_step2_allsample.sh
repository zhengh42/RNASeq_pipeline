#!/bin/bash

organismname="homosapiens"			### "homosapiens/" for human and "musmusculus" for mouse
transcriptome="gencodev32"	     		### "gencodev32" or "gencodev32noncodev5" for human, and "gencodevM24" for mouse
workDir="/srv/gevaertlab/data/Hong/RNASeq"
projectid="rnaseqproject"
refDir="/srv/gevaertlab/reference"
RscriptsDir="/srv/gevaertlab/data/Hong/RNASeq/RNASeq_pipeline/github/RNASeq_pipeline/R"

summary_kallisto=0				### 1 means running the step
run_deseq2fromkallisto=0

###-------------------------------------------------------------------------------------------------------------------
### Kallisto summary
###-------------------------------------------------------------------------------------------------------------------
if [[ $summary_kallisto == 1 ]]; then
	bash run_stats.sh $workDir $projectid $transcriptome
	Rscript --vanilla run_tximport.R $workDir $projectid $organismname $transcriptome $refDir
fi
###-------------------------------------------------------------------------------------------------------------------
###-------------------------------------------------------------------------------------------------------------------



###-------------------------------------------------------------------------------------------------------------------
### DESeq from Kallisto
###-------------------------------------------------------------------------------------------------------------------
if [[ $run_deseq2fromkallisto == 1 ]]; then
	Rscript --vanilla run_deseq2.R $workDir $projectid $organismname $transcriptome $RscriptsDir "tximport" "kallisto"
fi
###-------------------------------------------------------------------------------------------------------------------
###-------------------------------------------------------------------------------------------------------------------



