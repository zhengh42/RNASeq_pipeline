#!/bin/bash

workDir=$1
projectid=$2
transcriptome=$3
#workDir=/srv/gevaertlab/data/Hong/RNASeq
#projectid=202002_juneho_micecellline

### kallisto
ls $workDir/$projectid/kallisto/*${transcriptome}/*json | sed 's/.gen.*//;s/.*\///' > $workDir/$projectid/kallisto/kallisto.reads.sampleID
cat $workDir/$projectid/kallisto/*${transcriptome}/*json | egrep 'n_processed' | awk '{print $2}' | tr -d ',' > $workDir/$projectid/kallisto/kallisto.reads.total
cat $workDir/$projectid/kallisto/*${transcriptome}/*json | egrep 'n_pseudoaligned' | awk '{print $2}' | tr -d ',' > $workDir/$projectid/kallisto/kallisto.reads_a.aligned
cat $workDir/$projectid/kallisto/*${transcriptome}/*json | egrep 'p_pseudoaligned' | awk '{print $2}' | tr -d ',' > $workDir/$projectid/kallisto/kallisto.reads_a.aligned.percentage
cat $workDir/$projectid/kallisto/*${transcriptome}/*json | egrep 'n_unique' | awk '{print $2}' | tr -d ',' > $workDir/$projectid/kallisto/kallisto.reads_a.unique
cat $workDir/$projectid/kallisto/*${transcriptome}/*json | egrep 'p_unique' | awk '{print $2}' | tr -d ',' > $workDir/$projectid/kallisto/kallisto.reads_a.unique.percentage

paste $workDir/$projectid/kallisto/kallisto.reads.sampleID $workDir/$projectid/kallisto/kallisto.reads.total $workDir/$projectid/kallisto/kallisto.reads_* | sed '1i sample\tReadsTotal\tReadsAlignedG\tPCTAlignedG\tReadsUniqueG\tPCTUniqueG' > $workDir/$projectid/kallisto/summary.readsstats.txt


