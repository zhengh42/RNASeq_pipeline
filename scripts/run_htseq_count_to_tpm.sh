salmon_dir=$1
HTSeq_dir=$2
sampleID=$3
glef.py=/srv/gevaertlab/data/Hong/RNASeq/RNASeq_pipeline/scripts/glef.py
tx2gene=/srv/gevaertlab/data/Hong/RNASeq/RNASeq_pipeline/prepare_reference/tx2gene.txt
count2TPM.R=/srv/gevaertlab/data/Hong/RNASeq/RNASeq_pipeline/prepare_reference/count2TPM.R

cat $sampleID | while read LINE; do
      less $HTSeq_dir/$LINE.gene.tab | grep -v '^_'| sed '1i Name\tcount' > $HTSeq_dir/$LINE.gene.tab.tmp;
      $glef.py --ginput $HTSeq_dir/$LINE.gene.tab.tmp --tinput $salmon_dir/$LINE/quant.sf --tgmap $tx2gene --output $HTSeq_dir/$LINE.gene.effectiveLen.count
      Rscript $count2TPM.R $HTSeq_dir/$LINE.gene.effectiveLen.count $HTSeq_dir/$LINE.gene.effectiveLen.count.TPM
done
