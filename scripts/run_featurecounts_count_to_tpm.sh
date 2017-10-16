salmon_dir=$1
featureCounts_dir=$2
sampleID=$3
glef.py=/srv/gevaertlab/data/Hong/RNASeq/RNASeq_pipeline/scripts/glef.py
tx2gene=/srv/gevaertlab/data/Hong/RNASeq/RNASeq_pipeline/prepare_reference/tx2gene.txt
count2TPM.R=/srv/gevaertlab/data/Hong/RNASeq/RNASeq_pipeline/prepare_reference/count2TPM.R

cat $sampleID | while read LINE; do
      less $featureCounts_dir/$LINE.out | grep -v '^_'| sed '1i Name\tcount' > $featureCounts_dir/$LINE.out.tmp;
      $glef.py --ginput $featureCounts_dir/$LINE.out.tmp --tinput $salmon_dir/$LINE/quant.sf --tgmap $tx2gene --output $featureCounts_dir/$LINE.gene.effectiveLen.count
      Rscript $count2TPM.R $featureCounts_dir/$LINE.gene.effectiveLen.count $featureCounts_dir/$LINE.gene.effectiveLen.count.TPM
done

