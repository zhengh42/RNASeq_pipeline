work_dir=/scratch/users/zhengh42/data/RNASeq/RNASeq_pipeline
kallisto=/home/users/zhengh42/tools/kallisto/kallisto_linux-v0.44.0/kallisto
salmon=/home/users/zhengh42/tools/Salmon/Salmon-0.9.1_linux_x86_64/bin/salmon
star=/home/users/zhengh42/tools/STAR/STAR-2.5.4a/bin/Linux_x86_64/STAR
rsem=/home/users/zhengh42/tools/RSEM/RSEM-1.3.0
subreadbuildindex=/home/users/zhengh42/tools/Subread/subread-1.6.1-Linux-x86_64/bin/subread-buildindex
hisat2build=/home/users/zhengh42/tools/hisat2-2.1.0/hisat2-build
gffread=/home/users/zhengh42/tools/gffread/gffread/gffread
refgenome=GRCh38.primary_assembly.genome.fa

#----------------------------------------------------------------------------------------------
# GENCODE
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.annotation.gtf.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.transcripts.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/GRCh38.primary_assembly.genome.fa.gz
gzip -d gencode.v27.annotation.gtf.gz gencode.v27.transcripts.fa.gz GRCh38.primary_assembly.genome.fa.gz

gencodeV=gencode.v27

less $gencodeV.annotation.gtf | grep -v ^# | awk '$3~/transcript/' | cut -f9 | awk 'OFS="\t"{print $2,$4,$6,$12}' | sed 's/[";]//g' | awk 'OFS="\t"{$5="other";$6="other";if($4~/non_coding|3prime_overlapping_ncRNA|antisense|bidirectional_promoter_lncRNA|lincRNA|macro_lncRNA|sense_intronic|sense_overlapping/)$6="lncRNA";if($4~/protein_coding/)$6="proteincoding";if($3~/protein_coding/)$5="proteincoding";if($3~/non_coding|3prime_overlapping_ncRNA|antisense|bidirectional_promoter_lncRNA|lincRNA|macro_lncRNA|sense_intronic|sense_overlapping/)$5="lncRNA";  print $0}'| sed '1i gene_id\ttranscript\tgenetype1\ttranscripttype1\tgenetype\ttranscripttype' > gene_transcript_type
less gene_transcript_type | sed 1d | cut -f1,3,5 | sort | uniq | sed '1i gene_id\tdetail\ttype'> gene_type
less gene_transcript_type | awk '$5=="lncRNA"' > gene_transcript_type.lncRNA
less gene_transcript_type | awk '$5=="proteincoding"' > gene_transcript_type.proteincoding

less $gencodeV.transcripts.fa | awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }'  | sed 's/^/%/' | sed '$!N;s/\n/\t/g' | sed 's/^%//' > $gencodeV.transcripts.fa.tmp
less $gencodeV.transcripts.fa.tmp | awk -F "|" 'OFS="\t"{print $1,$0}' | sed 's/>//' | Extract.pl gene_transcript_type.lncRNA 2 - 1 | cut -f2- | sed 's/%/\n/'  > $gencodeV.lncRNA.transcripts.fa
less $gencodeV.transcripts.fa.tmp | awk -F "|" 'OFS="\t"{print $1,$0}' | sed 's/>//' | Extract.pl gene_transcript_type.proteincoding 2 - 1 | cut -f2- | sed 's/%/\n/'  > $gencodeV.proteincoding.transcripts.fa

less $gencodeV.annotation.gtf | grep -v ^# | awk '{print $10,$0}' |sed 's/"//;s/"//;s/;/\t/' | Extract.pl gene_transcript_type.lncRNA 1 - 1 | cut -f2- |sed 's/^\s//' > $gencodeV.lncRNA.annotation.gtf
less $gencodeV.annotation.gtf | grep -v ^# | awk '{print $10,$0}' |sed 's/"//;s/"//;s/;/\t/' | Extract.pl gene_transcript_type.proteincoding 1 - 1 | cut -f2- | sed 's/^\s//' > $gencodeV.proteincoding.annotation.gtf

#----------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
# GENCODE+NONCODE
wget http://www.noncode.org/datadownload/NONCODEv5_human_hg38_lncRNA.gtf.gz
gzip -d NONCODEv5_human_hg38_lncRNA.gtf.gz
wget http://www.noncode.org/datadownload/NONCODEv5_human.fa.gz
gzip -d NONCODEv5_human.fa.gz

$gffread NONCODEv5_human_hg38_lncRNA.gtf -g Homo_sapiens_assembly38.fasta -w NONCODEv5_human_hg38_lncRNA.fa
$gffread gencode.v27.annotation.gtf -g $refgenome -w gencode.v27.gffread.fa
cat gencode.v27.transcripts.fa NONCODEv5_human_hg38_lncRNA.fa > GENCODE_NONCODE.transcripts.tmp

less ../output/Salmon/TCGA-AX-A3FS-01A-11R-A22K-07_simul/quant.sf  | cut -f1 | sed 1d | awk '{print $1"\t"$1}' | sed 's/|ENSG[[:graph:]]\+//' > gencode.v27.annotation.gtf.txids
less gencode.v27.annotation.gtf | grep -v ^# | awk '$3!="gene"' | awk '{print $12,$0}' | sed 's/^"//;s/";\s\+/\t/' > gencode.v27.annotation.gtf.tmp
less gencode.v27.annotation.gtf.tmp | match.pl gencode.v27.annotation.gtf.txids 1 - 1 | paste - gencode.v27.annotation.gtf.tmp | cut -f2,4- > gencode.v27.annotation.gtf.tmp1
less NONCODEv5_human_hg38_lncRNA.gtf  | grep -v ^# | awk '{print $12,$0}' | sed 's/^"//;s/";\s\+/\t/' > NONCODEv5_human_hg38_lncRNA.gtf.tmp
cat gencode.v27.annotation.gtf.tmp1 NONCODEv5_human_hg38_lncRNA.gtf.tmp | Exclude.pl GENCODE_NONCODE.salmon_ambiguousStrandKept.idx/duplicate_clusters.tsv 2 - 1 | cut -f2- | awk '$7!="." && $1!~/_/'> GENCODE_NONCODE.annotation.gtf

cat gencode.v27.annotation.gtf.tmp1 NONCODEv5_human_hg38_lncRNA.gtf.tmp | Exclude.pl GENCODE_NONCODE.salmon_ambiguousStrandKept.idx/duplicate_clusters.tsv 2 - 1 |  awk '$8!="." && $2!~/_/ && $4=="transcript"' | cut -f1 > GENCODE_NONCODE.transcripts.ids
less GENCODE_NONCODE.transcripts.tmp |  awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' | sed 'N;s/\n/\t/g;s/>/>\t/' | Extract.pl  GENCODE_NONCODE.transcripts.ids 1 - 2 | sed 's/>\t/>/;s/\t/\n/' > GENCODE_NONCODE.transcripts.fa
#----------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
PRE=gencode.v27
#PRE=GENCODE_NONCODE
### Kallisto index
$kallisto index -i ${PRE}.kallisto.idx ${PRE}.transcripts.fa
$kallisto index -i ${PRE}.kallisto.lncRNA.idx ${PRE}.lncRNA.transcripts.fa
$kallisto index -i ${PRE}.kallisto.proteincoding.idx ${PRE}.proteincoding.transcripts.fa

### Salmon index
$salmon index -i ${PRE}.salmon.idx -t ${PRE}.transcripts.fa
$salmon index -i ${PRE}.salmon.keepDuplicates.idx -t ${PRE}.transcripts.fa --keepDuplicates

$salmon index -i ${PRE}.salmon.lncRNA.keepDuplicates.idx -t ${PRE}.lncRNA.transcripts.fa --keepDuplicates
$salmon index -i ${PRE}.salmon.proteincoding.keepDuplicates.idx -t ${PRE}.proteincoding.transcripts.fa --keepDuplicates

### STAR index
sjdbOverhang=100 # ideal read length - 1
out_dir=${PRE}.STAR.idx.sjdbOverhang${sjdbOverhang}
mkdir -p $work_dir/reference/$out_dir
$star --runThreadN 4 --runMode genomeGenerate --genomeDir $work_dir/reference/$out_dir --genomeFastaFiles $refgenome --sjdbGTFfile ${PRE}.annotation.gtf --sjdbOverhang $sjdbOverhang

out_dir=${PRE}.STAR.lncRNA.idx.sjdbOverhang${sjdbOverhang}
mkdir -p $work_dir/reference/$out_dir
$star --runThreadN 4 --runMode genomeGenerate --genomeDir $work_dir/reference/$out_dir --genomeFastaFiles $refgenome --sjdbGTFfile ${PRE}.lncRNA.annotation.gtf --sjdbOverhang $sjdbOverhang

### rsem index
$rsem/rsem-prepare-reference -p 4 --gtf ${PRE}.annotation.gtf  $refgenome ${PRE}.rsem.idx
$rsem/rsem-prepare-reference -p 4 --gtf ${PRE}.annotation.gtf  --bowtie $refgenome ${PRE}.rsembowtie.idx

### Subread index
$subreadbuildindex -o ${PRE}.subread.idx $refgenome

### HISAT2 index
$hisat2build -f $refgenome ${PRE}.hisat2.idx
#----------------------------------------------------------------------------------------------
