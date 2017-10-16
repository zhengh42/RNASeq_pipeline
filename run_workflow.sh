work_dir=/srv/gevaertlab/data/Hong/RNASeq/RNASeq_pipeline
ref_dir=/srv/gevaertlab/data/Hong/RNASeq/RNASeq_pipeline/prepare_reference
id=$1 # ID in manifest file, i.e. 5cbff956-2abb-4bad-bc92-9b9aea5787e4
ida=$2 # Manully assigned sample ID, i.e. TCGA_BRCA1

#######
# Download bam files from GDC data portal
#######
cd $work_dir/TCGA_RNASeq/seq
token=[the token from GDC]
id=$1
source /home/zhengh42/tools/gdc-client/gdc-client/venv/bin/activate
/home/zhengh42/tools/gdc-client/gdc-client/bin/gdc-client download $id -t $token 1>$work_dir/logs/$id.download.log 2>&1

#######
# Convert to fastq
#######
pre=`basename ./$id/*bam | cut -f1 -d "_"`

docker -v $work_dir/TCGA_RNASeq/seq:/mnt run zhengh42/picard:2.13.2 java -jar picard.jar SamToFastq \
        INPUT=/mnt/$id/${pre}_gdc_realn_rehead.bam \
        FASTQ=/mnt/${ida}_1.fq \
        SECOND_END_FASTQ=/mnt/${ida}_2.fq \
        VALIDATION_STRINGENCY=LENIENT \
        UNPAIRED_FASTQ=/mnt/${ida}_unpaired.fq \
	1>$work_dir/logs/$ida.bamtofastq.log 2>&1

gzip ${ida}_1.fq ${ida}_2.fq
rm ./$id/${pre}_gdc_realn_rehead.ba*

#######
# QC by trim-galore
#######
docker run -v $work_dir/TCGA_RNASeq/seq:/home zhengh42/trim-galore:0.4.4  \
        trim_galore -q 15 --stringency 3 --gzip --length 15 --paired 
        /home/${ida}_1.fq.gz /home/${ida}_2.fq.gz --fastqc --output_dir /home 
        1> $work_dir/logs/$ida.trim_galore.log 2>&1

#######
# Kallisto
#######
### with sequence bias estimation
docker run -v $work_dir/TCGA_RNASeq/seq:/home/seq -v $ref_dir:/home/ref -v $work_dir/TCGA_RNASeq/kallisto/biasC:/home/out \
        zhengh42/kallisto:0.43.0 \
        kallisto quant -i /home/ref/gencode.v25.kallisto.idx -t 4 -b 100 \
        -o /home/out/$ida /home/seq/${ida}_1_val_1.fq.gz /home/seq/${ida}_2_val_2.fq.gz  --bias  \
        1> $work_dir/logs/$ida.kallisto.log 2>&1

### without sequence bias estimation
docker run -v $work_dir/TCGA_RNASeq/seq:/home/seq -v $ref_dir:/home/ref -v $work_dir/TCGA_RNASeq/kallisto/nobiasC:/home/out \
        zhengh42/kallisto:0.43.0 \
        kallisto quant -i /home/ref/gencode.v25.kallisto.idx -t 4 -b 100 \
        -o /home/out/$ida /home/seq/${ida}_1_val_1.fq.gz /home/seq/${ida}_2_val_2.fq.gz  \
        1> $work_dir/logs/$ida.kallisto.log 2>&1

#######
# Salmon
#######
### with sequence and gc bias estimation
docker run -v $work_dir/TCGA_RNASeq/seq:/home/seq -v $ref_dir:/home/ref -v $work_dir/TCGA_RNASeq/salmon/biasC:/home/out \
        zhengh42/salmon:0.8.2 \
        salmon quant -i /home/ref/gencode.v25.salmon.idx -l A -p 4 --numBootstraps 100 \
        -o /home/out/$ida -1 /home/seq/${ida}_1_val_1.fq.gz -2 /home/seq/${ida}_2_val_2.fq.gz --seqBias --gcBias \
        1> $work_dir/logs/$ida.salmon.log 2>&1

### without sequence and gc bias estimation
docker run -v $work_dir/TCGA_RNASeq/seq:/home/seq -v $ref_dir:/home/ref -v $work_dir/TCGA_RNASeq/salmon/nobiasC:/home/out \
        zhengh42/salmon:0.8.2 \
        salmon quant -i /home/ref/gencode.v25.salmon.idx -l A -p 4 --numBootstraps 100 \
        -o /home/out/$ida -1 /home/seq/${ida}_1_val_1.fq.gz -2 /home/seq/${ida}_2_val_2.fq.gz  \
        1> $work_dir/logs/$ida.salmon.log 2>&1

#######
# STAR
#######
sjdbOverhang=47 # read length - 1
docker run -v $work_dir/TCGA_RNASeq/seq:/home/seq -v $ref_dir:/home/ref -v $work_dir/TCGA_RNASeq/STAR:/home/out \
	zhengh42/star:2.5.3a --genomeDir /home/ref/STAR.index.sjdbOverhang${sjdbOverhang} --runThreadN 4 \
	--readFilesIn /home/seq/${ida}_1_val_1.fq.gz -2 /home/seq/${ida}_2_val_2.fq.gz --readFilesCommand zcat \
	--outFileNamePrefix /home/out/$ida. --outSAMtype BAM  SortedByCoordinate \
	--twopassMode Basic --sjdbOverhang 47 \
	--quantMode TranscriptomeSAM GeneCounts --chimSegmentMin 10 \
	--limitBAMsortRAM 56482990933 \		
	--outSAMunmapped Within  --outFilterType BySJout  --outSAMattributes NH HI AS NM MD

#######
# RSEM
#######
docker run -v $work_dir/TCGA_RNASeq/STAR:/home/bam -v $ref_dir:/home/ref -v $work_dir/TCGA_RNASeq/RSEM:/home/out  \
	zhengh42/rsem:1.3.0 rsem-calculate-expression -p 4 --forward-prob 0.5 --paired-end \
	--bam /home/bam/$ida.Aligned.toTranscriptome.out.bam --no-bam-output \
	/home/ref/RSEM/gencode.v25 /home/out/$ida --estimate-rspd --append-names

#######
# HTSeq
#######
docker run -v $work_dir/TCGA_RNASeq/STAR:/home/bam -v $ref_dir:/home/ref -v $work_dir/TCGA_RNASeq/HTSeq:/home/out  \
	zhengh42/htseq:0.7.2 htseq-count -f bam /home/bam/$ida.Aligned.sortedByCoord.out.bam \
	/home/ref/gencode.v25.annotation.gtf -s no --order pos > $work_dir/TCGA_RNASeq/HTSeq/$ida.gene.tab

#######
# featureCounts
#######
docker run -v $work_dir/TCGA_RNASeq/STAR:/home/bam -v $ref_dir:/home/ref -v $work_dir/TCGA_RNASeq/featureCounts:/home/out  \
	zhengh42/featurecounts:1.5.2 featureCounts -a /home/ref/gencode.v25.annotation.gtf -o /home/out/$ida.out \
	/home/bam/$ida.Aligned.sortedByCoord.out.bam -T 4 -p 


