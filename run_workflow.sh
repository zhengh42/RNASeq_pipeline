work_dir=/scratch/users/zhengh42/data/RNASeq/RNASeq_pipeline
kallisto=/home/users/zhengh42/tools/kallisto/kallisto_linux-v0.44.0/kallisto
salmon=/home/users/zhengh42/tools/Salmon/Salmon-0.9.1_linux_x86_64/bin/salmon
star=/home/users/zhengh42/tools/STAR/STAR-2.5.4a/bin/Linux_x86_64/STAR
Subread=/home/users/zhengh42/tools/Subread/subread-1.6.1-Linux-x86_64/bin/subread-align
rsem=/home/users/zhengh42/tools/RSEM/RSEM-1.3.0
featureCounts=/home/users/zhengh42/tools/Subread/subread-1.6.1-Linux-x86_64/bin/featureCounts
htseqcount=/home/users/zhengh42/python/pythonv/bin/htseq-count
hisat2=/home/users/zhengh42/tools/hisat2-2.1.0/hisat2

#SEQ_DIR=$work_dir/datasets/reversestranded/TCGA
#SEQ_DIR=$work_dir/datasets/reversestranded/simul
#SEQ_DIR=$work_dir/datasets/unstranded/TCGA
#SEQ_DIR=$work_dir/datasets/unstranded/simul
#SEQ_DIR=$work_dir/datasets/reversestranded/PRJEB11797
#SEQ_DIR=$work_dir/datasets/reversestranded/PRJEB11797_simul
KALLISTO_INDEX=$work_dir/reference/gencode.v27.kallisto.idx
KALLISTO_DIR=$work_dir/output/Kallisto
SALMON_INDEX=$work_dir/reference/gencode.v27.salmon.keepDuplicates.idx/ # to make the comparision consistent with Kallisto
SALMON_DIR=$work_dir/output/Salmon
STAR_INDEX=$work_dir/reference/gencode.v27.STAR.idx.sjdbOverhang100
STAR_DIR=$work_dir/output/STAR
RSEM_INDEX=$work_dir/reference/gencode.v27.rsem.idx
RSEMBOWTIE_INDEX=$work_dir/reference/gencode.v27.rsembowtie.idx
RSEM_DIR=$work_dir/output/rsem
SUBREAD_INDEX=$work_dir/reference/gencode.v27.subread.idx
SUBREAD_DIR=$work_dir/output/Subread
HISAT2_INDEX=$work_dir/reference/gencode.v27.hisat2.idx
HISAT2_DIR=$work_dir/output/HISAT2
FEATURECOUNTS_DIR=$work_dir/output/featureCounts
HTSEQ_DIR=$work_dir/output/HTSeq
GTF=$work_dir/reference/gencode.v27.annotation.gtf

FASTQ=$1
IDA=$2
IDB=$3
strand=$4
#IDAsuffix1="_1"
#IDAsuffix2="_2"
IDAsuffix1="_1_val_1"
IDAsuffix2="_2_val_2"

#---------------------------------------------------------------------------------------------------------------
#######
# unzip fastq files
#######
tar -zxf ${SEQ_DIR}/${FASTQ} --directory ${SEQ_DIR}

#######
# QC
#######
echo "Goodluck! trim_galore started on `date`"
$trim_galore -q 20 --stringency 3 --gzip --length 20 --paired  ${SEQ_DIR}/${IDA}_1.fastq.gz ${SEQ_DIR}/${IDA}_2.fastq.gz --output_dir ${SEQ_DIR}
echo "Goodluck! trim_galore finished on `date`"

#######
# delete raw fastq files
#######
rm ${SEQ_DIR}/${IDA}_1.fastq ${SEQ_DIR}/${IDA}_2.fastq
#---------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------
#######
#  Salmon
#######
echo "Goodluck! Salmon started on `date`"
$salmon quant -i ${SALMON_INDEX} -o ${SALMON_DIR}/${IDB} -1 ${SEQ_DIR}/${IDA}${IDAsuffix1}.fq.gz -2 ${SEQ_DIR}/${IDA}${IDAsuffix2}.fq.gz -p 4 -l A 
echo "Goodluck! Salmon finished on `date`"

echo "Goodluck! Salmon_bias started on `date`"
$salmon quant -i ${SALMON_INDEX} -o ${SALMON_DIR}/${IDB}.bias -1 ${SEQ_DIR}/${IDA}${IDAsuffix1}.fq.gz -2 ${SEQ_DIR}/${IDA}${IDAsuffix2}.fq.gz -p 4 -l A --seqBias --gcBias
echo "Goodluck! Salmon_bias finished on `date`"

#######
# Kallisto
#######
echo "Goodluck! Kallisto started on `date`"
if [  $strand == "unstranded" ]; then
$kallisto quant -i ${KALLISTO_INDEX} -o ${KALLISTO_DIR}/${IDB} ${SEQ_DIR}/${IDA}${IDAsuffix1}.fq.gz ${SEQ_DIR}/${IDA}${IDAsuffix2}.fq.gz -t 4
$kallisto quant -i ${KALLISTO_INDEX} -o ${KALLISTO_DIR}/${IDB}.bias ${SEQ_DIR}/${IDA}${IDAsuffix1}.fq.gz ${SEQ_DIR}/${IDA}${IDAsuffix2}.fq.gz -t 4 --bias
fi

if [  $strand == "reverse" ]; then
$kallisto quant -i ${KALLISTO_INDEX} -o ${KALLISTO_DIR}/${IDB} ${SEQ_DIR}/${IDA}${IDAsuffix1}.fq.gz ${SEQ_DIR}/${IDA}${IDAsuffix2}.fq.gz -t 4 --rf-stranded
$kallisto quant -i ${KALLISTO_INDEX} -o ${KALLISTO_DIR}/${IDB}.bias ${SEQ_DIR}/${IDA}${IDAsuffix1}.fq.gz ${SEQ_DIR}/${IDA}${IDAsuffix2}.fq.gz -t 4 --rf-stranded --bias
fi

if [  $strand == "stranded" ]; then
$kallisto quant -i ${KALLISTO_INDEX} -o ${KALLISTO_DIR}/${IDB} ${SEQ_DIR}/${IDA}${IDAsuffix1}.fq.gz ${SEQ_DIR}/${IDA}${IDAsuffix2}.fq.gz -t 4 --fr-stranded
$kallisto quant -i ${KALLISTO_INDEX} -o ${KALLISTO_DIR}/${IDB}.bias ${SEQ_DIR}/${IDA}${IDAsuffix1}.fq.gz ${SEQ_DIR}/${IDA}${IDAsuffix2}.fq.gz -t 4 --fr-stranded --bias
fi
echo "Goodluck! Kallisto finished on `date`"
#---------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------
#######
# Subread
#######
if [  $strand == "unstranded" ]; then
$Subread -i $SUBREAD_INDEX -T 4 --multiMapping -B 4 -t 0  -S ff \
        -r ${SEQ_DIR}/${IDA}${IDAsuffix1}.fq.gz -R ${SEQ_DIR}/${IDA}${IDAsuffix2}.fq.gz \
        -o $SUBREAD_DIR/${IDB}.bam
fi

if [  $strand == "reverse" ]; then
$Subread -i $SUBREAD_INDEX -T 4 --multiMapping -B 4 -t 0  -S fr \
        -r ${SEQ_DIR}/${IDA}${IDAsuffix1}.fq.gz -R ${SEQ_DIR}/${IDA}${IDAsuffix2}.fq.gz \
        -o $SUBREAD_DIR/${IDB}.bam
fi

if [  $strand == "stranded" ]; then
$Subread -i $SUBREAD_INDEX -T 4 --multiMapping -B 4 -t 0  -S rf \
        -r ${SEQ_DIR}/${IDA}${IDAsuffix1}.fq.gz -R ${SEQ_DIR}/${IDA}${IDAsuffix2}.fq.gz \
        -o $SUBREAD_DIR/${IDB}.bam
fi

########
# HISAT2
########
if [  $strand == "unstranded" ]; then
$hisat2 -x $HISAT2_INDEX -1 ${SEQ_DIR}/${IDA}${IDAsuffix1}.fq.gz -2 ${SEQ_DIR}/${IDA}${IDAsuffix2}.fq.gz -S $HISAT2_DIR/${IDB}.sam
fi

if [  $strand == "reverse" ]; then
$hisat2 -x $HISAT2_INDEX -1 ${SEQ_DIR}/${IDA}${IDAsuffix1}.fq.gz -2 ${SEQ_DIR}/${IDA}${IDAsuffix2}.fq.gz -S $HISAT2_DIR/${IDB}.sam --rna-strandness RF
fi

if [  $strand == "stranded" ]; then
$hisat2 -x $HISAT2_INDEX -1 ${SEQ_DIR}/${IDA}${IDAsuffix1}.fq.gz -2 ${SEQ_DIR}/${IDA}${IDAsuffix2}.fq.gz -S $HISAT2_DIR/${IDB}.sam --rna-strandness FR
fi

samtools view -bS $HISAT2_DIR/${IDB}.sam > $HISAT2_DIR/${IDB}.bam

rm  $HISAT2_DIR/${IDB}.sam

#######
# STAR
#######
echo "Goodluck! STAR started on `date`"
$star --genomeDir $STAR_INDEX --runThreadN 4 \
        --readFilesIn ${SEQ_DIR}/${IDA}${IDAsuffix1}.fq.gz ${SEQ_DIR}/${IDA}${IDAsuffix2}.fq.gz --readFilesCommand zcat \
        --outFileNamePrefix $STAR_DIR/${IDB}. --outSAMtype BAM  Unsorted \
       --twopassMode Basic --sjdbOverhang 100 \
	--quantMode TranscriptomeSAM GeneCounts --chimSegmentMin 10 \
       --limitBAMsortRAM 56482990933 \
       --outSAMunmapped Within  --outFilterType BySJout  --outSAMattributes NH HI AS NM MD
echo "Goodluck! STAR finished on `date`"

samtools view -bS  $STAR_DIR/${IDB}.Chimeric.out.sam > $STAR_DIR/${IDB}.Chimeric.out.bam
rm -rf $STAR_DIR/${IDB}._STARgenome/ $STAR_DIR/${IDB}._STARpass1/ $STAR_DIR/${IDB}.Chimeric.out.sam
gzip $STAR_DIR/${IDB}.Chimeric.out.junction $STAR_DIR/${IDB}.SJ.out.tab

#---------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------
#######
# RSEMAfterBowtie
########
echo "Goodluck! RSEM started on `date`"
if [  $strand == "unstranded" ]; then
$rsem/rsem-calculate-expression -p 4 --forward-prob 0.5 --paired-end \
       --paired-end ${SEQ_DIR}/${IDA}${IDAsuffix1}.fq.gz ${SEQ_DIR}/${IDA}${IDAsuffix2}.fq.gz --no-bam-output \
        $RSEMBOWTIE_INDEX $RSEM_DIR/${IDB}.bowtie --estimate-rspd --append-names
fi

if [  $strand == "reverse" ]; then
$rsem/rsem-calculate-expression -p 4 --forward-prob 0 --paired-end \
       --paired-end ${SEQ_DIR}/${IDA}${IDAsuffix1}.fq.gz ${SEQ_DIR}/${IDA}${IDAsuffix2}.fq.gz --no-bam-output \
        $RSEMBOWTIE_INDEX $RSEM_DIR/${IDB}.bowtie --estimate-rspd --append-names
fi

if [  $strand == "stranded" ]; then
$rsem/rsem-calculate-expression -p 4 --forward-prob 1 --paired-end \
       --paired-end ${SEQ_DIR}/${IDA}${IDAsuffix1}.fq.gz ${SEQ_DIR}/${IDA}${IDAsuffix2}.fq.gz --no-bam-output \
        $RSEMBOWTIE_INDEX $RSEM_DIR/${IDB}.bowtie --estimate-rspd --append-names
fi
echo "Goodluck! RSEM finished on `date`"

#######
# RSEMAfterSTAR
#######
echo "Goodluck! RSEM started on `date`"

if [  $strand == "unstranded" ]; then
$rsem/rsem-calculate-expression -p 4 --forward-prob 0.5 --paired-end \
       --bam $STAR_DIR/${IDB}.Aligned.toTranscriptome.out.bam --no-bam-output \
        $RSEM_INDEX $RSEM_DIR/${IDB}.STAR --estimate-rspd --append-names
fi

if [  $strand == "reverse" ]; then
$rsem/rsem-calculate-expression -p 4 --forward-prob 0 --paired-end \
       --bam $STAR_DIR/${IDB}.Aligned.toTranscriptome.out.bam --no-bam-output \
        $RSEM_INDEX $RSEM_DIR/${IDB}.STAR --estimate-rspd --append-names
fi

if [  $strand == "stranded" ]; then
$rsem/rsem-calculate-expression -p 4 --forward-prob 1 --paired-end \
       --bam $STAR_DIR/${IDB}.Aligned.toTranscriptome.out.bam --no-bam-output \
        $RSEM_INDEX $RSEM_DIR/${IDB}.STAR --estimate-rspd --append-names
fi

echo "Goodluck! RSEM finished on `date`"
#---------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------
#######
# simulation
#######
sample=$1
theta0=$2
N=60000000
#out_dir=$work_dir/datasets/reversestranded/simul
#out_dir=$work_dir/datasets/unstranded/simul
#out_dir=$work_dir/datasets/reversestranded/PRJEB11797_simul
RSEM_INDEX=$work_dir/reference/gencode.v27.rsem.idx
RSEM_DIR=$work_dir/output/rsem/real
estimated_model_sample=$RSEM_DIR/${sample}.STAR.stat/$sample.STAR.model
estimated_isoform_results=$RSEM_DIR/${sample}.STAR.isoforms.results

$rsem/rsem-simulate-reads $RSEM_INDEX $estimated_model_sample $estimated_isoform_results $theta0 $N $out_dir/${sample}_simul

gzip $out_dir/${sample}_simul_1.fq
gzip $out_dir/${sample}_simul_2.fq
#---------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------
#######
# featureCountsAfterSubread
#######
echo "Goodluck! featureCounts started on `date`"
if [  $strand == "unstranded" ]; then
$featureCounts -a $GTF -o $FEATURECOUNTS_DIR/${IDB}.Subread.gene.out $SUBREAD_DIR/${IDB}.bam -T 4 -p
$featureCounts -a $GTF -o $FEATURECOUNTS_DIR/${IDB}.Subread.geneMO.out $SUBREAD_DIR/${IDB}.bam -T 4 -p -M -O --fraction
fi

if [  $strand == "reverse" ]; then
$featureCounts -a $GTF -o $FEATURECOUNTS_DIR/${IDB}.Subread.gene.out $SUBREAD_DIR/${IDB}.bam -T 4 -p -s 2
$featureCounts -a $GTF -o $FEATURECOUNTS_DIR/${IDB}.Subread.geneMO.out $SUBREAD_DIR/${IDB}.bam -T 4 -p -M -O --fraction -s 2
fi

if [  $strand == "stranded" ]; then
$featureCounts -a $GTF -o $FEATURECOUNTS_DIR/${IDB}.Subread.gene.out $SUBREAD_DIR/${IDB}.bam -T 4 -p -s 1
$featureCounts -a $GTF -o $FEATURECOUNTS_DIR/${IDB}.Subread.geneMO.out $SUBREAD_DIR/${IDB}.bam -T 4 -p -M -O --fraction -s 1
fi
echo "Goodluck! featureCounts finished on `date`"

#######
# featureCountsAfterSTAR
#######
echo "Goodluck! featureCounts started on `date`"
if [  $strand == "unstranded" ]; then
$featureCounts -a $GTF -o $FEATURECOUNTS_DIR/${IDB}.STAR.gene.out $STAR_DIR/${IDB}.Aligned.out.bam -T 4 -p
$featureCounts -a $GTF -o $FEATURECOUNTS_DIR/${IDB}.STAR.geneMO.out $STAR_DIR/${IDB}.Aligned.out.bam -T 4 -p -M -O --fraction
fi

if [  $strand == "reverse" ]; then
$featureCounts -a $GTF -o $FEATURECOUNTS_DIR/${IDB}.STAR.gene.out $STAR_DIR/${IDB}.Aligned.out.bam -T 4 -p -s 2
$featureCounts -a $GTF -o $FEATURECOUNTS_DIR/${IDB}.STAR.geneMO.out $STAR_DIR/${IDB}.Aligned.out.bam -T 4 -p -M -O --fraction -s 2
fi

if [  $strand == "stranded" ]; then
$featureCounts -a $GTF -o $FEATURECOUNTS_DIR/${IDB}.STAR.gene.out $STAR_DIR/${IDB}.Aligned.out.bam -T 4 -p -s 1
$featureCounts -a $GTF -o $FEATURECOUNTS_DIR/${IDB}.STAR.geneMO.out $STAR_DIR/${IDB}.Aligned.out.bam -T 4 -p -M -O --fraction -s 1
fi
echo "Goodluck! featureCounts finished on `date`"

#######
# featureCountsAfterHISAT2
#######
echo "Goodluck! featureCounts started on `date`"
if [  $strand == "unstranded" ]; then
$featureCounts -a $GTF -o $FEATURECOUNTS_DIR/${IDB}.HISAT2.gene.out $HISAT2_DIR/${IDB}.bam -T 4 -p
$featureCounts -a $GTF -o $FEATURECOUNTS_DIR/${IDB}.HISAT2.geneMO.out $HISAT2_DIR/${IDB}.bam -T 4 -p -M -O --fraction
fi

if [  $strand == "reverse" ]; then
$featureCounts -a $GTF -o $FEATURECOUNTS_DIR/${IDB}.HISAT2.gene.out $HISAT2_DIR/${IDB}.bam -T 4 -p -s 2
$featureCounts -a $GTF -o $FEATURECOUNTS_DIR/${IDB}.HISAT2.geneMO.out $HISAT2_DIR/${IDB}.bam -T 4 -p -M -O --fraction -s 2
fi

if [  $strand == "stranded" ]; then
$featureCounts -a $GTF -o $FEATURECOUNTS_DIR/${IDB}.HISAT2.gene.out $HISAT2_DIR/${IDB}.bam -T 4 -p -s 1
$featureCounts -a $GTF -o $FEATURECOUNTS_DIR/${IDB}.HISAT2.geneMO.out $HISAT2_DIR/${IDB}.bam -T 4 -p -M -O --fraction -s 1
fi
echo "Goodluck! featureCounts finished on `date`"
#---------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------
#######
# HTSeqAfterSTAR
########
echo "Goodluck! HTSeq started on `date`"
if [  $strand == "unstranded" ]; then
$htseqcount -f bam $STAR_DIR/${IDB}.Aligned.out.bam $GTF -s no  > $HTSEQ_DIR/${IDB}.STAR.tab
$htseqcount -f bam $STAR_DIR/${IDB}.Aligned.out.bam $GTF -s no   --nonunique all > $HTSEQ_DIR/${IDB}.STAR.MO.tab
fi

if [  $strand == "reverse" ]; then
$htseqcount -f bam $STAR_DIR/${IDB}.Aligned.out.bam $GTF -s reverse  > $HTSEQ_DIR/${IDB}.STAR.tab
$htseqcount -f bam $STAR_DIR/${IDB}.Aligned.out.bam $GTF -s reverse  --nonunique all > $HTSEQ_DIR/${IDB}.STAR.MO.tab
fi

if [  $strand == "stranded" ]; then
$htseqcount -f bam $STAR_DIR/${IDB}.Aligned.out.bam $GTF -s yes  > $HTSEQ_DIR/${IDB}.STAR.tab
$htseqcount -f bam $STAR_DIR/${IDB}.Aligned.out.bam $GTF -s yes  --nonunique all > $HTSEQ_DIR/${IDB}.STAR.MO.tab
fi
echo "Goodluck! HTSeq finished on `date`"

#######
# HTSeqAfterHISAT2
########
echo "Goodluck! HTSeq started on `date`"
if [  $strand == "unstranded" ]; then
$htseqcount -f bam $HISAT2_DIR/${IDB}.bam $GTF -s no  > $HTSEQ_DIR/${IDB}.HISAT2.tab
$htseqcount -f bam $HISAT2_DIR/${IDB}.bam $GTF -s no   --nonunique all > $HTSEQ_DIR/${IDB}.HISAT2.MO.tab
fi

if [  $strand == "reverse" ]; then
$htseqcount -f bam $HISAT2_DIR/${IDB}.bam $GTF -s reverse  > $HTSEQ_DIR/${IDB}.HISAT2.tab
$htseqcount -f bam $HISAT2_DIR/${IDB}.bam $GTF -s reverse  --nonunique all > $HTSEQ_DIR/${IDB}.HISAT2.MO.tab
fi

if [  $strand == "stranded" ]; then
$htseqcount -f bam $HISAT2_DIR/${IDB}.bam $GTF -s yes  > $HTSEQ_DIR/${IDB}.HISAT2.tab
$htseqcount -f bam $HISAT2_DIR/${IDB}.bam $GTF -s yes  --nonunique all > $HTSEQ_DIR/${IDB}.HISAT2.MO.tab
fi
echo "Goodluck! HTSeq finished on `date`"

######
# HTSeqAfterSubread
#######
echo "Goodluck! HTSeq started on `date`"
if [  $strand == "unstranded" ]; then
$htseqcount -f bam $SUBREAD_DIR/${IDB}.bam $GTF -s no  > $HTSEQ_DIR/${IDB}.Subread.tab
$htseqcount -f bam $SUBREAD_DIR/${IDB}.bam $GTF -s no   --nonunique all > $HTSEQ_DIR/${IDB}.Subread.MO.tab
fi

if [  $strand == "reverse" ]; then
$htseqcount -f bam $SUBREAD_DIR/${IDB}.bam $GTF -s reverse  > $HTSEQ_DIR/${IDB}.Subread.tab
$htseqcount -f bam $SUBREAD_DIR/${IDB}.bam $GTF -s reverse  --nonunique all > $HTSEQ_DIR/${IDB}.Subread.MO.tab
fi

if [  $strand == "stranded" ]; then
$htseqcount -f bam $SUBREAD_DIR/${IDB}.bam $GTF -s yes  > $HTSEQ_DIR/${IDB}.Subread.tab
$htseqcount -f bam $SUBREAD_DIR/${IDB}.bam $GTF -s yes  --nonunique all > $HTSEQ_DIR/${IDB}.Subread.MO.tab
fi
echo "Goodluck! HTSeq finished on `date`"
#---------------------------------------------------------------------------------------------------------------
