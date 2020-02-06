#!/bin/bash
###-------------------------------------------------------------------------------------------------------------------
### set up
###-------------------------------------------------------------------------------------------------------------------

#--------------------- sample information --------------------------------------------------
organism="hg"			### "hg" for human and "mm" for mouse
transcriptome="gencodev32"	### "gencodev32" or "gencodev32noncodev5" for human, and "gencodevM24" for mouse
paired="yes"			### "yes" for paired-end and "no" for single-end
readlength=100			### for choosing which STAR index to use. Possible values: 50, 100, 150

workDir=/srv/gevaertlab/data/Hong/RNASeq
projectid=rnaseqproject
sampleid=$1

#-------------------- analysis to run ------------------------------------------------------
run_trim_galore=0	### 1 means running the step
run_strandInfo=0
run_kallisto=1
run_star=0
run_rseqc=0

#-------------------- tools ----------------------------------------------------------------
toolsDir=/srv/gevaertlab/tools
trim_galore=$toolsDir/trim_galore/TrimGalore-0.4.4/trim_galore # cutadapt also needs to be installed and useable in the working environment
star=$toolsDir/STAR/STAR-2.7.3a/bin/Linux_x86_64/STAR
kallisto=$toolsDir/kallisto/kallisto_linux-v0.46.0/kallisto
rsem=$toolsDir/RSEM/RSEM-1.3.2

#-------------------- references -----------------------------------------------------------
refDir=/srv/gevaertlab/reference

if [ $organism == "hg" ]; then
        refgenome=$refDir/homosapiens/GRCh38.primary_assembly.genome.fa
        organismname="homosapiens"
fi

if [ $organism == "mm" ]; then
        refgenome=$refDir/musmusculus/GRCm38.primary_assembly.genome.fa
	organismname="musmusculus"
fi

KALLISTO_INDEX=$refDir/$organismname/kallisto.${transcriptome}.idx

###-------------------------------------------------------------------------------------------------------------------
### QC (basic)
###-------------------------------------------------------------------------------------------------------------------

if [ $run_trim_galore == 1 ]; then
	echo "Goodluck! trim_galore started on `date`"

	# get suffix of fastq files
	fqsuffix=`ls $workDir/$projectid/seq/${sampleid}* | sed 's/.*fq.gz/fq.gz/' | sed 's/.*fastq.gz/fastq.gz/' | head -n1`

	if [ ! -d $workDir/$projectid/seq ]; then mkdir -p $workDir/$projectid/seq; fi
	if [  $paired == "yes" ]; then
		$trim_galore -q 20 --stringency 1 --gzip --length 20 --fastqc --paired \
			$workDir/$projectid/seq/${sampleid}_1.$fqsuffix $workDir/$projectid/seq/${sampleid}_2.$fqsuffix --output_dir $workDir/$projectid/seq
	fi

	if [  $paired == "no" ]; then
		$trim_galore -q 20 --stringency 1 --gzip --length 20 --fastqc $workDir/$projectid/seq/${sampleid}.$fqsuffix --output_dir $workDir/$projectid/seq
	fi
	echo "Goodluck! trim_galore finished on `date`"
fi

###-------------------------------------------------------------------------------------------------------------------
### infer strandness
###-------------------------------------------------------------------------------------------------------------------

if [ $run_strandInfo == 1 ]; then
	if [ ! -d $workDir/$projectid/kallisto ]; then mkdir -p $workDir/$projectid/kallisto;fi
        echo "Goodluck! Kallisto started on `date`"
	
	if [  $paired == "yes" ]; then
		zcat $workDir/$projectid/seq/${sampleid}_1_val_1.fq.gz | head -n 40000 > $workDir/$projectid/seq/${sampleid}_1.test.fq
		zcat $workDir/$projectid/seq/${sampleid}_2_val_2.fq.gz | head -n 40000 > $workDir/$projectid/seq/${sampleid}_2.test.fq

		$kallisto quant -i $KALLISTO_INDEX -o $workDir/$projectid/kallisto/${sampleid}.un $workDir/$projectid/seq/${sampleid}_1.test.fq $workDir/$projectid/seq/${sampleid}_2.test.fq 
		$kallisto quant -i $KALLISTO_INDEX -o $workDir/$projectid/kallisto/${sampleid}.rf $workDir/$projectid/seq/${sampleid}_1.test.fq $workDir/$projectid/seq/${sampleid}_2.test.fq --rf-stranded
		$kallisto quant -i $KALLISTO_INDEX -o $workDir/$projectid/kallisto/${sampleid}.fr $workDir/$projectid/seq/${sampleid}_1.test.fq $workDir/$projectid/seq/${sampleid}_2.test.fq --fr-stranded
	
		paste $workDir/$projectid/kallisto/${sampleid}.un/abundance.tsv $workDir/$projectid/kallisto/${sampleid}.rf/abundance.tsv $workDir/$projectid/kallisto/${sampleid}.fr/abundance.tsv | cut -f1,4,9,14  | awk 'BEGIN{sum1=0;sum2=0;sun3=0}{sum1+=$2;sum2+=$3;sum3+=$4}END{print sum1,sum2,sum3}' > $workDir/$projectid/kallisto/${sampleid}.libtypetesting
		less $workDir/$projectid/kallisto/${sampleid}.libtypetesting | awk '{print $2/$1,$3/$1,$3/$2}' | awk '{if($1<0.3 && $3>3)print "stranded";else if($1>3 && $2>3)print "reverse";else print "unstranded"}' > $workDir/$projectid/kallisto/${sampleid}.libtype
		rm $workDir/$projectid/seq/${sampleid}_1.test.fq $workDir/$projectid/seq/${sampleid}_2.test.fq
		rm -rf $workDir/$projectid/kallisto/${sampleid}.un $workDir/$projectid/kallisto/${sampleid}.rf $workDir/$projectid/kallisto/${sampleid}.fr
	fi

	if [  $paired == "no" ]; then
                zcat $workDir/$projectid/seq/${sampleid}_trimmed.fq.gz | head -n 40000 > $workDir/$projectid/seq/${sampleid}_trimmed.test.fq

                $kallisto quant -i $KALLISTO_INDEX -o $workDir/$projectid/kallisto/${sampleid}.un $workDir/$projectid/seq/${sampleid}_trimmed.test.fq --single -l 190 -s 20
                $kallisto quant -i $KALLISTO_INDEX -o $workDir/$projectid/kallisto/${sampleid}.rf $workDir/$projectid/seq/${sampleid}_trimmed.test.fq --single -l 190 -s 20 --rf-stranded
                $kallisto quant -i $KALLISTO_INDEX -o $workDir/$projectid/kallisto/${sampleid}.fr $workDir/$projectid/seq/${sampleid}_trimmed.test.fq --single -l 190 -s 20 --fr-stranded

                paste $workDir/$projectid/kallisto/${sampleid}.un/abundance.tsv $workDir/$projectid/kallisto/${sampleid}.rf/abundance.tsv $workDir/$projectid/kallisto/${sampleid}.fr/abundance.tsv | cut -f1,4,9,14  | awk 'BEGIN{sum1=0;sum2=0;sun3=0}{sum1+=$2;sum2+=$3;sum3+=$4}END{print sum1,sum2,sum3}' > $workDir/$projectid/kallisto/${sampleid}.libtypetesting
                less $workDir/$projectid/kallisto/${sampleid}.libtypetesting | awk '{print $2/$1,$3/$1,$3/$2}' | awk '{if($1<0.3 && $3>3)print "stranded";else if($1>3 && $2>3)print "reverse";else print "unstranded"}' > $workDir/$projectid/kallisto/${sampleid}.libtype
                rm $workDir/$projectid/seq/${sampleid}_trim.test.fq 
                rm -rf $workDir/$projectid/kallisto/${sampleid}.un $workDir/$projectid/kallisto/${sampleid}.rf $workDir/$projectid/kallisto/${sampleid}.fr
	fi

	strand=`less $workDir/$projectid/kallisto/${sampleid}.libtype | awk '{print $1}'`
	echo "The sample is $strand"
	echo "Goodluck! Kallisto finished on `date`"

fi

###-------------------------------------------------------------------------------------------------------------------
### kallisto
###-------------------------------------------------------------------------------------------------------------------

if [ $run_kallisto == 1 ]; then
        if [ ! -d $workDir/$projectid/kallisto ]; then mkdir -p $workDir/$projectid/kallisto;fi
        echo "Goodluck! Kallisto started on `date`"
	strand=`less $workDir/$projectid/kallisto/${sampleid}.libtype | awk '{print $1}'`

        if [  $paired == "yes" ]; then
		 if [  $strand == "unstranded" ]; then
			$kallisto quant -i $KALLISTO_INDEX -o $workDir/$projectid/kallisto/${sampleid}.$transcriptome $workDir/$projectid/seq/${sampleid}_1_val_1.fq.gz $workDir/$projectid/seq/${sampleid}_2_val_2.fq.gz 
		fi

		if [  $strand == "reverse" ]; then
			$kallisto quant -i $KALLISTO_INDEX -o $workDir/$projectid/kallisto/${sampleid}.$transcriptome $workDir/$projectid/seq/${sampleid}_1_val_1.fq.gz $workDir/$projectid/seq/${sampleid}_2_val_2.fq.gz --rf-stranded
		fi

		if [  $strand == "stranded" ]; then
			kallisto quant -i $KALLISTO_INDEX -o $workDir/$projectid/kallisto/${sampleid}.$transcriptome $workDir/$projectid/seq/${sampleid}_1_val_1.fq.gz $workDir/$projectid/seq/${sampleid}_2_val_2.fq.gz --fr-stranded
		fi
	fi

	if [  $paired == "no" ]; then
                 if [  $strand == "unstranded" ]; then
                        $kallisto quant -i $KALLISTO_INDEX -o $workDir/$projectid/kallisto/${sampleid}.$transcriptome $workDir/$projectid/seq/${sampleid}_trimmed.fq.gz --single -l 190 -s 20
                fi

                if [  $strand == "reverse" ]; then
                        $kallisto quant -i $KALLISTO_INDEX -o $workDir/$projectid/kallisto/${sampleid}.$transcriptome $workDir/$projectid/seq/${sampleid}_trimmed.fq.gz --single -l 190 -s 20 --rf-stranded
                fi

                if [  $strand == "stranded" ]; then
                        $kallisto quant -i $KALLISTO_INDEX -o $workDir/$projectid/kallisto/${sampleid}.$transcriptome $workDir/$projectid/seq/${sampleid}_trimmed.fq.gz --single -l 190 -s 20 --fr-stranded
                fi
        fi

	gzip $workDir/$projectid/kallisto/${sampleid}/abundance.tsv
	echo "Goodluck! Kallisto finished on `date`"
fi


###-------------------------------------------------------------------------------------------------------------------
### STAR
###-------------------------------------------------------------------------------------------------------------------

if [ $run_star == 1 ]; then
	if [ ! -d $workDir/$projectid/star ]; then mkdir -p $workDir/$projectid/star;fi
        echo "Goodluck! STAR started on `date`"

        if [  $paired == "yes" ]; then
		$star --genomeDir $refDir/$organismname/star.${transcriptome}.sjdbOverhang${readlength}.idx --runThreadN 2 \
               --readFilesIn $workDir/$projectid/seq/${sampleid}_1_val_1.fq.gz $workDir/$projectid/seq/${sampleid}_2_val_2.fq.gz --readFilesCommand zcat \
               --outFileNamePrefix $workDir/$projectid/star/${sampleid}.${transcriptome}. --outSAMtype BAM  Unsorted \
               --twopassMode Basic --sjdbOverhang $readlength \
               --quantMode TranscriptomeSAM GeneCounts \
               --limitBAMsortRAM 56482990933 \
               --outFilterType BySJout \
               --outSAMunmapped Within  \
               --outSAMattributes NH HI AS nM NM MD XS \
               --outSAMstrandField intronMotif \
               --chimSegmentMin 12 \
               --chimJunctionOverhangMin 12 \
               --chimOutJunctionFormat 1 \
               --alignSJDBoverhangMin 10 \
               --alignMatesGapMax 100000 \
               --alignIntronMax 100000 \
               --alignSJstitchMismatchNmax 5 -1 5 5 \
               --outSAMattrRGline ID:GRPundef \
               --chimMultimapScoreRange 3 \
               --chimScoreJunctionNonGTAG -4 \
               --chimMultimapNmax 20 \
               --chimNonchimScoreDropMin 10 \
               --peOverlapNbasesMin 12 \
               --peOverlapMMp 0.1
	fi

	if [  $paired == "no" ]; then
		$star --genomeDir $refDir/$organismname/star.${transcriptome}.sjdbOverhang${readlength}.idx --runThreadN 2 \
               --readFilesIn $workDir/$projectid/seq/${sampleid}_trimmed.fq.gz  --readFilesCommand zcat \
               --outFileNamePrefix $workDir/$projectid/star/${sampleid}.${transcriptome} --outSAMtype BAM  Unsorted \
               --twopassMode Basic --sjdbOverhang $readlength \
               --quantMode TranscriptomeSAM GeneCounts \
               --limitBAMsortRAM 56482990933 \
               --outFilterType BySJout \
               --outSAMunmapped Within  \
               --outSAMattributes NH HI AS nM NM MD XS \
               --outSAMstrandField intronMotif \
               --chimSegmentMin 12 \
               --chimJunctionOverhangMin 12 \
               --chimOutJunctionFormat 1 \
               --alignSJDBoverhangMin 10 \
               --alignMatesGapMax 100000 \
               --alignIntronMax 100000 \
               --alignSJstitchMismatchNmax 5 -1 5 5 \
               --outSAMattrRGline ID:GRPundef \
               --chimMultimapScoreRange 3 \
               --chimScoreJunctionNonGTAG -4 \
               --chimMultimapNmax 20 \
               --chimNonchimScoreDropMin 10 \
               --peOverlapNbasesMin 12 \
               --peOverlapMMp 0.1
        fi

	echo "Goodluck! STAR finished on `date`"
fi

###-------------------------------------------------------------------------------------------------------------------
### QC (advanced)
###-------------------------------------------------------------------------------------------------------------------
if [ $run_rseqc == 1 ]; then
	genebed=$refDir/$organismname/gencodev${transcriptome}.bed
	if [ ! -d $workDir/$projectid/rseqc ]; then mkdir -p $workDir/$projectid/rseqc;fi
	bam_stat.py  -i $workDir/$projectid/star/${sampleid}.${transcriptome}.Aligned.out.bam  > $workDir/$projectid/rseqc/${sampleid}.bam_stat
	read_distribution.py  -i $workDir/$projectid/star/${sampleid}.${transcriptome}.Aligned.out.bam -r $genebed  > $workDir/$projectid/rseqc/${sampleid}.read_distribution
	inner_distance.py -i $workDir/$projectid/star/${sampleid}.${transcriptome}.Aligned.out.bam -r $genebed -o $workDir/$projectid/rseqc/${sampleid}
	junction_annotation.py -i $workDir/$projectid/star/${sampleid}.${transcriptome}.Aligned.out.bam -r $genebed -o $workDir/$projectid/rseqc/${sampleid}.junction_annotation
	junction_saturation.py -i $workDir/$projectid/star/${sampleid}.${transcriptome}.Aligned.out.bam -r $genebed -o $workDir/$projectid/rseqc/${sampleid}.junction_saturation
	read_duplication.py -i $workDir/$projectid/star/${sampleid}.${transcriptome}.Aligned.out.bam -o $workDir/$projectid/rseqc/${sampleid}.read_duplication
fi
