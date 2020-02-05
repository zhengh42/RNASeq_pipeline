###-------------------------------------------------------------------------------------------------------------------
### SET UP
### refDir: the directory for building and storing the references and indexes
### organism: "hg" for homo sapiens and "mm" for mus musculus
### transcriptome: "gencode" for GENCODE, or "gencode_noncode" for GENCODE and NONCODE combined. "gencode_noncode" option is only valid for homo sapiens at the moment.
### gencodeversion: for example 32 or homo sapiens and M24 for mus musculus
### readlength: for choosing which STAR index to use (setting -sjdbOverhang)
###-------------------------------------------------------------------------------------------------------------------
refDir=/srv/gevaertlab/reference
if [ ! -d $refDir/homosapiens ]; then mkdir -p $workDir/homosapiens;fi
if [ ! -d $refDir/musmusculus ]; then mkdir -p $workDir/musmusculus;fi

organism="hg"
transcriptome="gencode_noncode"
gencodeversion=32
readlength=75

# set the number to 1 if 
gencodedownload=0
noncodedownload=0
buildKALLISTOindex=0
buildSTARindex=0
buildRSEMindex=0
###-------------------------------------------------------------------------------------------------------------------
###-------------------------------------------------------------------------------------------------------------------


###-------------------------------------------------------------------------------------------------------------------
### TOOLS
### tools to be installed by the user: STAR, Kallisto, RSEM
### Other tools can be found in smalltools/
###-------------------------------------------------------------------------------------------------------------------
toolsDir=/srv/gevaertlab/tools
star=$toolsDir/STAR/STAR-2.7.3a/bin/Linux_x86_64/STAR
kallisto=$toolsDir/kallisto/kallisto_linux-v0.46.0/kallisto
rsem=$toolsDir/RSEM/RSEM-1.3.2
gffread=$toolsDir/smalltools/gffread/gffread
gtfToGenePred=$toolsDir/smalltools/gtfToGenePred
gtf2bed=$toolsDir/smalltools/gtf2bed
liftOver=$toolsDir/smalltools/liftOver
hg19ToHg38overchain=$toolsDir/smalltools/hg19ToHg38.over.chain
sequence_cleaner=$toolsDir/smalltools/sequence_cleaner.py
Exclude=$toolsDir/smalltools/Exclude.pl
###-------------------------------------------------------------------------------------------------------------------
###-------------------------------------------------------------------------------------------------------------------


###-------------------------------------------------------------------------------------------------------------------
### genome
### human: wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh38.primary_assembly.genome.fa.gz
### mice: wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M24/GRCm38.primary_assembly.genome.fa.gz
###-------------------------------------------------------------------------------------------------------------------
if [[ $organism == "hg" ]]; then
	organismname="homosapiens"
	refgenome=$refDir/$organismname/GRCh38.primary_assembly.genome.fa
fi

if [[ $organism == "mm" ]]; then
	organismname="musmusculus"
        refgenome=$refDir/musmusculus/GRCm38.primary_assembly.genome.fa
fi
###-------------------------------------------------------------------------------------------------------------------
###-------------------------------------------------------------------------------------------------------------------


###-------------------------------------------------------------------------------------------------------------------
### GENCODE download
###-------------------------------------------------------------------------------------------------------------------
if [[ $gencodedownload == 1 ]]; then
	if [[ $organism == "hg" ]]; then
                wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${gencodeversion}/gencode.v${gencodeversion}.annotation.gtf.gz
	fi

	if [[ $organism == "mm" ]]; then
                wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_${gencodeversion}/gencode.v${gencodeversion}.annotation.gtf.gz
	fi

	mv gencode.v${gencodeversion}.*.gz $refDir/$organismname/
	gzip -d $refDir/$organismname/gencode.v${gencodeversion}*.gz
	mv $refDir/$organismname/gencode.v${gencodeversion}.annotation.gtf $refDir/$organismname/gencodev${gencodeversion}.annotation.gtf
	$gffread $refDir/$organismname/gencodev${gencodeversion}.annotation.gtf -g $refgenome -w $refDir/$organismname/gencodev${gencodeversion}.fa
	less $refDir/$organismname/gencodev${gencodeversion}.annotation.gtf | egrep -v '^#' | awk '$3~/gene/' | cut -f9 | awk '{print $2,$6}' | tr -d '";' | sed 's/\s\+/\t/' > $refDir/$organismname/gencodev${gencodeversion}.gene_id_name
	$gtf2bed $refDir/$organismname/gencodev${gencodeversion}.annotation.gtf > $refDir/$organismname/gencodev${gencodeversion}.bed

fi
###-------------------------------------------------------------------------------------------------------------------
###-------------------------------------------------------------------------------------------------------------------


###-------------------------------------------------------------------------------------------------------------------
### NONCODE download and merge with GENCODE
###-------------------------------------------------------------------------------------------------------------------
if [[ $noncodedownload == 1 ]]; then
	wget http://www.noncode.org/datadownload/NONCODEv5_human_hg38_lncRNA.gtf.gz
	mv NONCODEv5_human.* $refDir/$organismname/
	gzip -d NONCODEv5_human*.gz

	# change chromosome ids to keep them consistent with the referece genome file
	sed 's/chr4_GL000008v2_random/GL000008.2/;s/chr14_GL000009v2_random/GL000009.2/;s/chr14_GL000194v1_random/GL000194.1/;s/chrUn_GL000195v1/GL000195.1/;s/chr17_GL000205v2_random/GL000205.2/;s/chr5_GL000208v1_random/GL000208.1/;s/chrUn_GL000213v1/GL000213.1/;s/chrUn_GL000214v1/GL000214.1/;s/chrUn_GL000216v2/GL000216.2/;s/chrUn_GL000218v1/GL000218.1/;s/chrUn_GL000219v1/GL000219.1/;s/chrUn_GL000220v1/GL000220.1/;s/chr3_GL000221v1_random/GL000221.1/;s/chrUn_GL000224v1/GL000224.1/;s/chr14_GL000225v1_random/GL000225.1/;s/chrUn_GL000226v1/GL000226.1/;s/chrUn_KI270302v1/KI270302.1/;s/chrUn_KI270303v1/KI270303.1/;s/chrUn_KI270304v1/KI270304.1/;s/chrUn_KI270305v1/KI270305.1/;s/chrUn_KI270310v1/KI270310.1/;s/chrUn_KI270311v1/KI270311.1/;s/chrUn_KI270312v1/KI270312.1/;s/chrUn_KI270315v1/KI270315.1/;s/chrUn_KI270316v1/KI270316.1/;s/chrUn_KI270317v1/KI270317.1/;s/chrUn_KI270320v1/KI270320.1/;s/chrUn_KI270322v1/KI270322.1/;s/chrUn_KI270329v1/KI270329.1/;s/chrUn_KI270330v1/KI270330.1/;s/chrUn_KI270333v1/KI270333.1/;s/chrUn_KI270334v1/KI270334.1/;s/chrUn_KI270335v1/KI270335.1/;s/chrUn_KI270336v1/KI270336.1/;s/chrUn_KI270337v1/KI270337.1/;s/chrUn_KI270338v1/KI270338.1/;s/chrUn_KI270340v1/KI270340.1/;s/chrUn_KI270362v1/KI270362.1/;s/chrUn_KI270363v1/KI270363.1/;s/chrUn_KI270364v1/KI270364.1/;s/chrUn_KI270366v1/KI270366.1/;s/chrUn_KI270371v1/KI270371.1/;s/chrUn_KI270372v1/KI270372.1/;s/chrUn_KI270373v1/KI270373.1/;s/chrUn_KI270374v1/KI270374.1/;s/chrUn_KI270375v1/KI270375.1/;s/chrUn_KI270376v1/KI270376.1/;s/chrUn_KI270378v1/KI270378.1/;s/chrUn_KI270379v1/KI270379.1/;s/chrUn_KI270381v1/KI270381.1/;s/chrUn_KI270382v1/KI270382.1/;s/chrUn_KI270383v1/KI270383.1/;s/chrUn_KI270384v1/KI270384.1/;s/chrUn_KI270385v1/KI270385.1/;s/chrUn_KI270386v1/KI270386.1/;s/chrUn_KI270387v1/KI270387.1/;s/chrUn_KI270388v1/KI270388.1/;s/chrUn_KI270389v1/KI270389.1/;s/chrUn_KI270390v1/KI270390.1/;s/chrUn_KI270391v1/KI270391.1/;s/chrUn_KI270392v1/KI270392.1/;s/chrUn_KI270393v1/KI270393.1/;s/chrUn_KI270394v1/KI270394.1/;s/chrUn_KI270395v1/KI270395.1/;s/chrUn_KI270396v1/KI270396.1/;s/chrUn_KI270411v1/KI270411.1/;s/chrUn_KI270412v1/KI270412.1/;s/chrUn_KI270414v1/KI270414.1/;s/chrUn_KI270417v1/KI270417.1/;s/chrUn_KI270418v1/KI270418.1/;s/chrUn_KI270419v1/KI270419.1/;s/chrUn_KI270420v1/KI270420.1/;s/chrUn_KI270422v1/KI270422.1/;s/chrUn_KI270423v1/KI270423.1/;s/chrUn_KI270424v1/KI270424.1/;s/chrUn_KI270425v1/KI270425.1/;s/chrUn_KI270429v1/KI270429.1/;s/chrUn_KI270435v1/KI270435.1/;s/chrUn_KI270438v1/KI270438.1/;s/chrUn_KI270442v1/KI270442.1/;s/chrUn_KI270448v1/KI270448.1/;s/chrUn_KI270465v1/KI270465.1/;s/chrUn_KI270466v1/KI270466.1/;s/chrUn_KI270467v1/KI270467.1/;s/chrUn_KI270468v1/KI270468.1/;s/chrUn_KI270507v1/KI270507.1/;s/chrUn_KI270508v1/KI270508.1/;s/chrUn_KI270509v1/KI270509.1/;s/chrUn_KI270510v1/KI270510.1/;s/chrUn_KI270511v1/KI270511.1/;s/chrUn_KI270512v1/KI270512.1/;s/chrUn_KI270515v1/KI270515.1/;s/chrUn_KI270516v1/KI270516.1/;s/chrUn_KI270517v1/KI270517.1/;s/chrUn_KI270518v1/KI270518.1/;s/chrUn_KI270519v1/KI270519.1/;s/chrUn_KI270521v1/KI270521.1/;s/chrUn_KI270522v1/KI270522.1/;s/chrUn_KI270528v1/KI270528.1/;s/chrUn_KI270529v1/KI270529.1/;s/chrUn_KI270530v1/KI270530.1/;s/chrUn_KI270538v1/KI270538.1/;s/chrUn_KI270539v1/KI270539.1/;s/chrUn_KI270544v1/KI270544.1/;s/chrUn_KI270548v1/KI270548.1/;s/chrUn_KI270579v1/KI270579.1/;s/chrUn_KI270580v1/KI270580.1/;s/chrUn_KI270581v1/KI270581.1/;s/chrUn_KI270582v1/KI270582.1/;s/chrUn_KI270583v1/KI270583.1/;s/chrUn_KI270584v1/KI270584.1/;s/chrUn_KI270587v1/KI270587.1/;s/chrUn_KI270588v1/KI270588.1/;s/chrUn_KI270589v1/KI270589.1/;s/chrUn_KI270590v1/KI270590.1/;s/chrUn_KI270591v1/KI270591.1/;s/chrUn_KI270593v1/KI270593.1/;s/chr1_KI270706v1_random/KI270706.1/;s/chr1_KI270707v1_random/KI270707.1/;s/chr1_KI270708v1_random/KI270708.1/;s/chr1_KI270709v1_random/KI270709.1/;s/chr1_KI270710v1_random/KI270710.1/;s/chr1_KI270711v1_random/KI270711.1/;s/chr1_KI270712v1_random/KI270712.1/;s/chr1_KI270713v1_random/KI270713.1/;s/chr1_KI270714v1_random/KI270714.1/;s/chr2_KI270715v1_random/KI270715.1/;s/chr2_KI270716v1_random/KI270716.1/;s/chr9_KI270717v1_random/KI270717.1/;s/chr9_KI270718v1_random/KI270718.1/;s/chr9_KI270719v1_random/KI270719.1/;s/chr9_KI270720v1_random/KI270720.1/;s/chr11_KI270721v1_random/KI270721.1/;s/chr14_KI270722v1_random/KI270722.1/;s/chr14_KI270723v1_random/KI270723.1/;s/chr14_KI270724v1_random/KI270724.1/;s/chr14_KI270725v1_random/KI270725.1/;s/chr14_KI270726v1_random/KI270726.1/;s/chr15_KI270727v1_random/KI270727.1/;s/chr16_KI270728v1_random/KI270728.1/;s/chr17_KI270729v1_random/KI270729.1/;s/chr17_KI270730v1_random/KI270730.1/;s/chr22_KI270731v1_random/KI270731.1/;s/chr22_KI270732v1_random/KI270732.1/;s/chr22_KI270733v1_random/KI270733.1/;s/chr22_KI270734v1_random/KI270734.1/;s/chr22_KI270735v1_random/KI270735.1/;s/chr22_KI270736v1_random/KI270736.1/;s/chr22_KI270737v1_random/KI270737.1/;s/chr22_KI270738v1_random/KI270738.1/;s/chr22_KI270739v1_random/KI270739.1/;s/chrY_KI270740v1_random/KI270740.1/;s/chrUn_KI270741v1/KI270741.1/;s/chrUn_KI270742v1/KI270742.1/;s/chrUn_KI270743v1/KI270743.1/;s/chrUn_KI270744v1/KI270744.1/;s/chrUn_KI270745v1/KI270745.1/;s/chrUn_KI270746v1/KI270746.1/;s/chrUn_KI270747v1/KI270747.1/;s/chrUn_KI270748v1/KI270748.1/;s/chrUn_KI270749v1/KI270749.1/;s/chrUn_KI270750v1/KI270750.1/;s/chrUn_KI270751v1/KI270751.1/;s/chrUn_KI270752v1/KI270752.1/;s/chrUn_KI270753v1/KI270753.1/;s/chrUn_KI270754v1/KI270754.1/;s/chrUn_KI270755v1/KI270755.1/;s/chrUn_KI270756v1/KI270756.1/;s/chrUn_KI270757v1/KI270757.1/' $refDir/$organismname/NONCODEv5_human_hg38_lncRNA.gtf > $refDir/$organismname/NONCODEv5_human_hg38_lncRNA.modifieda.gtf
	
	# merge GENCODE and NONCODE transcriptome fasta and remove redundant sequences
	less $refDir/$organismname/NONCODEv5_human_hg38_lncRNA.modifieda.gtf | awk '$1!~/alt/' | awk -F "\t" 'OFS="\t"{if($7==".")$7="+";print $0}' > $refDir/$organismname/NONCODEv5_human_hg38_lncRNA.modifiedb.gtf
	$gffread $refDir/$organismname/NONCODEv5_human_hg38_lncRNA.modifiedb.gtf  -g $refgenome -w $refDir/$organismname/NONCODEv5_human_hg38_lncRNA.gffread.fa
	cat $refDir/$organismname/gencodev${gencodeversion}.fa $refDir/$organismname/NONCODEv5_human_hg38_lncRNA.gffread.fa > $refDir/$organismname/gencodev${gencodeversion}noncodev5.transcripts.fa
	python $sequence_cleaner $refDir/$organismname/gencodev${gencodeversion}noncodev5.transcripts.fa 0 100 > $refDir/$organismname/gencodev${gencodeversion}noncodev5.transcripts.clean.fa

	less $refDir/$organismname/gencodev${gencodeversion}noncodev5.transcripts.clean.fa  | egrep '^>' | sed 's/>//' > $refDir/$organismname/gencodev${gencodeversion}noncodev5.transcripts.clean.ids
	less $refDir/$organismname/gencodev${gencodeversion}noncodev5.transcripts.clean.ids | egrep '_' | cut -f2- -d '_' | sed 's/_NONHSAT/\nNONHSAT/g;s/_ENST/\nENST/g' > $refDir/$organismname/gencodev${gencodeversion}noncodev5.transcripts.clean.ids2remove
	less $refDir/$organismname/gencodev${gencodeversion}noncodev5.transcripts.clean.fa | sed 's/_NONHSAT.*//;s/_ENST.*//' > $refDir/$organismname/gencodev${gencodeversion}noncodev5.fa
	
	# generate GTF file for merge GENCODE+NONCODE
	less $refDir/$organismname/gencodev${gencodeversion}.annotation.gtf | egrep -v ^# | awk '$3!~/gene|CDS|UTR|start_codon|stop_codon|Selenocysteine/' > $refDir/$organismname/gencodev${gencodeversion}.annotation.gtf.tmp1
	cat $refDir/$organismname/gencodev${gencodeversion}.annotation.gtf.tmp1 $refDir/$organismname/NONCODEv5_human_hg38_lncRNA.modifiedb.gtf > $refDir/$organismname/gencodev${gencodeversion}noncodev5.gtf
	less $refDir/$organismname/gencodev${gencodeversion}noncodev5.gtf | awk 'OFS="\t"{print $12,$0}' | sed 's/"//;s/";//' | $Exclude $refDir/$organismname/gencodev${gencodeversion}noncodev5.transcripts.clean.ids2remove 1 - 1 | cut -f2- > $refDir/$organismname/gencodev${gencodeversion}noncodev5.annotation.gtf
	
	# get the gene type information
	less $refDir/$organismname/gencodev${gencodeversion}noncodev5.annotation.gtf | awk '$3=="transcript"' | grep 'ENST' | awk '{print $10,$14}' | sort | uniq | sed 's/[";]//g'  | sed 's/\s\+/\t/' > $refDir/$organismname/gencodev${gencodeversion}noncodev5.genetype1
	less $refDir/$organismname/gencodev${gencodeversion}noncodev5.annotation.gtf | awk '$3=="transcript"' | grep -v 'ENST' | awk '{print $10}' | sort | uniq | sed 's/[";]//g'  | awk 'OFS="\t"{print $1,"lncRNA"}' > $refDir/$organismname/gencodev${gencodeversion}noncodev5.genetype2
	cat $refDir/$organismname/gencodev${gencodeversion}noncodev5.genetype1 $refDir/$organismname/gencodev${gencodeversion}noncodev5.genetype2 > $refDir/$organismname/gencodev${gencodeversion}noncodev5.genetype

fi
###-------------------------------------------------------------------------------------------------------------------
###-------------------------------------------------------------------------------------------------------------------


###-------------------------------------------------------------------------------------------------------------------
### Build Kallisto index
###-------------------------------------------------------------------------------------------------------------------
if [[ $buildKALLISTOindex == 1 && $transcriptome == "gencode" ]]; then
	$kallisto index -i $refDir/$organismname/kallisto.gencodev${gencodeversion}.idx  $refDir/$organismname/gencodev${gencodeversion}.fa
fi

if [[ $buildKALLISTOindex == 1 && $transcriptome == "gencode_noncode" ]]; then
        $kallisto index -i $refDir/$organismname/kallisto.gencodev${gencodeversion}noncodev5.idx  $refDir/$organismname/gencodev${gencodeversion}noncodev5.fa
fi
###-------------------------------------------------------------------------------------------------------------------
###-------------------------------------------------------------------------------------------------------------------


###-------------------------------------------------------------------------------------------------------------------
### Build STAR index
###-------------------------------------------------------------------------------------------------------------------
if [[ $buildSTARindex == 1 && $transcriptome == "gencode" ]]; then
	sjdbOverhang=$readlength
	mkdir -p $refDir/$organismname/star.gencodev${gencodeversion}.sjdbOverhang${sjdbOverhang}.idx
	$star --runThreadN 4 --runMode genomeGenerate \
       		--genomeDir $refDir/$organismname/star.gencodev${gencodeversion}.sjdbOverhang${sjdbOverhang}.idx \
       		--genomeFastaFiles $refgenome \
       		--sjdbGTFfile $refDir/$organismname/gencodev${gencodeversion}.annotation.gtf \
       		--sjdbOverhang $sjdbOverhang
fi

if [[ $buildSTARindex == 1 && $transcriptome == "gencode_noncode" ]]; then
        sjdbOverhang=$readlength
        mkdir -p $refDir/$organismname/star.gencodev${gencodeversion}noncodev5.sjdbOverhang${sjdbOverhang}.idx
        $star --runThreadN 4 --runMode genomeGenerate \
                --genomeDir $refDir/$organismname/star.gencodev${gencodeversion}noncodev5.sjdbOverhang${sjdbOverhang}.idx \
                --genomeFastaFiles $refgenome \
                --sjdbGTFfile $refDir/$organismname/gencodev${gencodeversion}noncodev5.annotation.gtf \
                --sjdbOverhang $sjdbOverhang
fi
###-------------------------------------------------------------------------------------------------------------------
###-------------------------------------------------------------------------------------------------------------------


###-------------------------------------------------------------------------------------------------------------------
### Build RSEM index
###-------------------------------------------------------------------------------------------------------------------
if [[ $buildRSEMindex == 1 && $transcriptome == "gencode" ]]; then
	$rsem/rsem-prepare-reference -p 4 --gtf $refDir/$organismname/gencodev${gencodeversion}.annotation.gtf $refgenome $refDir/$organismname/rsem.gencodev${gencodeversion}.idx
fi

if [[ $buildRSEMindex == 1 && $transcriptome == "gencode_noncode" ]]; then
	$rsem/rsem-prepare-reference -p 4 --gtf $refDir/$organismname/gencodev${gencodeversion}noncodev5.annotation.gtf $refgenome $refDir/$organismname/rsem.gencodev${gencodeversion}noncodev5.idx
fi
###-------------------------------------------------------------------------------------------------------------------
###-------------------------------------------------------------------------------------------------------------------

