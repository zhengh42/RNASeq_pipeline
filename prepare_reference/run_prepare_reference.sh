gencode.transcripts.fa=gencode.v25.transcripts.fa.gz
ref.genome=GRCh38.primary_assembly.genome.fa
gencode.gtf=gencode.v25.annotation.gtf
gtfToGenePred=../scripts/gtfToGenePred
match.pl=../scripts/match.pl
refFlatparser.pl=../scripts/refFlatparser.pl
fasta2unique_kmers.py=../scripts/fasta2unique_kmers.py

#######
# download data
#######
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.annotation.gtf.gz
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.transcripts.fa.gz
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh38.primary_assembly.genome.fa.gz
gzip *.gz

#######
# build index for Kallisto
#######
docker run -v $PWD:/mnt zhengh42/kallisto:0.43.0 kallisto index -i /mnt/gencode.v25.kallisto.idx /mnt/$gencode.transcripts.fa

#######
# build index for Salmon
#######
docker run -v $PWD:/mnt zhengh42/salmon:0.8.2 salmon index -i /mnt/gencode.v25.salmon.idx -t /mnt/$gencode.transcripts.fa

#######
# build index for STAR
#######
sjdbOverhang=47 # read length - 1
out_dir=STAR.index.sjdbOverhang${sjdbOverhang}
docker run -v $PWD:/mnt zhengh42/star:2.5.3a --runThreadN 16 --runMode genomeGenerate --genomeDir /mnt/$out_dir --genomeFastaFiles /mnt/$ref.genome --sjdbGTFfile /mnt/$gencode.gtf --sjdbOverhang $sjdbOverhang

#######
# build index for RSEM
#######
docker run -v $PWD:/mnt zhengh42/rsem:1.3.0 rsem-prepare-reference -p 16 --gtf /mnt/$gencode.gtf  /mnt/$ref.genome /mnt/RSEM/gencode.v25

#######
# Generate RefFlat file
#######
$gtfToGenePred $gencode.gtf gencode.gtf.tmp
less $gencode.gtf | awk '$3=="transcript"' | cut -f9 | awk '{print $2,$4}' | sed 's/"//g;s/;//g;s/\s\+/\t/' > gene_trans.id
less gencode.gtf.tmp | $match.pl gene_trans.id 2 - 1 | paste - gencode.gtf.tmp | cut -f1,2,4- >  gencode.gtf.refFlat

#######
# Assign gene type
#######
less $gencode.gtf | grep -v ^# | awk '$3~/transcript/' | cut -f9 | awk 'OFS="\t"{print $2,$4,$6,$12}' | sed 's/[";]//g' | awk 'OFS="\t"{$5="other";$6="other";if($4~/non_coding|3prime_overlapping_ncRNA|antisense|bidirectional_promoter_lncRNA|lincRNA|macro_lncRNA|sense_intronic|sense_overlapping|TEC/)$6="lncRNA";if($4~/protein_coding/)$6="proteincoding";if($3~/protein_coding/)$5="proteincoding";if($3~/non_coding|3prime_overlapping_ncRNA|antisense|bidirectional_promoter_lncRNA|lincRNA|macro_lncRNA|sense_intronic|sense_overlapping|TEC/)$5="lncRNA";  print $0}'| sed '1i gene_id\ttranscript\tgenetype1\ttranscripttype1\tgenetype\ttranscripttype' > gene_transcript_type

less gene_transcript_type | sed 1d | cut -f1,3,5 | sort | uniq | sed '1i gene_id\tdetail\ttype'> gene_type

#######
# Get mapping file for transcripts and genes
#######
less $gencode.transcripts.fa | egrep '^>' | sed 's/>//' | awk -F "|" 'OFS="|"{print $0,$2}' | sed 's/||/|\t/' > tx2gene.txt

#######
# Get gene features
#######
# Number of transcripts, number of exons, transcript length
$refFlatparser.pl gencode.gtf.refFlat
$match.pl gene_type 1 gencode.gtf.refFlat.geneinfo  1 | paste - gencode.gtf.refFlat.geneinfo | cut -f2- > gencode.gtf.refFlat.geneinfo.txt

# unique kmers
# source /labs/gevaertlab/tools/cgat/conda-install/bin/activate cgat-s
python $fasta2unique_kmers.py --input-fasta="$gencode.transcripts.fa" --method="gene" --genemap="tx2gene.txt" --kmer-size=31 > gencode.transcripts.fa.uniquekmer

less gencode.transcripts.fa.uniquekmer | egrep 'fraction_unique|ENSG' | sed '1s/^/gene_/' | $match.pl - 1 gencode.gtf.refFlat.geneinfo.txt 3 | paste gencode.gtf.refFlat.geneinfo.txt - | cut -f1-6,8- > gencode.gtf.refFlat.geneinfo_b.txt
