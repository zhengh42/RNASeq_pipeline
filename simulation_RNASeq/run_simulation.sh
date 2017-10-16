work_dir=$PWD
TCGA_RNASeq_dir=../TCGA_RNASeq # TCGA data directory
rsem_reference_name=../prepare_reference/RSEM/gencode.v25 # pre-built reference index for RSEM
file=$1 # sample name
theta0=$2 # noise level (0.05, 0.15, 0.3)
N=60000000 # total number of reads to be simulated
out_dir=$work_dir/seq  # directory for storing simulated output files
estimated_model_file=$file.stat/$file.model # the model of real data
estimated_isoform_results=$file.isoforms.results # the quantification result of real data

#######
# run simualtions
#######
docker run -v $PWD/seq:/mnt/seq -v $TCGA_RNASeq_dir/RSEM:/mnt/model -v ../prepare_reference:/mnt/ref zhengh42/rsem:1.3.0 rsem-simulate-reads /mnt/ref/RSEM/gencode.v25 /mnt/model/$estimated_model_file /mnt/model/$estimated_isoform_results $theta0 $N /mnt/seq/all_${theta0}_${file}

#######
# compress output fastq files
#######
gzip $out_dir/all_${theta0}_${file}_1.fq
gzip $out_dir/all_${theta0}_${file}_2.fq
