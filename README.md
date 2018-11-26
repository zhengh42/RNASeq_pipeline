# RNASeq_pipeline

run_workflow.sh describes each step of the analysis. 

### Prepare reference

Before running the workflow, download and index reference files first. See prepare_reference/run_reference.sh

### Download RNA-Seq datasets

Both un-stranded and reverse-stranded RNA-Seq data from TCGA samples were downloaded from ISB Cancer Genomics Cloud (ISB-CGC). The other reverse-stranded dataset was downloaded from NCBI Sequence Read Archive (SRA) under the accession PRJEB11797. 

### Process real RNA-Seq datasets

Reads QC were performed with Trim galore (27), with the setting “-q 20 --stringency 3 --gzip --length 20 --paired”. Afterwards the reads were mapped to the human transcriptome (both GENCODE and GENCODE combined with NONCODE) by STAR, and were further processed by RSEM (28) (version 1.3.0) to obtain gene and transcript expression. Stand-specific option was set as `-- forward-prob 0.5` for un-stranded samples and `--forward-prob 0` for reverse-stranded samples. Refer to "QC", "STAR", and "RSEMAfterSTAR" sections in run_workflow.sh.

### Simulation of RNA-Seq reads

Simulation of RNA-Seq reads was performed with the RSEM command rsem-simulate-reads, using the model information and quantification results of real samples. The total number of simulated reads for each sample is 60 million. The simulated samples with pre-defined gene expression levels serve as the “ground truth” for the evaluation of other pipelines. Refer to "simulation" section in run_workflow.sh.

### Process simulated RNA-Seq datasets

1. Pseudoalignment methods

-  Kallisto
-  Salmon

2. Alignment-based methods

    - Alignment

        - STAR
        - Subread
        - HISAT2

    - Quantification

        - HTSeq
            - HTSeq after STAR
            - HTSeq after Subread
            - HTSeq after HISAT2

        - featureCounts
            - featureCounts after STAR
            - featureCounts after Subread
            - featureCounts after HISAT2

