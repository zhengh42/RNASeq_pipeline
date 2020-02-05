# RNASeq pipeline

This repository describes common pipelines for RNA sequencing analysis.

Outline:

1. Transcriptome and references

  -  Human (homo sapiens)
    - GENCODE: the commonly-used transcriptome (https://www.gencodegenes.org/)
    - GENCODE + NONCODE: for investigation long non-coding RNAs not included in GENCODE. NONCODE (http://www.noncode.org/) is an integrated knowledge database dedicated to non-coding RNAs (excluding tRNAs and rRNAs).
  -  Mice (mus musmusculus)
    - GENCODE
    
2. Pipelines    

  - Quality control
    - Trim Galore
      - Unaligned reads
      - Low-quality bases
      - Adaptor sequences
      - Very short reads

    - RSeQC
      - Mapping rates
      - Unique reads
      - Reads distribution

  - Human gene expression quantification
    - Kalliso
    - STAR and RSME
    - STAR and featureCounts

  - Other gene expression quantification
    - Exogenous viruses
    - Endogenous retrovirus 
  
  - Differential expression analysis
    - DESeq2
    
  - Module network analysis
    - AMARETTO
    - WGCNA
  
## Prepare reference

Before running the pipelines, download genome and transcriptome refernces and prepare indexes for each tool first. See prepare_reference.sh for detailed instructions.

GENCODE references can be downloaded directly from their website. 

If including NONCODE in the analysis and combining with GENCODE transcriptome, the redundant records between the two need to be removed. The pprepare_reference script takes care of that, genetating gencodev32noncodev5.fa and gencodev32noncodev5.annotation.gtf for downstream analysis. For convenience the two files can be found in Stanford Medicine Box: https://stanfordmedicine.app.box.com/folder/102547052245

## Quality control

## Gene expression quantification




