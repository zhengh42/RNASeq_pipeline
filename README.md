# RNASeq_pipeline

## Process each sample
run_workflow.sh describes the workflow for analysis of each sample from TCGA pan-cancer RNASeq dataset.

0. Before running the workflow, download and index reference files first. See prepare_reference/run_prepare_reference.sh
1. Download bam files from GDC data portal. Sample selection is described in TCGA_RNASeq/run_samplingTCGA.R
2. Convert to fastq
3. QC by trim-galore
4. Kallisto
4.1 with sequence bias estimation
4.2 without sequence bias estimation
5. Salmon
5.1 with sequence and gc bias estimation
5.2 without sequence and gc bias estimation
6. STAR
7. RSEM
8. HTSeq
9. featureCounts

Simulation for RNASeq reads is described in simulation_RNASeq/run_simulation.sh
Processing for simulated dataset is the same as TCGA pan-cancer dataset (starting at step 3)

## Further processing, batch-level
Kallisto and Salmon measure the expression level of each transcript by default. To get gene-level expression results, the package tximport was used. The script is found in scripts/tximport.R

HTSeq and featureCounts output read count for each gene. TPM values were generated from read counts with scripts/run_HTSeq_count_to_tpm.sh and scripts/run_featurecounts_count_to_tpm.sh.


