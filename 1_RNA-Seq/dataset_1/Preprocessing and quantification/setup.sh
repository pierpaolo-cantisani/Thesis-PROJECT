#!/bin/bash


### Creating directories ###
mkdir -p ~/Multi-omics-analysis/01_rnaseq/{fastq,fastqc_results/not_trimmed,fastq_trimmed,fastqc_results/trimmed,salmon_output,references,DE_results}
mkdir -p ~/Multi-omics-analysis/02_wgbs/{DM_results,references}
mkdir -p ~/Multi-omics-analysis/03_integration_analysis/integration_results


### Downloading and indexing reference transcriptome for salmon quantification ###

#Downloading
wget -P ~/Multi-omics-analysis/01_rnaseq/references/ \
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.transcripts.fa.gz

#Environment for salmon:
conda env create -f ~/Multi-omics-analysis/envs/rna_seq.yml
#Activating "rna_seq" conda environment before:
source activate rna_seq
# Indexing
salmon index \
  -t ~/Multi-omics-analysis/01_rnaseq/references/gencode.v49.transcripts.fa.gz \
  -i ~/Multi-omics-analysis/01_rnaseq/references/salmon_index_hg38_49 \
  --threads 2
conda deactivate

#Otherwise do :
#conda activate rna_seq
#salmon index \
#  -t ~/Multi-omics-analysis/01_rnaseq/references/gencode.v49.transcripts.fa.gz \
#  -i ~/Multi-omics-analysis/01_rnaseq/references/salmon_index_hg38_49 \
#  --threads 2


### Downloading the reference for the DESeq2 analysis:
wget -P ~/Multi-omics-analysis/01_rnaseq/references/ \
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz


### Donloading reference for the DMR analysis
wget -P ~/Multi-omics-analysis/02_wgbs/references/ \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
gunzip ~/Multi-omics-analysis/02_wgbs/references/hg19ToHg38.over.chain.gz