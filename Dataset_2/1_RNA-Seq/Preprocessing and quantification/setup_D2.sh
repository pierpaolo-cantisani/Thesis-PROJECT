#!/bin/bash

### Creating directories ###
mkdir -p /mnt/data/home/ubuntu/dataset2/{fastq,not_trimmed,trimmed,fastq_trimmed,references,salmon_output}


#Environment for salmon:
conda env create -f rna_seq.yml
#Conda activation
source ~/.bashrc
conda activate rna_seq



### Downloading and indexing reference transcriptome for salmon quantification ###

#Downloading
wget -P references/ \
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.transcripts.fa.gz



# Indexing
salmon index \
  -t references/gencode.v49.transcripts.fa.gz \
  -i references/salmon_index_hg38_49 \
  --threads 2
conda deactivate

#Otherwise do :
#conda activate rna_seq
#salmon index \
#  -t ~/Multi-omics-analysis/01_rnaseq/references/gencode.v49.transcripts.fa.gz \
#  -i ~/Multi-omics-analysis/01_rnaseq/references/salmon_index_hg38_49 \
#  --threads 2
