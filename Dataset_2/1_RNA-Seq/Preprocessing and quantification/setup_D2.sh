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



## Index fot Dataset2: all unspliced

# 1. GTF GENCODE v49
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz
gunzip gencode.v49.annotation.gtf.gz

# 2. Primary assembly genome (standard chromosomes + main scaffold, no haplotypes)
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz


# Now the R script (in the correct conda environment)
Rscript EisaR.R 2>&1 | tee EisaR.log

# Indexing
salmon index -t gentrome_with_decoy.fa \
             -i salmon_idx_spliced_intron_hg38 \
             -k 31 \
             --threads 2        


# Can add decoy later
#cat gentrome_spliced_intron.fa GRCh38.primary_assembly.genome.fa \
#  > gentrome_with_decoy.fa

# decoys.txt = names of sequences
#grep "^>" GRCh38.primary_assembly.genome.fa | cut -d " " -f1 | sed 's/>//' \
#  > decoys.txt