#!/bin/bash

#Externally activate environment:
#conda activate rna_seq

### 1. Downloading and merging fastq files ###
# SRR codes are in "SRR_Acc_List.txt", obtained from the SRA website.

# Downloading
while read id; do
  echo "Downloading and compressing sample: $id"

  #Downloading
  prefetch "$id" -O fastq
  fasterq-dump \
   "fastq/${id}" --outdir fastq
  #Compressing
  gzip fastq/${id}.fastq

  # Removing the prefetch directory
  rm -r fastq/${id}
done < SRR_Acc_List.txt



### 2. Quality Control

#Performing QC:
fastqc -t 2 -o fastqc_results/not_trimmed fastq/*.fastq.gz
#Merging QC results:
multiqc fastqc_results/not_trimmed -n full_qc_report_not_trimmed -o fastqc_results/


### 3. Trimming and second Quality Control
for ((i=1988; i<=2026; i++)); do
  #Skipping missing samples
  case $i in
    1989|2007|2009|2010|2011)
      echo "Skipping SRR634${i}"
      continue
      ;;
  esac
  
  echo "Trimming SRR: $i"

  fastp \
    -w 2 \
    --in1 "fastq/SRR634${i}.fastq.gz" \
    --out1 "fastq_trimmed/SRR634${i}_trimmed.fastq.gz" \
    --html "fastq_trimmed/SRR634${i}_fastp.html" \
    --json "fastq_trimmed/SRR634${i}_fastp.json" \
    --qualified_quality_phred 20 \
    --length_required 40 && \
  rm "fastq/SRR634${i}.fastq.gz"

done

#Performing QC:
fastqc -t 2 -o fastqc_results/trimmed fastq_trimmed/*.fastq.gz
#Merging QC results:
multiqc fastqc_results/trimmed -n full_qc_report_trimmed -o fastqc_results/



### 4. Mapping and quantification ###
##Using Salmon

for ((i=1988; i<=2026; i++)); do
  #Skipping missing samples
  case $i in
    1989|2007|2009|2010|2011)
      echo "Skipping SRR634${i}"
      continue
      ;;
  esac

  echo "Quantifying sample: ${i}"

  salmon quant -l A \
        -i references/salmon_index_hg38_49 \
        -r "fastq_trimmed/SRR634${i}_trimmed.fastq.gz" \
        -p 2 \
        --validateMappings \
        -o "salmon_output/sample_${i}_quant_D2"

done
