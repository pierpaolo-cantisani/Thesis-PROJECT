#!/bin/bash

#Externally activate environment:
#conda activate rna_seq

### 1. Downloading and merging fastq files ###
 SRR codes are in "SRR_Acc_List.txt", obtained from the SRA website.

# Downloading
while read id; do
  echo "Downloading and compressing sample: $id"

  #Downloading
  prefetch "$id" -O ~/Multi-omics-analysis/01_rnaseq/fastq
  fasterq-dump \
   "$HOME/Multi-omics-analysis/01_rnaseq/fastq/${id}" --outdir ~/Multi-omics-analysis/01_rnaseq/fastq
  #Compressing
  gzip ~/Multi-omics-analysis/01_rnaseq/fastq/${id}.fastq

  # Removing the prefetch directory
  rm -r ~/Multi-omics-analysis/01_rnaseq/fastq/${id}
done < ~/Multi-omics-analysis/01_rnaseq/SRR_Acc_List.txt


# Merging
#The experiment is composed of 48 runs, 4 for each of 12 samples. Merging each quadruple of files:
mapfile -t srr_codes < ~/Multi-omics-analysis/01_rnaseq/SRR_Acc_List.txt

for ((i=0; i<${#srr_codes[@]}; i+=4)); do
   j=$((i/4 + 1))
   echo "Merging sample $j"

   zcat \
      "$HOME/Multi-omics-analysis/01_rnaseq/fastq/${srr_codes[i]}.fastq.gz" \
      "$HOME/Multi-omics-analysis/01_rnaseq/fastq/${srr_codes[i+1]}.fastq.gz" \
      "$HOME/Multi-omics-analysis/01_rnaseq/fastq/${srr_codes[i+2]}.fastq.gz" \
      "$HOME/Multi-omics-analysis/01_rnaseq/fastq/${srr_codes[i+3]}.fastq.gz" \
      | gzip > "$HOME/Multi-omics-analysis/01_rnaseq/fastq/sample_${j}_merged.fastq.gz"
  
   #Removing non-merged files
   rm "$HOME/Multi-omics-analysis/01_rnaseq/fastq/${srr_codes[i]}.fastq.gz" "$HOME/Multi-omics-analysis/01_rnaseq/fastq/${srr_codes[i+1]}.fastq.gz" "$HOME/Multi-omics-analysis/01_rnaseq/fastq/${srr_codes[i+2]}.fastq.gz" "$HOME/Multi-omics-analysis/01_rnaseq/fastq/${srr_codes[i+3]}.fastq.gz"
done


### 2. Quality Control

#Performing QC:
fastqc -t 2 -o ~/Multi-omics-analysis/01_rnaseq/fastqc_results/not_trimmed ~/Multi-omics-analysis/01_rnaseq/fastq/*.fastq.gz
#Merging QC results:
multiqc ~/Multi-omics-analysis/01_rnaseq/fastqc_results/not_trimmed -n full_qc_report_not_trimmed -o ~/Multi-omics-analysis/01_rnaseq/fastqc_results/


### 3. Trimming and second Quality Control
for ((i=1; i<13; i+=1)); do
  echo "Trimming sample: $i"

  fastp \
    -w 2 \
    --in1 "$HOME/Multi-omics-analysis/01_rnaseq/fastq/sample_${i}_merged.fastq.gz" \
    --out1 "$HOME/Multi-omics-analysis/01_rnaseq/fastq_trimmed/sample_${i}_trimmed.fastq.gz" \
    --html "$HOME/Multi-omics-analysis/01_rnaseq/fastq_trimmed/sample_${i}_fastp.html" \
    --json "$HOME/Multi-omics-analysis/01_rnaseq/fastq_trimmed/sample_${i}_fastp.json" \
    --qualified_quality_phred 20 \
    --length_required 40

  rm "$HOME/Multi-omics-analysis/01_rnaseq/fastq/sample_${i}_merged.fastq.gz"

done

#Performing QC:
fastqc -t 2 -o ~/Multi-omics-analysis/01_rnaseq/fastqc_results/trimmed ~/Multi-omics-analysis/01_rnaseq/fastq_trimmed/*.fastq.gz
#Merging QC results:
multiqc ~/Multi-omics-analysis/01_rnaseq/fastqc_results/trimmed -n full_qc_report_trimmed -o ~/Multi-omics-analysis/01_rnaseq/fastqc_results/



### 4. Mapping and quantification ###
##Using Salmon

for ((i=1; i<13; i+=1)); do
  echo "Quantifying sample: ${i}"

  salmon quant -l A \
        -i ~/Multi-omics-analysis/01_rnaseq/references/salmon_index_hg38_49 \
        -r "$HOME/Multi-omics-analysis/01_rnaseq/fastq_trimmed/sample_${i}_trimmed.fastq.gz" \
        -p 2 \
        --validateMappings \
        -o "$HOME/Multi-omics-analysis/01_rnaseq/salmon_output/sample_${i}_quant"

done

#At this point it's important to check the results: statistics and file format
#less ~/Multi-omics-analysis/01_rnaseq/salmon_output/sample_1_quant/lib_format_counts.json
#head -n 20 ~/Multi-omics-analysis/01_rnaseq/salmon_output/sample_1_quant/quant.sf | column -t | less

#Renaming directories for differential analysis:
n=1
for ((i=1; i<13; i+=2)); do

  mv "$HOME/Multi-omics-analysis/01_rnaseq/salmon_output/sample_${i}_quant" "$HOME/Multi-omics-analysis/01_rnaseq/salmon_output/rep_${n}_TB_quant"
  mv "$HOME/Multi-omics-analysis/01_rnaseq/salmon_output/sample_$((i+1))_quant" "$HOME/Multi-omics-analysis/01_rnaseq/salmon_output/rep_${n}_NI_quant"
  ((n++))

done