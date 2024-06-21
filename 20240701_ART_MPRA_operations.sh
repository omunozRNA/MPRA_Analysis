#!/usr/bin/env bash
#BSUB -J OM_ART_sort
#BSUB -e logs/UMI_%I.%J.err
#BSUB -o logs/UMI_%I.%J.out
#BSUB -R "rusage[mem=160] span[hosts=1]"
#BSUB -n 16

# Code Description: The shell script will complete several steps for the purposes of Differential Expression Analysis
# The following steps are as follows:
# Step 1: load/unload relevant modules in bodhi HPC. (May or may not be reproducible on different servers, modify protocol accordingly)
# Step 2: Build an Index (essentially a 'custom_genome') with bowtie
# Step 3: Remove sequencing adaptor sequences from data, trimming the appropriate sequences will aid in the removal of any unexpected artifacts and should help with alignment
# Step 4: Alignment -> bowtie2 will take care of calling reads through 
# Step 5: Sorting my sam file and turning into bam files 
# Step 6: Deduplicating reads

module load bowtie
module load bowtie2
module load samtools
module load bedtools
module load anaconda
module unload perl
source activate ReadProcessing

cd /beevol/home/munozosc/20240701_Colony_LibPrep_Data/bowtie/

# Step 2: Building the custom index
bowtie2-build -f --threads 16 /beevol/home/munozosc/20240508_Colony_LibPrep_Data/bowtie/ART_seqs.fa Bowtie_MPRA_Ref

mkdir -p logs

# Step 3: Cutadapting...
# Storing callable prefixes into ARTs
ARTs=(
  "ART8_SMG1i_1" "ART8_SMG1i_2" "ART8_SMG1i_3" 
  "ART8_DMSO_1" "ART8_DMSO_2" "ART8_DMSO_3" 
  "ARTi8_SMG1i_1" "ARTi8_SMG1i_2" "ARTi8_SMG1i_3" 
  "ARTi8_DMSO_1" "ARTi8_DMSO_2" "ARTi8_DMSO_3" 
  "ARTx8_SMG1i_1" "ARTx8_SMG1i_2" "ARTx8_SMG1i_3" 
  "ARTx8_DMSO_1" "ARTx8_DMSO_2" "ARTx8_DMSO_3"
)

# Location of the fastqs to be processed
fastqs="MPRA_libs/"
trim_log="logs/cutadapt.log"
ind="/beevol/home/munozosc/20240701_Colony_LibPrep_Data/bowtie/Bowtie_MPRA_Ref"
fwdplus="CCGTCGTGGTCCTTGTACGTCGAATTC"
fwdminus="GAATTCGACGTACAAGGACCACGACGG"
revplus="AAACTGGGGCACAGCCTCGAG"
revminus="CTCGAGGCTGTGCCCCAGTTT"
m=10
O=6

# Step 3.1: Trimming with cutadapt. Removing the 5' ends first
for i in "${ARTs[@]}"; do
    r1="${fastqs}${i}1.fq.gz"
    r2="${fastqs}${i}2.fq.gz"
    r1_trim5="${fastqs}${i}1.5trimmed.fq"
    r2_trim5="${fastqs}${i}2.5trimmed.fq"
    
    cutadapt \
    --discard-untrimmed \
    -m $m \
    -O $O \
    -g $fwdplus \
    -G $revplus \
    -o $r1_trim5 \
    -p $r2_trim5 \
    $r1 \
    $r2 \
    > $trim_log
done

# Step 3.2: Trimming with cutadapt, removing the 3' ends second
for i in "${ARTs[@]}"; do
    r1_trim5="${fastqs}${i}1.5trimmed.fq"
    r2_trim5="${fastqs}${i}2.5trimmed.fq"
    r1_trim="${fastqs}${i}1.53trimmed.fq"
    r2_trim="${fastqs}${i}2.53trimmed.fq"
    
    cutadapt \
    --discard-untrimmed \
    -m $m \
    -O $O \
    -g $revminus \
    -G $fwdminus \
    -o $r1_trim \
    -p $r2_trim \
    $r1_trim5 \
    $r2_trim5 \
    > $trim_log
done

# Step 4: Aligning to custom genome
for i in "${ARTs[@]}"; do
    r1_trim="${fastqs}${i}1.53trimmed.fq"
    r2_trim="${fastqs}${i}2.53trimmed.fq"
    sam="${i}.sam"
    
    bowtie2 -q --end-to-end --no-discordant --no-unal -p 16 -D 15 -x ${ind} -1 ${r1_trim} -2 ${r2_trim} -S ${sam}
done

# Step 5: Sorting SAM files and converting to BAM
for i in "${ARTs[@]}"; do
    input_file="${i}.sam"
    output_file="${i}.bam"
    
    samtools sort -@ 4 "$input_file" -o "$output_file"
done

# Step 6: Deduplicating reads
module load picard

for i in "${ARTs[@]}"; do
    input_file="${i}.bam"
    output_file="${i}.deduplicated.bam"
    metrics_file="${i}.metrics.txt"
    
    java -jar $PICARD MarkDuplicates \
    I="$input_file" \
    O="$output_file" \
    M="$metrics_file" \
    REMOVE_DUPLICATES=true
    
    echo "UMI-based deduplication complete. Output saved to $output_file."
done

