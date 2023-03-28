#!/bin/bash

# Set the path to the reference genome
REF=reference.fasta

# Index the reference genome
bwa index reference.fasta

# Create an output directory for the results
mkdir -p output

# Loop through each plate
for PLATE in Plate42B2 Plate135H10
do
    # Create a directory for the plate results
    mkdir -p output/$PLATE

    # Set the input filenames for this plate
    R1=data/$PLATE.R1.fastq.gz
    R2=data/$PLATE.R2.fastq.gz

    # Align the reads to the reference genome
    bwa mem $REF $R1 $R2 | samtools sort -O BAM -o output/$PLATE/aln.sorted.bam -

    # Index the alignment file
    samtools index output/$PLATE/aln.sorted.bam

    # Call variants
    bcftools mpileup -f $REF output/$PLATE/aln.sorted.bam | bcftools call --ploidy 1 -mv -Ob -o output/$PLATE/variants.bcf

    # Index the variants file
    bcftools index output/$PLATE/variants.bcf

    # Create consensus sequence
    bcftools consensus -f $REF output/$PLATE/variants.bcf > output/$PLATE/consensus.fasta

done
