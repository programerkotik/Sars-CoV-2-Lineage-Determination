# Sars-CoV-2 Lineage Determination
This document describes the workflow used for Sars-CoV-2 lineage determination in the course Analytical Methods in Cancer Genomics 2023.

## Workflow
The workflow consisted of three main steps:

1. Alignment of reads to reference genome
2. Variant calling
3. Lineage determination

### 1. Alignment of reads to reference genome
The first step was to align the reads to the reference genome using the BWA tool. The reference genome was first indexed using the command ```bwa index reference.fasta```. Then, the reads were aligned to the reference genome using the command ```bwa mem $REF $R1 $R2 | samtools sort -O BAM -o output/$PLATE/aln.sorted.bam -```, where $REF is the reference Sars-CoV-2 genome file, $R1 and $R2 are the input read files, and $PLATE is the plate name. The BAM file was then indexed using the command ```samtools index output/$PLATE/aln.sorted.bam```.

### 2. Variant calling
The next step was to call variants using the bcftools tool. The command used was ```bcftools mpileup -f $REF output/$PLATE/aln.sorted.bam | bcftools call --ploidy 1 -mv -Ob -o output/$PLATE/variants.bcf```. The resulting BCF file was then indexed using the command ```bcftools index output/$PLATE/variants.bcf```. Finally, the consensus sequence was generated using the command ```bcftools consensus -f $REF output/$PLATE/variants.bcf > output/$PLATE/consensus.fasta```.

### 3. Lineage determination
To determine the lineage of the sample, the Nextclade tool was used. The consensus sequence was uploaded to the Nextclade web application, and the lineage was determined as B.1.1.7 (Alpha) for the first dataset and BA.1.18 (Omicron) for the second dataset. The results were saved in the csv file nextclade.csv.

## Usage
To run the full workflow the shell analysis.sh script can be used. Requirements are: 
- [BWA](https://github.com/lh3/bwa)
- [SAMtools](http://www.htslib.org/)
- [BCFTOOLS](https://samtools.github.io/bcftools/)

#### My notes
1. Are there a way to find out from the data information about single-end or paired-end reads (in my case it is paired-end because we have 2 read groups?), short or long reads, and the sequencing technology used?
My idea is that I can look at the header but it provides us with only *seemingly* irrelavant information about machine name, flowcell id, lane, and tile + coords. To determine which aligner to use, I would need to know if the reads are single-end or paired-end, and if they are short or long reads. 
2. ```zcat``` doesn't work on mac, so I used ```gzcat``` instead. installed through ```brew install coreutils```
3. Found out the length of reads using ```gzcat Plate135H10.R1.fastq.gz | awk '{if(NR%4==2) print length($0)}' | sort | uniq -c``` command and found out that the reads are up to 10 kbp long.