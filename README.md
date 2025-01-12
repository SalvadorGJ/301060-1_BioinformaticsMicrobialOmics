*2024W 301060-1 Bioinformatics for microbial omics - Computational genomics, pangenomics, metagenomics and multi-omics analysis*

# Assembly of E.coli genome from surveillance project

**Aim**: Assemble the genome sequence of a potentially pathogenic E.coli strain from this surveillance project: [Greig et al., Gigascience 2019](https://pubmed.ncbi.nlm.nih.gov/31433830/)

## Input data

FASTQ files were retrieved from Sequence Read Archive (SRA) by fasterq-dump (from the sratoolkit v3.1.1 suite) using the following accessions numbers:

SRR15361790: paired-end data from Illumina HiSeq 2500.
SRR17645346: long-read data from the MinION and base-called using Albacore. 

## Processing of Illumina sequence data

Quality reports of the raw FASTQ files were made with FastQC v0.12.1. The FASTQ reads were then processed using Trimmomatic v0.39 to (1) remove adapters, (2) remove bases with a PHRED score of <30 from the leading and trailing ends, (3) trim when the average quality per base drops below 20 within a 30-base wide sliding window, (4) and discard reads <90 bp after quality trimming. FastQC reports were generated again for the processed FASTQ files.

## Genome assembly and evaluation

We used SPAdes v4.0.0 to assemble a first draft using only the processed Illumina reads as input (spades_filtered_illumina-only) with k-mer lengths of 21, 33, 55, 65, 77, 83, and 91. Then we used again SPAdes v4.0.0 to generate a hybrid assembly using both the processed Illumina reads and the raw ONT long-reads as input (spades_filtered_hybrid) with k-mer lengths of 21, 33, 55, 65, 77, 83, and 91. Finally, we used Flye v2.9.5 to assemble a third draft using only the ONT long-reads as input (flye_nanopore-only). To compare the generated assemblies, statistics of each draft were computed using QUAST v5.3.0.
