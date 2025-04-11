#!/bin/bash
#SBATCH --account=pas2905
#SBATCH --job-name=hisat2_align_r54
#SBATCH --mem=64gb
#SBATCH --time=40:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=R-%x.%j.out
#SBATCH --error=R-%x.%j.err


module load hisat2/2.1.0



hisat2 -x hisat2/rat_tran_r54 \
         -1 fastq/X12656/X12656_*_R1_001.fastq.gz \
         -2 fastq/X12656/X12656_*_R2_001.fastq.gz \
         -p 8 -S hisat2_align/X12656_hisat2_align.sam