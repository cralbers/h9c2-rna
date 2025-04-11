#!/bin/bash
#SBATCH --account=pas2905
#SBATCH --job-name=bam_sort_dmso2
#SBATCH --mem=64gb
#SBATCH --time=40:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=R-%x.%j.out
#SBATCH --error=R-%x.%j.err


module load samtools

samtools view -b -S /dmso2_star/dmso2Aligned.out.sam > dmso2Aligned.out.bam
samtools sort -m 1000000000 -o /dmso2_star/dmso2.sorted /dmso2_star/dmso2Aligned.out.bam 