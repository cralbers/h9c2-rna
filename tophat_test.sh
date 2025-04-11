#!/bin/bash
#SBATCH --account=pas2905
#SBATCH --job-name=tophat
#SBATCH --time=40:00:00
#SBATCH --mem=8000mb
#SBATCH --mail-type=ALL

cd /users/PAS2905/coraalbers/H9c2_RNA_seq

module load tophat/2.1.1
module load bowtie2
module load samtools
module load cufflinks/2.2.1

tophat -o X12656_DMSO_S1_thout /users/PAS2905/coraalbers/H9c2_RNA_seq/rat_genome/Rattus_norvegicus/Ensembl/RGSC3.4/Sequence/Bowtie2Index/genome fastq/X12656_SmithS_H9c2-Cells_DMSO_48-Hrs_V1Q_1_S1_L002_R1_001.fastq.gz fastq/X12656_SmithS_H9c2-Cells_DMSO_48-Hrs_V1Q_1_S1_L002_R2_001.fastq.gz


