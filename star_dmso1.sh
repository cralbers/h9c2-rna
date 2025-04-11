#!/bin/bash
#SBATCH --account=pas2905
#SBATCH --job-name=star_dmso1
#SBATCH --mem=64gb
#SBATCH --time=40:00:00
#SBATCH --mail-type=ALL
#SBATCH --output=R-%x.%j.out
#SBATCH --error=R-%x.%j.err




cd /users/PAS2905/coraalbers/H9c2_RNA_seq

module load star/2.5.2a

STAR --genomeDir /users/PAS2905/coraalbers/h9c2_RNA_seq/rat_index \
--runThreadN 6 \
--readFilesIn fastq/X12656_SmithS_H9c2-Cells_DMSO_48-Hrs_V1Q_1_S1_L002_R1_001.fastq.gz fastq/X12656_SmithS_H9c2-Cells_DMSO_48-Hrs_V1Q_1_S1_L002_R2_001.fastq.gz \
--readFilesCommand zcat \
