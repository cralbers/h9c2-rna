#!/bin/bash
#SBATCH --account=pas2905
#SBATCH --job-name=fastqc
#SBATCH --mem=64gb
#SBATCH --time=40:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=R-%x.%j.out
#SBATCH --error=R-%x.%j.err

module load fastqc

cd X12661
fastqc -o X12661_qc X12661_SmithS_H9c2-Cells_DMSO_48-Hrs_V1Q_1_S5_L002_R1_001.fastq.gz X12661_SmithS_H9c2-Cells_DMSO_48-Hrs_V1Q_1_S5_L002_R2_001.fastq.gz
cd ../

cd X12662
fastqc -o X12662_qc X12662_SmithS_H9c2-Cells_DMSO_48-Hrs_V1Q_1_S5_L002_R1_001.fastq.gz X12662_SmithS_H9c2-Cells_DMSO_48-Hrs_V1Q_1_S5_L002_R2_001.fastq.gz
cd ../

cd X12663
fastqc -o X12663_qc X12663_SmithS_H9c2-Cells_DMSO_48-Hrs_V1Q_1_S5_L002_R1_001.fastq.gz X12663_SmithS_H9c2-Cells_DMSO_48-Hrs_V1Q_1_S5_L002_R2_001.fastq.gz
cd ../

cd X12664
fastqc -o X12664_qc X12664_SmithS_H9c2-Cells_DMSO_48-Hrs_V1Q_1_S5_L002_R1_001.fastq.gz X12664_SmithS_H9c2-Cells_DMSO_48-Hrs_V1Q_1_S5_L002_R2_001.fastq.gz
cd ../

cd X12665
fastqc -o X12665_qc X12665_SmithS_H9c2-Cells_DMSO_48-Hrs_V1Q_1_S5_L002_R1_001.fastq.gz X12665_SmithS_H9c2-Cells_DMSO_48-Hrs_V1Q_1_S5_L002_R2_001.fastq.gz
cd ../

cd X12666
fastqc -o X12666_qc X12666_SmithS_H9c2-Cells_DMSO_48-Hrs_V1Q_1_S5_L002_R1_001.fastq.gz X12666_SmithS_H9c2-Cells_DMSO_48-Hrs_V1Q_1_S5_L002_R2_001.fastq.gz
cd ../

cd X12667
fastqc -o X12667_qc X12667_SmithS_H9c2-Cells_DMSO_48-Hrs_V1Q_1_S5_L002_R1_001.fastq.gz X12667_SmithS_H9c2-Cells_DMSO_48-Hrs_V1Q_1_S5_L002_R2_001.fastq.gz
cd ../

cd X12668
fastqc -o X12668_qc X12668_SmithS_H9c2-Cells_DMSO_48-Hrs_V1Q_1_S5_L002_R1_001.fastq.gz X12668_SmithS_H9c2-Cells_DMSO_48-Hrs_V1Q_1_S5_L002_R2_001.fastq.gz
cd ../

cd X12669
fastqc -o X12669_qc X12669_SmithS_H9c2-Cells_DMSO_48-Hrs_V1Q_1_S5_L002_R1_001.fastq.gz X12669_SmithS_H9c2-Cells_DMSO_48-Hrs_V1Q_1_S5_L002_R2_001.fastq.gz
cd ../

cd X12670
fastqc -o X12670_qc X12670_SmithS_H9c2-Cells_DMSO_48-Hrs_V1Q_1_S5_L002_R1_001.fastq.gz X12670_SmithS_H9c2-Cells_DMSO_48-Hrs_V1Q_1_S5_L002_R2_001.fastq.gz
cd ../

cd X12671
fastqc -o X12671_qc X12671_SmithS_H9c2-Cells_DMSO_48-Hrs_V1Q_1_S5_L002_R1_001.fastq.gz X12671_SmithS_H9c2-Cells_DMSO_48-Hrs_V1Q_1_S5_L002_R2_001.fastq.gz
cd ../

cd X12672
fastqc -o X12672_qc X12672_SmithS_H9c2-Cells_DMSO_48-Hrs_V1Q_1_S5_L002_R1_001.fastq.gz X12672_SmithS_H9c2-Cells_DMSO_48-Hrs_V1Q_1_S5_L002_R2_001.fastq.gz
cd ../

cd X12673
fastqc -o X12673_qc X12673_SmithS_H9c2-Cells_DMSO_48-Hrs_V1Q_1_S5_L002_R1_001.fastq.gz X12673_SmithS_H9c2-Cells_DMSO_48-Hrs_V1Q_1_S5_L002_R2_001.fastq.gz
cd ../
