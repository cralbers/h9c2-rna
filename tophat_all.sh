#!/bin/bash
#SBATCH --account=pas2905
#SBATCH --job-name=tophat_mapping
#SBATCH --mem=64gb
#SBATCH --time=40:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=R-%x.%j.out
#SBATCH --error=R-%x.%j.err


module load tophat/2.1.1
module load bowtie2/2.4.1
module load samtools/1.10
module load cufflinks/2.2.1

tophat -o X12656_DMSO_S1_thout \
/users/PAS2905/coraalbers/H9c2_RNA_seq/rat_genome/Rattus_norvegicus/Ensembl/RGSC3.4/Sequence/Bowtie2Index/genome \
fastq/X12656_SmithS_H9c2-Cells_DMSO_48-Hrs_V1Q_1_S1_L002_R1_001.fastq.gz fastq/X12656_SmithS_H9c2-Cells_DMSO_48-Hrs_V1Q_1_S1_L002_R2_001.fastq.gz


for fn in fastq/X126{56..73};
do
samp=`basename ${fn}`
echo "Processing sample ${samp}"
tophat -i rat_transcriptome/rat_salmon_index -l A \
         -1 ${fn}/${samp}_*_R1_001.fastq.gz \
         -2 ${fn}/${samp}_*_R2_001.fastq.gz \
         -p 8 -o tophat/${samp}_thout -G 

done 

