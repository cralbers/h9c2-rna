#!/bin/bash
#SBATCH --account=pas2905
#SBATCH --job-name=hisat2_align
#SBATCH --mem=64gb
#SBATCH --time=40:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=R-%x.%j.out
#SBATCH --error=R-%x.%j.err


module load hisat2/2.1.0



for fn in users/PAS2905/coraalbers/h9c2_RNA_seq/fastq/X126{56..73};
do
samp=`basename ${fn}`
echo "Processing sample ${samp}"
hisat2 -x rat_tran \
         -1 ${fn}/${samp}_*_R1_001.fastq.gz \
         -2 ${fn}/${samp}_*_R2_001.fastq.gz \
         -p 8 -S ${samp}_hisat2_align

done 

