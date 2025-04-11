#!/bin/bash
#SBATCH --account=pas2905
#SBATCH --job-name=fastqc_trimmed
#SBATCH --mem=64gb
#SBATCH --time=40:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=R-%x.%j.out
#SBATCH --error=R-%x.%j.err

module load fastqc


for fn in fastq/X126{58..73};
do
samp=`basename ${fn}`
echo "Processing sample ${samp}"
mkdir ${samp}_qc_trimmed
fastqc -o ${samp}_qc_trimmed ${samp}.trimmed_1P.fastq.gz ${samp}.trimmed_2P.fastq.gz 
done 

