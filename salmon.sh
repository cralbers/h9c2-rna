#!/bin/bash
#SBATCH --account=pas2905
#SBATCH --job-name=salmon
#SBATCH --mem=64gb
#SBATCH --time=40:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=R-%x.%j.out
#SBATCH --error=R-%x.%j.err

module load salmon

for fn in fastq/X126{56..73};
do
samp=`basename ${fn}`
echo "Processing sample ${samp}"
salmon quant -i rat_transcriptome/rat_salmon_index -l A \
         -1 ${fn}/${samp}_*_R1_001.fastq.gz \
         -2 ${fn}/${samp}_*_R2_001.fastq.gz \
         -p 8 -o quants/${samp}_quant

done 

