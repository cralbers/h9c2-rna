#!/bin/bash
#SBATCH --account=pas2905
#SBATCH --job-name=trimmomatic
#SBATCH --mem=64gb
#SBATCH --time=40:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=R-%x.%j.out
#SBATCH --error=R-%x.%j.err

module load trimmomatic


for fn in fastq/X126{56..73};
do
samp=`basename ${fn}`
echo "Processing sample ${samp}"
java -jar $TRIMMOMATIC \
    PE -trimlog ${fn}/${samp}_trimlog.txt \
    -summary ${fn}/${samp}_trimsummary.txt \
    ${fn}/${samp}_*_R1_001.fastq.gz \
    ${fn}/${samp}_*_R2_001.fastq.gz \
    -baseout ${fn}/${samp}.trimmed.fastq.gz \
    ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 \
    SLIDINGWINDOW:4:15

done 

