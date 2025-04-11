#!/bin/bash
#SBATCH --account=pas2905
#SBATCH --job-name=star-all
#SBATCH --mem=64gb
#SBATCH --time=40:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=R-%x.%j.out
#SBATCH --error=R-%x.%j.err

module load star/2.5.2a


for fn in fastq/X126{58..73};
do
samp=`basename ${fn}`
echo "Processing sample ${samp}"
mkdir STAR/${samp}
STAR --genomeDir /users/PAS2905/coraalbers/h9c2_RNA_seq/rat_star_index \
--runThreadN 6 \
--readFilesIn ${fn}/${samp}_*_R1_001.fastq.gz ${fn}/${samp}_*_R2_001.fastq.gz \
--readFilesCommand zcat \
--outFileNamePrefix STAR/${samp}/${samp}

done 

