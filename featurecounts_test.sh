#!/bin/bash
#SBATCH --account=pas2905
#SBATCH --job-name=featurecounts_test
#SBATCH --mem=64gb
#SBATCH --time=40:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=R-%x.%j.out
#SBATCH --error=R-%x.%j.err

module load samtools
module load subread


script dmso6_featureCounts_log.txt
featureCounts -p -O -T 6 -a /users/PAS2905/coraalbers/h9c2_RNA_seq/rat_genome/Rattus_norvegicus/Ensembl/RGSC3.4/Annotation/Genes/genes.gtf -o dmso6_featureCounts_output.txt dmso6.sorted.byname.bam
exit