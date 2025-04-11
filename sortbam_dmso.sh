#!/bin/bash
#SBATCH --account=pas2905
#SBATCH --job-name=bam_sort_dmso_byname_fc
#SBATCH --mem=64gb
#SBATCH --time=40:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=R-%x.%j.out
#SBATCH --error=R-%x.%j.err


module load samtools
module load subread


cd dmso2_star
samtools sort -m 1000000000 -o dmso2.sorted.byname.bam -n dmso2Aligned.out.bam 
script dmso2_featureCounts_log.txt
featureCounts -p -O -T 6 -a /users/PAS2905/coraalbers/h9c2_RNA_seq/rat_genome/Rattus_norvegicus/Ensembl/RGSC3.4/Annotation/Genes/genes.gtf -o dmso2_featureCounts_output.txt dmso2.sorted.byname.bam
exit
cd ../

cd dmso3_star
samtools sort -m 1000000000 -o dmso3.sorted.byname.bam -n dmso3Aligned.out.bam 
script dmso3_featureCounts_log.txt
featureCounts -p -O -T 6 -a /users/PAS2905/coraalbers/h9c2_RNA_seq/rat_genome/Rattus_norvegicus/Ensembl/RGSC3.4/Annotation/Genes/genes.gtf -o dmso3_featureCounts_output.txt dmso3.sorted.byname.bam
exit
cd ../

cd dmso4_star
samtools sort -m 1000000000 -o dmso4.sorted.byname.bam -n dmso4Aligned.out.bam 
script dmso4_featureCounts_log.txt
featureCounts -p -O -T 6 -a /users/PAS2905/coraalbers/h9c2_RNA_seq/rat_genome/Rattus_norvegicus/Ensembl/RGSC3.4/Annotation/Genes/genes.gtf -o dmso4_featureCounts_output.txt dmso4.sorted.byname.bam
exit
cd ../

cd dmso5_star
samtools sort -m 1000000000 -o dmso5.sorted.byname.bam -n dmso5Aligned.out.bam 
script dmso5_featureCounts_log.txt
featureCounts -p -O -T 6 -a /users/PAS2905/coraalbers/h9c2_RNA_seq/rat_genome/Rattus_norvegicus/Ensembl/RGSC3.4/Annotation/Genes/genes.gtf -o dmso5_featureCounts_output.txt dmso5.sorted.byname.bam
exit
cd ../

cd dmso6_star
samtools sort -m 1000000000 -o dmso6.sorted.byname.bam -n dmso6Aligned.out.bam 
script dmso6_featureCounts_log.txt
featureCounts -p -O -T 6 -a /users/PAS2905/coraalbers/h9c2_RNA_seq/rat_genome/Rattus_norvegicus/Ensembl/RGSC3.4/Annotation/Genes/genes.gtf -o dmso6_featureCounts_output.txt dmso6.sorted.byname.bam
exit
cd ../
