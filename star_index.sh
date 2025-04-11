#!/bin/bash
#SBATCH --account=pas2905
#SBATCH --job-name=star_rat_index
#SBATCH --time=40:00:00
#SBATCH --mem=16000mb
#SBATCH --mail-type=ALL


cd /users/PAS2905/coraalbers/H9c2_RNA_seq

module load star/2.5.2a

STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir rat_index \
--genomeFastaFiles /users/PAS2905/coraalbers/H9c2_RNA_seq/rat_genome/Rattus_norvegicus/Ensembl/RGSC3.4/Sequence/WholeGenomeFasta/genome.fa \
--sjdbGTFfile /users/PAS2905/coraalbers/H9c2_RNA_seq/rat_genome/Rattus_norvegicus/Ensembl/RGSC3.4/Annotation/Genes/genes.gtf
