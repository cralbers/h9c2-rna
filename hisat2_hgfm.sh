#!/bin/bash
#SBATCH --account=pas2905
#SBATCH --job-name=hisat2_hgfm
#SBATCH --mem=200gb
#SBATCH --time=40:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=R-%x.%j.out
#SBATCH --error=R-%x.%j.err





module load hisat2
module load python/3.6-conda5.2

hisat2-build -p 16 --exon Rattus_norvegicus.mRatBN7.2.113.exon --ss Rattus_norvegicus.mRatBN7.2.113.ss Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa rat_tran