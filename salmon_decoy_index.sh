#!/bin/bash
#SBATCH --account=pas2905
#SBATCH --job-name=salmon-decoy-index
#SBATCH --mem=64gb
#SBATCH --time=40:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=R-%x.%j.out
#SBATCH --error=R-%x.%j.err


module load salmon/1.2.1



salmon index -t rat_gentrome.fa.gz --decoys decoys.txt -p 12 -i rat_salmon_index
