#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=20
#SBATCH --time=8:00:00

CHR=$1

# Submit: for CHR in $(cat Chromosomes.list); do sbatch -J dnds_${CHR} 3B.Submit_Annotation.sh ${CHR}; done
# mamba activate R
Rscript 3.Annotate_Variants.R ${CHR}

sed '1d' raw_dnds/Annotated_Variants_${CHR}__20250301.txt | awk '{OFS="\t"}{print $1, $2, $2, $5, $6}' | sed 's/gene-//g' | sed 's/ID=//g' > raw_dnds/${CHR}.bed

