#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=2
#SBATCH --time=48:00:00

# mamba activate snps
# for CHR in $(cat Chromosomes.list); do sbatch -J Merge_${CHR} ~/EvoBioWolf/CUCKOO_gentes/gene_hunt/both_gwas/2B.Merge_VCF.sh ${CHR}; done
CHR=$1

bcftools merge ../canorus/vcfs/${CHR}.SNP.DP3.vcf.gz ../optatus/vcfs/${CHR}.SNP.DP3.vcf.gz | \
    bcftools view --force-samples --samples-file AllSamples.list -Ou | \
    bcftools view --min-alleles 2 --max-alleles 2 --min-ac 1 -e 'F_MISSING > 0.1' -Oz -o vcfs/${CHR}.SNP.DP3.vcf.gz
