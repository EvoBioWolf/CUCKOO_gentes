#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=20
#SBATCH --time=48:00:00

# mamba activate snps

wd=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/gwas_n3/both
cd $wd

mkdir -p autosomal_files 

#merge autosomes
bcftools concat --file-list Autosomes.list --threads 20 -Oz -o autosomal_files/autos.vcf.gz
bcftools index --threads 20 autosomal_files/autos.vcf.gz

#LD prune
~/modules/plink2 --threads 20 --vcf autosomal_files/autos.vcf.gz --allow-extra-chr --set-missing-var-ids @:# \
        --rm-dup --indep-pairwise 50 5 0.1 --maf 0.05 --hwe 1e-10 --max-alleles 2 --min-alleles 2 --out autosomal_files/autos_both_LD
        
#extract, also a vcf and run PCA 
~/modules/plink2 --threads 20 --vcf autosomal_files/autos.vcf.gz --allow-extra-chr --set-missing-var-ids @:# \
        --extract autosomal_files/autos_both_LD.prune.in \
        --make-bed --recode vcf bgz --pca --out autosomal_files/autos_both_LD
bcftools index --threads 20 autosomal_files/autos_both_LD.vcf.gz
sed -i 's/chr_//g' autosomal_files/autos_both_LD.bim
