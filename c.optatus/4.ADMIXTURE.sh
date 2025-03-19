#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=10
#SBATCH --time=200:00:00

#mamba activate snps
# Submit with for i in {2..10}; do sbatch -J BAD_BOY_SERGIO_${i} 4.ADMIXTURE.sh ${i}; done

wd=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/manyhost_hunt/optatus
cd $wd

K=$1

mkdir -p admixture
cd admixture

#Run Admixture
admixture -j7 --cv=5 ../autosomal_files/autos_optatus_LD.bed ${K} > autos_optatus_LD.log${K}.out

#Evaluate Runs
~/modules/evalAdmix/evalAdmix -plink ../autosomal_files/autos_optatus_LD -fname autos_optatus_LD.${K}.P -qname autos_optatus_LD.${K}.Q -P 10 -o eval_${K}
