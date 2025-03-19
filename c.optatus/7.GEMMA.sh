#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=5
#SBATCH --time=200:00:00

# Submit as: for EGG in $(cat Eggs.list); do sbatch -J gwas_${EGG} 7.GEMMA.sh ${EGG} ; done 
TARGET=$1

WD=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/gwas_n3/optatus
CHRS=~/EvoBioWolf/CUCKOO_gentes/c.optatus/Chromosomes.list

mkdir -p ${WD}/gemma/${TARGET} mkdir -p ${WD}/gwas
cd ${WD}/gemma/${TARGET}

for CHR in $(cat ${CHRS}); do

if [[ $CHR = 'chr_W' ]]; then
        #grab phenotypes, females only
        awk -v target="$TARGET" '{if ($2 == target) print 1; else print 0}' ../../Females.pop > ${CHR}.phenotypes

        #run gemma
        ~/EvoBioWolf/CUCKOO_gentes/gemma-0.98.6-pre1 -bfile ../../beds/${CHR} -k ../../beds/Covariance_Matrix_N50_females.txt -p ${CHR}.phenotypes -lm 1 -o ${CHR} -miss 0.1 -maf 0.01

else

        #grab phenotypes
        awk -v target="$TARGET" '{if ($2 == target) print 1; else print 0}' ../../Eggs.pop > ${CHR}.phenotypes

        #run gemma
        ~/EvoBioWolf/CUCKOO_gentes/gemma-0.98.6-pre1 -bfile ../../beds/${CHR} -k ../../beds/Covariance_Matrix_N50.txt -p ${CHR}.phenotypes -lm 1 -o ${CHR} -miss 0.1 -maf 0.01

fi

#save output
awk -v t=${TARGET} '{OFS="\t"}{print $1, $3, $3, $11, t}' output/${CHR}.assoc.txt | \
        sed '1d' | \
        bedtools intersect -a - -b ../../raw_dnds/${CHR}.bed -wao | \
        awk '{OFS="\t"}{print $1, $2, $4, $5, $9, $10}' | \
        sed 's/gene-//g' | sed 's/ID=//g' > ../../gwas/${TARGET}_${CHR}.gwas
done
