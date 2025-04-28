#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=25
#SBATCH --time=300:00:00

WD=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/phylogenetics/202503_snapp/1Mchains

cd ${WD}

RUN=$1
~/modules/beast/bin/beast -threads 25 -overwrite -beagle_SSE -seed 777 -java ${RUN}.xml
