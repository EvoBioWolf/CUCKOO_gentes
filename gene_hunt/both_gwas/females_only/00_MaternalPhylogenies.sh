#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=5
#SBATCH --time=24:00:00

#mamba activate snps

#mask with male-biased coverage
mask=/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/Coverage_Masks/GCA_017976375.1_bCucCan1.pri_genomic.CHR.N75-DoubleCoverage.mask.bed
wd=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/gwas_n3/both/females_only
raw_vcfs=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/merged/snps_only/

cd $wd
mkdir -p vcfs ml_trees

for CHR in chr_MT chr_W; do

#filter with NO SINGLETONS
bcftools view vcfs/${CHR}.SNP.DP3.vcf.gz --min-alleles 2 --max-alleles 2 --min-ac 2 --max-af 0.999 --types snps -Oz -o vcfs/${CHR}.SNP.DP3-AC2.vcf.gz
bcftools index --threads 5 vcfs/${CHR}.SNP.DP3-AC2.vcf.gz

#create tree NO SINGLETONS
python ~/modules/vcf2phylip.py -i vcfs/${CHR}.SNP.DP3-AC2.vcf.gz -f --output-folder ml_trees
iqtree --redo -keep-ident -T 5 -s ml_trees/${CHR}.SNP.DP3-AC2.min4.phy --seqtype DNA -m "MFP+ASC" -alrt 1000 -B 1000
iqtree --redo -keep-ident -T 5 -s ml_trees/${CHR}.SNP.DP3-AC2.min4.phy.varsites.phy --seqtype DNA -m "MFP+ASC" -alrt 1000 -B 1000

#create tree WITH SINGLETONS 
python ~/modules/vcf2phylip.py -i vcfs/${CHR}.SNP.DP3.vcf.gz -f --output-folder ml_trees
iqtree --redo -keep-ident -T 5 -s ml_trees/${CHR}.SNP.DP3.min4.phy --seqtype DNA -m "MFP+ASC" -alrt 1000 -B 1000
iqtree --redo -keep-ident -T 5 -s ml_trees/${CHR}.SNP.DP3.min4.phy.varsites.phy --seqtype DNA -m "MFP+ASC" -alrt 1000 -B 1000

done

