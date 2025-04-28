WD=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/phylogenetics/202503_snapp/snapp

#chr_W, run 2 iterations of this with a different SNP sample, N=1000 SNPs 
bcftools view ../vcfs/chr_W.SNP.DP3-AC1-MQ40-MM1.vcf.gz -Ov -o chr_W.SNP.DP3-AC1-MQ40-MM1.vcf

# First rep
ruby snapp_prep.rb --xml rep1_crown.xml -o rep1_crown -v chr_W.SNP.DP3-AC1-MQ40-MM1.vcf -t Females_N93.pop -c constraints_crown.txt -m 1000 -l 100000

# Second rep 
ruby snapp_prep.rb --xml rep2_crown.xml -o rep2_crown -v chr_W.SNP.DP3-AC1-MQ40-MM1.vcf -t Females_N93.pop -c constraints_crown.txt -m 1000 -l 100000

# huge snps 
ruby snapp_prep.rb --xml rep1_crown_30k.xml -o rep1_crown_30k -v chr_W.SNP.DP3-AC1-MQ40-MM1.vcf -t Females_N93.pop -c constraints_crown.txt -m 30000 -l 100000

# First rep, outgroup
ruby snapp_prep.rb --xml rep1_out.xml -o rep1_out -v chr_W.SNP.DP3-AC1-MQ40-MM1.vcf -t Females_N93.pop -c constraints_outgroup.txt -m 1000 -l 100000

# Second rep, outgroup
ruby snapp_prep.rb --xml rep2_out.xml -o rep2_out -v chr_W.SNP.DP3-AC1-MQ40-MM1.vcf -t Females_N93.pop -c constraints_outgroup.txt -m 1000 -l 100000

# Second rep, outgroup, huge snps 
ruby snapp_prep.rb --xml rep1_out_30k.xml -o rep1_out_30k -v chr_W.SNP.DP3-AC1-MQ40-MM1.vcf -t Females_N93.pop -c constraints_outgroup.txt -m 30000 -l 100000
