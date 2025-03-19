#### Plot the p-distance matrix on the subset GWAS individuals 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/gwas_n3/optatus')
.libPaths('~/mambaforge/envs/R/lib/R/library')
library(tidyverse)
library(viridis)
library(ggpubr)
library(meRo)
library(RColorBrewer)
library(ggpubr)
library(scales)

md = read_tsv('metadata_n50_optatus.txt')
eggcols = md %>% ungroup %>% select(Egg, col = HostColor) %>% unique 

#load in autosomal distance, calculated: ~/modules/VCF2Dis-1.50/bin/VCF2Dis -InPut autos_optatus_N50_LD.vcf.gz -OutPut autos_optatus_N50_LD.pdist 
auto = read.table('autosomal_files/autos_optatus_N50_LD.pdist',header=F)
#add proper names, because VCF2DIS truncates them
auto = auto %>% mutate(IDNum = gsub('_.*','',V1))
auto_id = left_join(auto,md %>% select(ID) %>% mutate(IDNum = gsub('_.*','',ID)))  
ord = auto_id$ID
autos = auto_id %>% select(!c(IDNum,ID,V1))
names(autos) = ord
rownames(autos) = ord

#the order we want, from the bed file 
ids = read.table('beds/chr_29.fam')

#extract columns in same order as geographic distance
common_names = intersect(rownames(autos),ids$V1)
auto_aligned = autos[common_names, common_names] #subset 'auto' matrix to match the 'geo_mat' matrix rows and columns
auto_mat = as.matrix(auto_aligned)

#visualize autosomal distance matrix in terms of ancestry K = 5
auto_mds = as.data.frame(cmdscale(auto_mat,2))
auto_mds$ID = rownames(auto_mds)
auto_mds_input = left_join(auto_mds,md)

kcols = auto_mds_input %>% select(AncestryA5) %>% unique %>% mutate(col = brewer.pal(4,'Set2'),
                                                                    shape = c(21,22,24,25))


auto_p = ggplot(auto_mds_input,aes(x=V1,y=V2,fill=AncestryA5,shape=AncestryA5)) + 
  geom_point(size=1.5) + theme_bw() + xlab('MDS1')+ylab('MDS2')+ggtitle('Autosomal Distance')+
  scale_shape_manual(values=kcols$shape,breaks=kcols$AncestryA5)+
  scale_fill_manual(values=kcols$col,breaks=kcols$AncestryA5)+
  theme(legend.position = "top", # Moves legend to top
        legend.text = element_text(size = 6), # Adjusts text size
        legend.title = element_text(size = 6))
auto_p

pdf('../../figures/20250301_GWAS_Covariance_Matrix_Coptatus.pdf',height=3,width=3)
auto_p
dev.off()

write.table(auto_mat,file='beds/Covariance_Matrix_N50.txt',quote=F,sep=' ',row.names=F,col.names=F)

##### Females
#the order we want, from the bed file 
ids = read.table('beds/chr_W.fam')

#extract columns in same order as geographic distance
common_names = intersect(rownames(autos),ids$V1)
auto_aligned = autos[common_names, common_names] #subset 'auto' matrix to match the 'geo_mat' matrix rows and columns
auto_mat = as.matrix(auto_aligned)

#visualize autosomal distance matrix in terms of ancestry K = 5
auto_mds = as.data.frame(cmdscale(auto_mat,2))
auto_mds$ID = rownames(auto_mds)
auto_mds_input = left_join(auto_mds,md)
chr_p = ggplot(auto_mds_input,aes(x=-V1,y=V2,fill=AncestryA5,shape=AncestryA5)) + 
  geom_point(size=1.5) + theme_bw() + xlab('MDS1')+ylab('MDS2')+ggtitle('chrW Distance')+
  scale_shape_manual(values=kcols$shape,breaks=kcols$AncestryA5)+
  scale_fill_manual(values=kcols$col,breaks=kcols$AncestryA5)+
  theme(legend.position = "top", # Moves legend to top
        legend.text = element_text(size = 6), # Adjusts text size
        legend.title = element_text(size = 6))
chr_p

pdf('../../figures/20250301_GWAS_Covariance_Matrix_Coptatus.pdf',height=2.5,width=5)
ggarrange(auto_p,chr_p)
dev.off()

write.table(auto_mat,file='beds/Covariance_Matrix_N50_females.txt',quote=F,sep=' ',row.names=F,col.names=F)

