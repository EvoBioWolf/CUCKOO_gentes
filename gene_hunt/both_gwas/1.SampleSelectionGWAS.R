#### Identify samples to use for GWAS, unrelated individuals - but include males 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/gwas_n3')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(viridis)
library(ggpubr)
library(sf)
library(ggspatial)
md <- read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt')
md <- md %>% mutate(Egg = ifelse(Egg == 'E6' & Haplogroup == 'M3','E6M3',Egg))
egg_count <- md %>% drop_na(Egg) %>%
  filter(Analysis_PopulationGenetics == 1) %>% 
  mutate(Egg = ifelse(Egg == 'E6' & Haplogroup == 'M3','E6M3',Egg)) %>% ungroup %>% 
  count(Egg) %>% 
  filter(n > 2)
eggcols = md %>% filter(Egg %in% egg_count$Egg) %>% select(Egg) %>% 
  unique %>% mutate(ord = as.numeric(gsub('M3','',gsub('E','',Egg)))) %>% arrange(ord) %>% 
  drop_na(Egg) %>% mutate(col = viridis(7,option='turbo')) 
targets <- md %>% filter(Egg %in% egg_count$Egg)

#save for all v 1 
allv1 = targets %>% select(ID,Egg)
write.table(allv1,file='AllSamples.pop',quote=F,sep='\t',row.names=F,col.names=F)
write.table(allv1$ID,file='AllSamples.list',quote=F,sep='\t',row.names=F,col.names=F)

#jitter points up to 1 lat/long for viewing
targets = targets %>% mutate(LatJit = jitter(Latitude,amount =2),
                   LonJit = jitter(Longitude,amount=2))
sites = st_as_sf(targets, coords = c("LonJit", "LatJit"), 
                 crs = 4326, agr = "constant")

#set up map and convert df to coordinate frame
world = map_data("world")

gwas_sample_plot = ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), col='grey90', fill='white') +
  geom_sf(data = sites, 
          aes(fill=as.factor(Egg)),
          size=3,show.legend = T,pch=21) +
  xlab('')+ylab('')+
  coord_sf(xlim = c(min(md$Longitude)-5, max(md$Longitude)+5), 
           ylim = c(min(md$Latitude)-5, max(md$Latitude)+5), expand = FALSE)+
  scale_fill_manual(values=eggcols$col,breaks=eggcols$Egg)+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),panel.background = element_rect(fill = "aliceblue"))+
  theme(legend.text = element_text(size = 6),legend.title = element_text(size = 6),legend.key.size = unit(0.1, 'cm'))+
  annotation_scale(line_width=0.5)+
  guides(fill=guide_legend(override.aes=list(shape=21)))
gwas_sample_plot

pdf('../figures/20250301_Canorus_SpatialDistribution-GWAS.pdf',height=4,width=7)
gwas_sample_plot
dev.off()
