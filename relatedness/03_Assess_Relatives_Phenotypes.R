#### Determine relatedness  
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/merged_full')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(viridis)
library(adegenet)
library(vcfR)
library(ggpubr)
library(scales)

# Analyze LD-pruned relatedness from vcftools --relatedness2, KING coefficient 
rel <- read.table('autos_canorus_LD.relatedness2',header=T) %>% select(1,2,7) %>% as_tibble
names(rel) <- c('ID_A','ID_B','PHI')
rel = rel %>% mutate(PHI = pmax(0, pmin(0.5, PHI)))
rel = rel %>% filter(ID_A != ID_B)

# Compare with plink IBS0
pk = read.table('autos_canorus_LD.ibs',sep=' ',header=TRUE); names(pk) = c('ID_A','ID_B','IBS0')
relpk = left_join(pk,rel)

# Values from here for designating relationships https://www.kingrelatedness.com/manual.shtml
fam = relpk %>% mutate(Relationship = ifelse(PHI > 0.354, 'First Degree',
                                             ifelse(PHI > 0.177, 'First Degree',
                                                    ifelse(PHI > 0.0884,'Second Degree',
                                                           ifelse(PHI > 0.0442, 'Third Degree',
                                                                  ifelse(PHI <= 0.0442,'Unrelated','Unassigned'))))))

fam %>% filter(PHI > 0.354) #usually > 0.354 is MZ twin / duplicate, but since our highest value is 0.393 and most around 0.37, seems more likely they are just first degree

# For visualizing KING vs IBS
# phi_ibs = fam %>% ggplot(aes(x=PHI,y=IBS0,col=Relationship))+
#   geom_point()+
#   scale_color_viridis(discrete=TRUE)+
#   theme_bw()
# pdf('../figures/Relatedness__PhivIBS_2023OCT27.pdf',height=5,width=6)
# phi_ibs
# dev.off()

# Add mtDNA differences, calculated externally which identifies pairwise SNPs between samples for mtDNA
mt = read.table('chrMT_Pairwise_Differences.txt')
names(mt) = c('ID_A','ID_B','Sites','SNPs')
# The script counted diploid genotypes as 2 SNPs, so divide by 2
mt <- mt %>% mutate(SNPs = SNPs / 2)
fam2 = left_join(fam,mt) %>% drop_na(SNPs)

# Remove redundant comparisons, e.g. ID_A vs ID_B or pairwise redundancies
famrm = fam2 %>% 
  select(-IBS0) %>% 
  rowwise() %>% 
  mutate(pair = sort(c(ID_A,ID_B)) %>% paste(collapse = ",")) %>%
  group_by(pair,Sites,SNPs) %>%
  distinct(pair, .keep_all = T) %>% 
  separate(pair,into=c('ID_A','ID_B'),remove=F,sep=',') %>% ungroup() %>% select(-pair)

# Merge with metadata
md = read_tsv('~/merondun/cuculus_host/relatedness/Full_Metadata.txt')
md = md %>% select(c(ID,HostParentShort,Habitat,Egg,KDist,Sampling_Year,Hap,Sex,Age))
mda = md

# Add an '_A' and '_B' to each metadata field 
names(mda) = paste0(names(mda),'_A')
mdb = md
names(mdb) = paste0(names(mdb),'_B')
fam3 = left_join(famrm,mda) %>% 
  left_join(.,mdb)

# Mislabeled sample, exclude
fam3 <- fam3 %>% filter(!grepl('148_',ID_A) & !grepl('148_',ID_B))

# Inspect, look to see which haps don't align, etc 
fam3 %>% filter(Relationship == 'First Degree' & Hap_A != Hap_B) %>% data.frame
fam3 %>% filter(Relationship == 'First Degree' & KDist_A != KDist_B) %>% data.frame
fam3 %>% filter(Relationship == 'First Degree' & Egg_A != Egg_B) %>% data.frame

# Only retain those with known egg morph, which are NESTLINGS, and which are caught within 2 years ( no parent / sibling)
famd = fam3 %>% drop_na(Egg_A) %>% drop_na(Egg_B)
famd = famd %>% filter(KDist_A == KDist_B) # Only compare within the same distance clade
famd = famd %>% filter(Age_A == 'Young' & Age_B == 'Young') # Only compare within YOUNG NESTLINGS!
famd %>% count(Age_A,Age_B)
famd = famd %>% filter(abs(Sampling_Year_A - Sampling_Year_B) <= 2) # Only compare within the same sampling year (e.g. SIBLINGS)
famd %>% count(Age_A,Age_B) #n = 2806

# Within first degree relatives, what's the distribution of the # of mtDNA snps? 
famd %>% filter(Relationship == 'First Degree') %>% 
  ggplot(aes(x=SNPs))+geom_histogram(show.legend = F)+theme_bw()+
  scale_x_continuous(breaks=pretty_breaks(n=14))

# If mtDNA haplotypes are identical, assign as maternal. 
famd = famd %>% mutate(Line = ifelse(SNPs <= 1,'Maternal','Paternal'))
famd %>% filter(Relationship == 'First Degree') %>%  count(Line)
famd %>% filter(Relationship != 'Unrelated') %>% count(Line)

#simply assign unrelated as paternal, move pie chart in final plot 
famd = famd %>% mutate(Line = ifelse(Relationship == 'Unrelated','Paternal',Line))

# and remove first degree relatives, OR NOT! 
#famd = famd %>% filter(Relationship != 'First Degree')

famd %>% filter(Relationship != 'Unrelated') %>% count(Relationship,Line)
#Without first
# Line         n
# <chr>    <int>
#   1 Maternal    63
# 2 Paternal    97

#With first
# Line         n
# <chr>    <int>
#   1 Maternal   128
# 2 Paternal   147

# Relationship  Line         n
# <chr>         <chr>    <int>
#   1 First Degree  Maternal    65
# 2 First Degree  Paternal    50
# 3 Second Degree Maternal    25
# 4 Second Degree Paternal    42
# 5 Third Degree  Maternal    38
# 6 Third Degree  Paternal    55

nrow(famd)
#[1] 2691 with no first degree, [1] 2806 with first degree, or [1] 2948 with no age and year restrictions
#how many unique individuals?
length(unique(c(famd$ID_A,famd$ID_B)))

#function to get matched data
get_matched_data <- function(df) {
  df %>% 
    mutate(
      Host = ifelse(HostParentShort_A == HostParentShort_B, 'Matched', 'Unmatched'),
      Habitat = ifelse(Habitat_A == Habitat_B, 'Matched', 'Unmatched'),
      Year = ifelse(abs(Sampling_Year_A - Sampling_Year_B) <= 2, 'Matched', 'Unmatched'),
      Egg = ifelse(Egg_A == Egg_B, 'Matched', 'Unmatched'),
      Distance = ifelse(KDist_A == KDist_B, 'Matched', 'Unmatched'),
      Haplogroup = ifelse(Hap_A == Hap_B, 'Matched', 'Unmatched')
    ) %>% 
    gather(key = "Variable", value = "Matched", Host, Habitat, Year, Egg, Distance,Haplogroup) %>% 
    group_by(Relationship,Line,Variable) %>% 
    count(Matched)
}

#Divide data into frames and store in a list
data_frames <- list(
  fem = famd %>% filter(Sex_A == 'F' & Sex_B == 'F'),
  mal = famd %>% filter(Sex_A == 'M' & Sex_B == 'M'),
  fm = famd %>% filter(Sex_A != Sex_B)
)

#function to each data frame and store results in a list
results_list <- lapply(names(data_frames), function(frame_name) {
  result_df <- get_matched_data(data_frames[[frame_name]])
  result_df$frame <- frame_name
  result_df
})

#combine
mm <- bind_rows(results_list)
mm$Relationship <- factor(mm$Relationship,levels=c('First Degree','Second Degree','Third Degree','Unrelated'))
mm <- mm %>%
  mutate(frame = case_when(
    frame == "fem" ~ "Female-Female",
    frame == "mal" ~ "Male-Male",
    frame == "fm" ~ "Intersexual",
    TRUE ~ frame)) %>% dplyr::rename(Count = n)

#proportions for pie charts
piece = mm %>% select(-frame) %>% group_by(Relationship,Line,Variable,Matched) %>% na.omit %>% summarize(Count = sum(Count)) %>% ungroup %>% group_by(Relationship,Line,Variable) %>% 
  mutate(Total = sum(Count),
         Proportion = Count/Total,
         Percent = paste0(round(Count/Total,3)*100,'% (',Total,')'))
piece = piece %>% filter(!grepl('Year|Distance|Hap',Variable))

#plot
relpie = piece %>% 
  ggplot(aes(x="",y=Proportion,fill=Matched))+
  geom_bar(stat='identity')+
  coord_polar("y", start=0)+xlab('')+ylab('')+
  facet_grid(Relationship~Variable+Line)+
  scale_fill_manual(values=c('forestgreen','grey95','black'))+
  theme_bw(base_size=7)+labs(fill='Phenotype Matching')+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
#add labels showing the % matched, and then the total sample size 
labs = piece %>% 
  filter(Matched == 'Matched')
piesR = relpie + 
  geom_label(data = labs, 
             aes(y=Inf,x=-Inf,label=Percent),fill='white',
             size=2,vjust=.3,hjust=.5,alpha=0.8)
piesR

pdf('../figures/20241210-Relatives_Phenotype_Matching-1Mismatch-FIRSTDEGREE.pdf',height=4.5,width=6)
piesR
dev.off()

famd %>% filter(Relationship != 'Unrelated' & Egg_A != Egg_B & Line == 'Paternal') %>%
  count(Relationship,Egg_A,Egg_B)

famd %>% filter(Relationship != 'Unrelated' & Egg_A == Egg_B & Line == 'Maternal') %>% 
  count(Relationship,Egg_A,Egg_B)

famd %>% count(Relationship,Egg_A,Egg_B,Line)

## Aggregate
df_agg <- famd %>% 
  count(Relationship,Egg_A,Egg_B,Line) %>% 
  mutate(
    Egg = pmin(Egg_A, Egg_B),  # ensures consistent grouping
    Match = ifelse(Egg_A == Egg_B, n, 0),
    Unmatch = ifelse(Egg_A != Egg_B, n, 0)
  ) %>%
  group_by(Relationship, Egg, Line) %>%
  summarise(
    Match = sum(Match),
    Unmatch = sum(Unmatch),
    .groups = "drop"
  )

# Supplementary table 
df_agg %>% data.frame

# Plot matched 
df_agg %>% 
  filter(Relationship != 'Unrelated') %>% 
  pivot_longer(c(Match,Unmatch)) %>% 
  ggplot(aes(y=Egg,x=value,fill=name))+
  xlab('Count')+
  geom_bar(col='black',stat='identity',position=position_dodge(width=0.5))+
  facet_grid(Relationship~Line,scales='free',space='free_y')+
  theme_bw()

# Perform Fisher's test for each Relationship category
results <- df_agg %>%
  group_by(Relationship) %>%
  summarise(
    p_value = list(
      fisher.test(matrix(c(sum(Match[Line == "Maternal"]),
                           sum(Unmatch[Line == "Maternal"]),
                           sum(Match[Line == "Paternal"]),
                           sum(Unmatch[Line == "Paternal"])),
                         nrow = 2))$p.value
    )
  ) %>%
  unnest(p_value)
print(results)
# Relationship   p_value
# <chr>            <dbl>
#   1 First Degree  0.435   
# 2 Second Degree 0.525   
# 3 Third Degree  0.000118
# 4 Unrelated     1   


#also add missingness, we will remove individual with more missing data if both are females or both are males 
miss <- read.table('autos_canorus_LD.imiss',head=T) %>% select(c(1,5))
names(miss) = c('ID_A','Missing_A')
miss$Missing_A = as.numeric(miss$Missing_A)
miss = miss %>% na.omit
fam4 = fam3 %>% left_join(.,miss) %>% 
  merge(.,miss %>% dplyr::rename(ID_B=ID_A,Missing_B=Missing_A)) 

#Remove unrelateds
reli = fam4 %>% filter(Relationship != 'Unrelated')
relatives = fam4

#Unique samples to test 
samples <- unique(c(reli$ID_A,reli$ID_B))
rm <- NULL 
choice <- reli
#Re-run this multiple times, until it stops removing any individuals. 
for (run in seq(1,5,1)) { 
  
  cat('Running ',run,'\n')
  #grab one sample at a time 
  for (samp in samples) {
    #grab all the records with this sample 
    sf.a <- choice[grepl(samp,choice$ID_A),] 
    sf.b <- choice[grepl(samp,choice$ID_B),] 
    sf <- rbind(sf.a,sf.b) %>% arrange(desc(PHI))
    if(nrow(sf) < 1) next
    #if one if male remove, otherwise, if sampled clade is higher, remove, otherwise if missing data is higher, remove, otherwise random (sample by ID# higher#)
    sf1 <- sf %>% mutate(Remove = ifelse(Sex_A == 'F' & Sex_B == 'M', ID_B,
                                         ifelse(Sex_A == 'M' & Sex_B == 'F', ID_A,
                                                ifelse(Missing_A > Missing_B, ID_A,
                                                       ifelse(Missing_A < Missing_B, ID_B,
                                                              'Problem')))))
    sf1 %>% select(contains(c('ID')))
    #take the ID from the sample to remove
    bad <- head(na.omit(sf1$Remove),1) 
    #unless there are no samples to remove.. 
    if(length(bad) < 1) next
    #remove that bad sample from the pool
    cat('Removing bad sample: ',bad,'\n')
    choice <- choice[!grepl(bad,choice$ID_A),]
    choice <- choice[!grepl(bad,choice$ID_B),]
    #and then restart the whole process 
    rm <- rbind(rm,bad)
    rm(bad)
  }
}

#these are the bad individuals we will remove
rms <- as.data.frame(rm)
names(rms) <- 'Remove'
row.names(rms) <- NULL

#remove those bad IDs from the full dataset, make sure to filter form both ID_A and ID_B columns 
keep <- relatives %>% filter(!ID_A %in% rms$Remove) %>% filter(!ID_B %in% rms$Remove)
#see how many individuals were in the full dataset (sanity check), and then how many are retained after filtering 
length(unique(c(relatives$ID_A,relatives$ID_B)))
length(unique(c(keep$ID_A,keep$ID_B)))
length(unique(c(relatives$ID_A,relatives$ID_B))) - length(unique(rms$Remove))

write.table(unique(c(keep$ID_A,keep$ID_B)),'Unrelated_2023OCT27.list',quote=F,row.names=F,sep='/t',col.names=F)
samps = read_tsv('Unrelated_2023OCT27.list',col_names = F)

#calculate pairwise geographic distance
md = read_tsv('../Cuckoo_Full_Metadata_2023OCT3.txt')
ll = md %>% select(ID,Latitude,Longitude)
ll1 = geosphere::distm(ll %>% select(Longitude,Latitude)) %>% as.data.frame
names(ll1) = ll$ID
ll1$ID = ll$ID
ll2 = ll1 %>% pivot_longer(!(ID),names_to='ID_B',values_to = 'GDistance') %>% dplyr::rename(ID_A=ID)
ll2 = ll2 %>% mutate(GDistance = GDistance/1000)

#Add genetic distance
dp = left_join(fam4,ll2)
