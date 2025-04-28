#!/usr/bin/env rscript
args <- commandArgs(trailingOnly = TRUE)

# Egg associations, MNLR, Model Selection c optatus
#### Find associations between MT/W haplotypes and features
setwd('~/EvoBioWolf/CUCKOO_gentes/gens_associations/mnlr/')
.libPaths('~/mambaforge/envs/rfs/lib/R/library')
library(tidyverse)
library(viridis)
library(RColorBrewer)
library(vegan)
library(FactoMineR)
library(factoextra)
library(ggrepel)
library(ggpubr)
library(ggdist)
library(nnet)
library(purrr)
library(broom)
library(meRo) #devtools::install_github('merondun/meRo',force=TRUE,upgrade='never')
library(caret)

sp = args[1]

set.seed(123)
md = read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt')

#only grab samples with known egg (hash out filter for related individuals)
md_egg = md %>%
  filter(SpeciesShort == sp & Analysis_FullGensAssociations == 1) %>%
  select(ID = ID, Egg, Haplogroup, Geography = GeographicGroup,Ancestry = AncestryA5)
md_egg %>% count(Egg)

md_egg = md_egg %>% group_by(Egg) %>% mutate(TotalEgg = n()) %>% ungroup
minobs=2
md_egg = md_egg %>% filter(TotalEgg >= minobs)
md_egg %>% count(Egg)

##### Initialize, summary stats on raw data #####
#ensure they are all factor variables
md_egg = md_egg %>% mutate_at(c('Egg','Haplogroup','Ancestry','Geography'),as.factor)

#### Count proportions first, count proportions for Egg and habitat and egg
ep = md_egg %>% group_by(Egg) %>% mutate(Total = n()) %>%
  pivot_longer(c(Ancestry,Geography,Haplogroup)) %>% ungroup %>%
  group_by(Egg,name,value) %>%
  summarize(Parts = n(),
            Proportion = Parts / Total,
            Percent = paste0(round( Parts / Total,3)*100,'% (',Total,')')) %>% unique

# Bind them together
ap = rbind(ep %>% ungroup %>% mutate(Response = Egg, variable = 'Egg') %>% select(-Egg))

# Just for ordering the covariates nicely
# Plot proportions
ord <- ap %>% select(name, value) %>%
  distinct() %>%
  mutate(ord = as.numeric(gsub("[^0-9.]", "", value))) %>%
  arrange(name, ord)
ap$value <- factor(ap$value,levels=ord$value)
egglev <- md %>% select(Egg,EggOrder) %>% distinct %>% arrange(EggOrder)
ap$Response <- factor(ap$Response,levels=egglev$Egg)
app = ap %>%
  arrange(value) %>%
  ggplot(aes(y=Response,x=value,fill=Proportion,label=Parts))+
  facet_grid(variable~name,scales='free',space='free')+
  scale_fill_gradient('Proportion Observed',low='yellow',high='red')+
  geom_tile()+
  geom_text(size=1)+
  ylab('')+xlab('')+
  theme_bw(base_size=6)+
  theme(legend.position = 'top')
app

ggsave(paste0('~/symlinks/host/figures/20250330_',sp,'_MNLR_Proportions.pdf'),app,height=3,width=7,dpi=300)

##### Model Selection #####
#assess covariate importance with model selection, using MNLR
vars = 'Egg'

# With bootstrapping
ctrl = trainControl(method = "boot",   # Use bootstrapping
                    number = 100,      # Replicates
                    summaryFunction = multiClassSummary,  #caret attempts to stratify levels
                    classProbs = TRUE, #For classification models, to save class probabilities
                    savePredictions = TRUE) #save predictions for each fold

# Change punctuation e.g. 'A. pal' to A_pal' for host fork
md_cv = md_egg %>% mutate(Egg = gsub('\\. ','_',Egg))

# For determining which predictors improve model:
adat = NULL; conf_dat = NULL; class_dat = NULL; new_preds = NULL; counter = 0

# Filter at the very base level, ensuring that across egg / Egg / habitat we have the same individuals with representative sampling
md_subbed = md_cv %>% mutate_at(c('Egg','Haplogroup','Ancestry','Geography'),as.factor)

for (rep in seq(1,100,1)){  # Create 100 replicate models
  for (var in vars) { counter = counter + 1;
  
  # Ensure that we have adequate levels, only a sanity since it is already filtered
  retained = md_subbed %>% ungroup %>% count(!!sym(var)) %>% filter(n >= 2) %>% pull(var)
  length(retained)
  # Number of samples to subsample
  subsamp = md_subbed %>% ungroup %>% count(!!sym(var)) %>% filter(n >= 2) %>% summarize(min = min(n)) %>% pull(min)
  
  cat('Downsampling to n = ',subsamp,', requiring min = ',2,' for variable: ',var,', ','replicate: ',rep,'\n')
  # Subsampling
  mdi = md_subbed %>%
    filter(!!sym(var) %in% retained) %>%
    group_by(!!sym(var)) %>%
    sample_n(min(n(), subsamp),replace = TRUE) %>%
    ungroup()
  mdi %>% count(!!sym(var))
  
  # Ensure the factor levels are dropped which aren't in the dataset after filtering
  mdi = droplevels(mdi)
  
  # First MNLR on combinations
  formula_1 = as.formula(paste(var, "~ Haplogroup + Ancestry + Geography"))
  m1 = train(formula_1, data = mdi, method = "multinom", trControl = ctrl, metric = "AUC", trace = FALSE)
  
  formula_2 = as.formula(paste(var, "~ Haplogroup + Geography"))
  m2 = train(formula_2, data = mdi, method = "multinom", trControl = ctrl, metric = "AUC", trace = FALSE)
  
  formula_3 = as.formula(paste(var, "~ Haplogroup "))
  m3 = train(formula_3, data = mdi, method = "multinom", trControl = ctrl, metric = "AUC", trace = FALSE)
  
  formula_4 = as.formula(paste(var, "~ Haplogroup + Ancestry"))
  m4 = train(formula_4, data = mdi, method = "multinom", trControl = ctrl, metric = "AUC", trace = FALSE)
  
  formula_5 = as.formula(paste(var, "~ Ancestry"))
  m5 = train(formula_5, data = mdi, method = "multinom", trControl = ctrl, metric = "AUC", trace = FALSE)
  
  formula_6 = as.formula(paste(var, "~ Ancestry + Geography"))
  m6 = train(formula_6, data = mdi, method = "multinom", trControl = ctrl, metric = "AUC", trace = FALSE)
  
  formula_7 = as.formula(paste(var, "~ Geography"))
  m7 = train(formula_7, data = mdi, method = "multinom", trControl = ctrl, metric = "AUC", trace = FALSE)
  
  models = c('m1','m2','m3','m4','m5','m6','m7')
  
  # Extract model fit
  for (model in models) {
    # Output model fit from confusion matrix
    mo = get(model)
    
    # Get AIC
    final_model = mo$finalModel;
    AIC = AIC(final_model)
    
    # Save the model results
    dat = data.frame(Model = model, Iteration = rep, Variable = var, Subsampled = subsamp, MinObs = 2,
                     decay = mo$results$decay, logLoss = mo$results$logLoss, logLossSD = mo$results$logLossSD, AIC = AIC, AUC = mo$results$AUC, Accuracy = mo$results$Accuracy, AccuracySD = mo$results$AccuracySD, Kappa = mo$results$Kappa)
    dat_best = dat %>% slice_min(logLoss)
    adat = rbind(adat,dat_best)
    
    # Also save training confusion matrix
    pred = mo$pred
    pred_best = pred %>% filter(decay == dat_best$decay)
    
    # Predict against real data
    predicted_classes = predict(mo, newdata = mdi, type = "raw")
    new = mdi %>% select(reference = !!sym(var)) %>% mutate(predicted = predicted_classes)
    
    conf_new = confusionMatrix(new$predicted, new$reference)
    conf_real = as.data.frame(conf_new$table) %>% mutate(
      Prediction = gsub('_','\\. ',Prediction),
      Reference = gsub('_','\\. ',Reference)) %>%
      mutate(Model = model, Iteration = counter, Variable = var, Subsample = subsamp, MinObs = 2, AUC=dat_best$AUC,logloss=dat_best$logLoss,Accuracy = dat_best$Accuracy,AccuracySD=dat_best$AccuracySD)
    new_preds = rbind(new_preds,conf_real)
    rm(conf_real,dat,dat_best)
    
  } # Exit model loop
  } # Exit response variable loop
} # Exit iteration loop

# Write the no-blue data
write.table(adat,paste0('20250330_Model_Selection_Boot-2Obs-100Reps_',sp,'.txt'),quote=F,sep='\t',row.names=F)
write.table(new_preds,paste0('20250330_ConfusionMatrix_Boot-2Obs-100Reps_',sp,'.txt'),quote=F,sep='\t',row.names=F)


