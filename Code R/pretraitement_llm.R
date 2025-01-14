# ------------------------------------------------------------------------------
# Prétraitement avant LLM
# 14/01/2025
# Leslie Cieters
# ------------------------------------------------------------------------------

# rm(list = ls())

## On cherche à automatiser les prompts LLM pour chaque association gène-maladie
## Le but de ce prétraitement est d'obtenir à partir d'une matrice d'assignation,
## un dataframe où chaque ligne correspond à une association

library(tidyverse)
library(reshape)

## Pour le moment nous travaillons sur le RData Associations_1 mais à changer 
## plus tard

load('Data/Associations_1.RData')
df <- associations_seuilpermut_1
colnames(df) <- c('disease','num','gene')

# Transformation du DataFrame
df_tidy <- df %>%
  select(disease,gene) %>% 
  separate_rows(gene, sep = ";") %>% # Sépare les gènes en lignes distinctes 
  mutate(gene = str_remove_all(gene, "\\([^)]*\\)")) %>%  # Enlever le contenu entre parenthèses
  mutate(gene = str_trim(gene))                           # Nettoyer les espaces inutiles

df_tidy <- unique(df_tidy) # retirer les lignes identiques

save(df_tidy,file='Data/asso_gene-maladie.RData')

# apply(df_tidy,1,FUN = function(row){
  # prompt <- sprintf("petit prompt",row[2],row[1],row[2],row[1])})
