# ------------------------------------------------------------------------------
# 18/12/2024
# Side Quest
# Préparation 'text to vec' pour Jean et Massamba 
# ---
# On a besoin d'un jeu de données contenant une colonne, chaque ligne de la 
# colonne correspond à une maladie et l'ensemble des ses phénotypes observés
# ------------------------------------------------------------------------------

# nettoyage environnement ------------------------------------------------------

rm(list = ls())

# Packages ---------------------------------------------------------------------

library(tidyverse)
library(reshape)

# Données ----------------------------------------------------------------------

load('Data/data_clean0.RData') # on travaille avec or_df0

# Code -------------------------------------------------------------------------

df_long <- melt(as.matrix(or_df0)) # transformation au format long 

colnames(df_long) <- c('maladie','phenotype','presence')

liste_phenotype <- df_long %>% 
  filter(presence == 1) %>% 
  group_by(maladie) %>% 
  summarise( liste_phenotype = paste(maladie,paste(phenotype,collapse = " "),sep = ' : '),.groups = 'drop')

liste_phenotype <- unique(liste_phenotype)[,2]


# Sauvegarder le jeu de données ------------------------------------------------

write_csv(x = liste_phenotype,file = 'Data/data_text2vec.csv')
