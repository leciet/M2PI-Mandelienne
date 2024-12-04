# rm(list = ls())

# 03/12/2024
# Leslie Cieters
# Assignation des maladies 
# -----------------------------------------------------------------------------


library(tidyverse)
library(reshape)


### On travaille pour le moment sur la matrice originel avec la distance de 
### Sorensen

load('distances.RData')

matrice_sorensen <- as.data.frame(as.matrix(origin_sorensen)) #matrice complète


### On récupère la partie d'intérêt dans la matrice complète

o_sorensen <- matrice_sorensen[6126:7089,1:6125] #dataframe de la matrice des distances

o_sorensen_mat <- as.matrix(o_sorensen) #matrice des distances


# Conversion de la matrice en format long
df <- melt(o_sorensen_mat)
colnames(df) <- c("mc", "ms", "distance")

df$ms <- sapply(strsplit(as.character(df$ms), " "), `[`, 1)

## 1. par un seuil global-------------------------------------------------------

### on fixe un seuil s pour la distance et on assigne aux m.c. les m.s. à une 
### distance inférieur ou égale à ce seuil -> utilisation de df

s <- 0.6 # à modifier si besoin 



df_filtre <- df %>% 
  filter(distance<=s) %>%   
  group_by(mc) %>% 
  summarise( mc = mc, liste = paste(ms,collapse = " ; "))


df_filtre <- unique(df_filtre)

str_count(df_filtre$liste,';')

df_filtre[103,]

## 2. par un seuil individuel --------------------------------------------------

q <- 0.005 # fixe le quantile qui servira de limite 

df_filtre <- df %>% 
  group_by(mc) %>% 
  filter(distance <= quantile(distance,q)) %>% 
  summarise( mc = mc, liste = paste(ms,collapse = " ; "),distance=distance,.groups = 'drop')


df_filtre <- unique(df_filtre)

str_count(df_filtre$liste,';')

df_filtre[103,]






