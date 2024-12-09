# 03/12/2024
# Leslie Cieters
# Matrices de distance depuis la matrice d'origine
# ------------------------------------------------------------------------------

rm(list = ls())

# Packages -------------------------------

library(tidyverse)
library(vegan)
library(ade4) #calcul des distances
library(reshape)

# Chargement des données -----------------

load('Data/data_clean0.RData')

### On travaille ici sur le dataframe new-full sans les deux maladies complexes  
### n'ayant aucune variabilité (0 pour tous les phénotypes)

# à ne lancer que si nécessaire car très long. Autrement la matrice de distance est importé plus bas
dist_or_sorensen <- dist.binary(new_full,method = 5)
dist_or_ochiai <- dist.binary(new_full,method = 7)

dist_or_sorensen_mx <- as.matrix(dist_or_sorensen)
dist_or_ochiai_mx <- as.matrix(dist_or_ochiai)

save(dist_or_sorensen_mx,file="Data/distance_origin_sorensen.RData")
save(dist_or_ochiai_mx,file = "Data/distance_origin_ochiai.RData")

load("Data/distance_origin_sorensen.RData")
load("Data/distance_origin_ochiai.RData")







