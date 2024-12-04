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

load('Visualisation factorielle/data_clean.RData')

### On travaille ici sur le dataframe new-full sans les deux maladies complexes  
### n'ayant aucune variabilité (0 pour tous les phénotypes)

# à ne lancer que si nécessaire car très long. Autrement la matrice de distance est importé plus bas
# d_sorensen <- dist.binary(new_full,method = 5)


load('distances.RData')

matrice_sorensen <- as.data.frame(as.matrix(origin_sorensen))

o_sorensen <- matrice_sorensen[6126:7089,1:6125]


o_sorensen_mat <- as.matrix(o_sorensen)

image(1:964, 1:6125, o_sorensen_mat, axes = T, xlab="", ylab="")

# Transform the matrix in long format
df <- melt(o_sorensen_mat)
colnames(df) <- c("y", "x", "value")
# 
# plot_sorensen <- ggplot(df, aes(x = x, y =y, fill = value)) +
#   geom_tile()+
#   theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,size=4),
#         axis.text.y = element_text(size = 4))
# 
# plot_sorensen






