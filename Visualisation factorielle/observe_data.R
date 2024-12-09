# Leslie Cieters
# 29/11/2024
# Test AFMH
# ------------------------------------------------------------------------------


rm(list=ls())
# Importation des packages
library(FactoMineR)
library(Factoshiny)
library(tidyverse)

# Importation des données
omim <- read.csv('Matrices de Jules/Profils_OMIM.csv.gz')
phecode <- read.csv('Matrices de Jules/Profils_Phecodes.csv.gz')


# Mise en forme et juxtaposition des tableaux de données
colnames(omim)[1] <- "Maladie"
colnames(phecode)[1] <- "Maladie"

omim$Maladie <- as.factor(omim$Maladie)
phecode$Maladie <- as.factor(phecode$Maladie)

data_full <- rbind(omim,phecode)

rownames(data_full) <- data_full[,1]
data_full <- data_full[,-1]


# vérification de la variabilité par phénotype
min(apply(data_full[c(1:6125),],MARGIN = 1,FUN = sum))

max(apply(data_full,MARGIN = 1,FUN = mean))

# # transformation en facteurs
# data_full <- data.frame(lapply(data_full,as.factor))



save(data_full,omim,phecode,new_full,file='data_clean.RData')

# checker pour la variabilité 
# envoyer le jeu de données propre
# faire tourner 

# res.AC <- CA(data_full[,-1],row.sup = c(6126:7090))
class(toto)

# Graphique standard
plot.CA(toto,col.row.sup = 'green',col.row = 'blue',autoLab = 'yes')

# Plot plus détaillé avec options supplémentaires
plot(toto, 
     axes = c(1, 2),          # Choix des axes à représenter
     col.row = "blue",        # Couleur des points de lignes
     col.col = "red",         # Couleur des points de colonnes
     title = "Analyse des Correspondances",
     autolab = "yes")         # Étiquetage automatique

# Options supplémentaires de visualisation
# Graphique des valeurs propres
plot(toto$eig[,1], type = "b", 
     main = "Valeurs propres", 
     xlab = "Dimensions", 
     ylab = "Valeur propre")

# Impression des informations principales
print(toto)

