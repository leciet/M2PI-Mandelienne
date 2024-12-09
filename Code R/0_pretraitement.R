# Leslie Cieters
# 29/11/2024
# Récupération des données 
# ------------------------------------------------------------------------------


rm(list=ls())


# ------------------------- Importation des packages ---------------------------
library(FactoMineR)
library(Factoshiny)
library(tidyverse)

#----------------------- Importation des données--------------------------------
omim <- read.csv('Data/Profils_OMIM.csv.gz')
phecode <- read.csv('Data/Profils_Phecodes.csv.gz')


#----------- Mise en forme et juxtaposition des tableaux de données-------------
colnames(omim)[1] <- "Maladie"
colnames(phecode)[1] <- "Maladie"

omim$Maladie <- as.factor(omim$Maladie)
phecode$Maladie <- as.factor(phecode$Maladie)

or_df0 <- rbind(omim,phecode)

rownames(or_df0) <- or_df0[,1]
or_df0 <- or_df0[,-1]


# ----------vérification de la variabilité par phénotype------------------------

min(apply(or_df0,MARGIN = 1,FUN = sum)) 
#il existe des mc ne possédant aucun des phénotypes

max(apply(or_df0,MARGIN = 1,FUN = sum))


save(or_df0,omim,phecode,file='Data/data_clean0.RData')

# checker pour la variabilité 


or_filter_df0 <- or_df0 %>% 
  filter(rowSums(or_df0)!=0)


