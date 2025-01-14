# Filter out rows with no variability
or_filter_df0 <- or_df0 %>% 
  filter(rowSums(or_df0)!=0)

phenotype_maladie_s_c <- or_filter_df0

# Save preprocessed data
save(phenotype_maladie_s_c,  file='Workflow final/phenotype_maladie_s_c.RData')

matrice_distance_Jules <- read.csv("distance_matrix_dist_Jules.csv")
library(dplyr)
matrice_distance_Jules_filtered <- matrice_distance_Jules %>%
  select(-c(Acute.bronchospasm, 
            Nonspecific.abnormal.findings.on.radiological.and.other.examination.of.other.intrathoracic.organs..echocardiogram..etc.))
matrice_distance_Jules_filtered <- t(matrice_distance_Jules_filtered)
matrice_distance_Jules_filtered  <- as.data.frame(matrice_distance_Jules_filtered )
colnames(matrice_distance_Jules_filtered) <- matrice_distance_Jules_filtered[1, ]
matrice_distance_Jules_filtered <- matrice_distance_Jules_filtered[-1, ]
dist_or_jaccard <- matrice_distance_Jules_filtered 
save(dist_or_jaccard ,  file='dist_or_jaccard.RData')

setwd("~/S9_ACO/Projet 2 mois/M2PI-Mendelienne/Data")
embedding_jean <- read.csv("embeddings_pheno_SBERT.csv")
coord_fact_embedding <- embedding_jean
setwd("~/S9_ACO/Projet 2 mois/M2PI-Mendelienne/Workflow final")
save(coord_fact_embedding, file = "coord_fact_embedding.RData")


matrice_similarite <- read.csv("agro_to.csv")
matrice_similarite_t <- matrice_similarite %>%
  select(-c(Acute.bronchospasm, 
            Nonspecific.abnormal.findings.on.radiological.and.other.examination.of.other.intrathoracic.organs..echocardiogram..etc.))
matrice_similarite_t <- t(matrice_similarite_t)
matrice_similarite_t  <- as.data.frame(matrice_similarite_t)
colnames(matrice_similarite_t ) <- matrice_similarite_t [1, ]
matrice_similarite_t  <- matrice_similarite_t[-1, ]
matrice_similarite_t_numeric <- matrice_similarite_t %>%
  mutate(across(everything(), ~ as.numeric(as.character(.))))
matrice_distance_apres_OT <- as.data.frame(1 - matrice_similarite_t_numeric)
dist_ot_jaccard <- matrice_distance_apres_OT
save(dist_ot_jaccard, file = "dist_ot_jaccard.RData")

setwd("~/S9_ACO/Projet 2 mois/M2PI-Mendelienne/Data")
Profils_OMIM <- read.csv("Profils_OMIM.csv")
Profils_Phecodes <- read.csv("Profils_Phecodes.csv")

setwd("~/S9_ACO/Projet 2 mois/M2PI-Mendelienne/Workflow final")
save(Profils_OMIM, file = "Profils_OMIM.RData")
save(Profils_Phecodes, file = "Profils_Phecodes.RData")

#################################################################

# Fonction pour créer la matrice d'association
create_association_matrix <- function(distance_matrix, threshold_table) {
  # Conversion en dataframe si nécessaire
  dist_df <- as.data.frame(distance_matrix)
  # Créer une matrice de même dimension que dist_df
  association_matrix <- matrix(0, nrow = nrow(dist_df), ncol = ncol(dist_df))
  # Pour chaque ligne
  for(i in 1:nrow(dist_df)) {
    # Récupérer le seuil correspondant à cette ligne
    current_threshold <- threshold_table$threshold[i]

    # Comparer chaque valeur de distance au seuil
    association_matrix[i,] <- ifelse(dist_df[i,] <= current_threshold, 1, 0)
  }
  # Convertir en dataframe et copier les noms de lignes/colonnes
  association_matrix <- as.data.frame(association_matrix)
  rownames(association_matrix) <- rownames(dist_df)
  colnames(association_matrix) <- colnames(dist_df)
  return(association_matrix)
}

# Utilisation de la fonction
association_df <- create_association_matrix(dist_or_sorensen_mx, tableau_permutation)
rowSums(association_df)

# Vérification
print(dim(association_df))  # Devrait être identique à dim(dist_or_sorensen)
head(association_df[,1:10]) # Afficher les 10 premières colonnes des premières lignes

rm(result_matrix_sorensen)
rm(dist_or_jaccard)
rm(dist_or_sorensen)
rm(dist_ot_jaccard)
################################################################
# Load required packages
library(vegan)  # For Procrustes analysis
library(ggplot2)
library(tidyverse)

# Préparer les configurations avec toutes les dimensions (384)
configs <- array(0, dim = c(nrow(maladies_acm), 384, 6))
configs[,,1] <- scale(as.matrix(maladies_acm))
configs[,,2] <- scale(as.matrix(maladies_sorensen))
configs[,,3] <- scale(as.matrix(maladies_ochiai))
configs[,,4] <- scale(as.matrix(maladies_jaccard))
configs[,,5] <- scale(as.matrix(maladies_ot_jaccard))
configs[,,6] <- scale(as.matrix(maladies_embedding))
configs <- as.data.frame(configs)
# GPA avec toutes les dimensions
gpa_full <- FactoMineR::GPA(configs, 
                            scale = TRUE, 
                            group = c(384, 384, 384, 384, 384, 384),
                            name.group = c("ACM", "Sorensen", "Ochiai", "Jaccard", "Jaccard_OT", "Embedding"))
# Visualiser les résultats
plot(gpa_full, axes = c(1,2))  # On peut toujours visualiser en 2D

# Sauvegarder les résultats
save(gpa_full, file = "gpa_full_384dim.RData")


############################################################################

# library(conclust)
# # Créer les contraintes cannot-link entre maladies complexes
# create_cannot_link_constraints <- function(mc_indices) {
#   cannot_link <- matrix(nrow = 0, ncol = 2)
#   for(i in 1:(length(mc_indices)-1)) {
#     for(j in (i+1):length(mc_indices)) {
#       cannot_link <- rbind(cannot_link, c(mc_indices[i], mc_indices[j]))
#     }
#   }
#   return(cannot_link)
# }

# Application
# cannot_link <- create_cannot_link_constraints(mc_indices)
# result_classif_sorensen <- cop_kmeans(data = data_classif_sorensen,
#                      k = length(mc_indices),
#                      must_link = NULL,
#                      cannot_link = cannot_link)
# result_classif_ochiai <- cop_kmeans(data = data_classif_ochiai,
#                                     k = length(mc_indices),
#                                     mc_indices = mc_indices)
# result_classif_jaccard <- cop_kmeans(data = data_classif_jaccard,
#                                      k = length(mc_indices),
#                                      mc_indices = mc_indices)
# result_classif_ot_jaccard <- cop_kmeans(data = data_classif_ot_jaccard,
#                                         k = length(mc_indices),
#                                         mc_indices = mc_indices)
# result_classif_embedding <- cop_kmeans(data = data_classif_embedding,
#                                        k = length(mc_indices),
#                                        mc_indices = mc_indices)
# result_classif_acm <- cop_kmeans(data = data_classif_acm,
#                                  k = length(mc_indices),
#                                  mc_indices = mc_indices)