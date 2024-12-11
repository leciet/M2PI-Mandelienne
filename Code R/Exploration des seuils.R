# Nemo Didier
# /12/2024
# Seuil pour les matrices d'assignation post acm 
# ------------------------------------------------------------------------------

# Espace de travail
rm(list = ls())

# Importation des données s
fichiers_rdata <- list.files(path = 'Data',pattern = "\\.RData$")
for (fichier in fichiers_rdata) {
  load(paste('Data/',fichier,sep = ''))
}

#####################################################################################Distance euclidienne

acm_eucli_df <- as.matrix(acm_eucli_df)
View(acm_eucli_df)
nom_col <- colnames(acm_eucli_df)
nom_col
############################################################Approche globale sur tout le tableau

# idée, méthode statistique Otsu, apprentissage supervisée SVM

otsu_threshold <- function(distances) {
  # Aplatir la matrice
  flat_distances <- as.vector(distances)
  # Créer l'histogramme
  hist_data <- hist(flat_distances, breaks = 256, plot = FALSE)
  bin_centers <- (hist_data$breaks[-1] + hist_data$breaks[-length(hist_data$breaks)]) / 2
  hist_counts <- hist_data$counts / sum(hist_data$counts)
  # Calculer les sommes cumulatives
  weight1 <- cumsum(hist_counts)
  weight2 <- rev(cumsum(rev(hist_counts)))
  # Calculer les moyennes
  mean1 <- cumsum(hist_counts * bin_centers) / weight1
  mean2 <- rev(cumsum(rev(hist_counts * bin_centers)) / weight2)
  # Calculer la variance inter-classe
  variance <- weight1[-length(weight1)] * weight2[-1] * 
    (mean1[-length(mean1)] - mean2[-1])^2
  # Trouver le seuil optimal
  optimal_idx <- which.max(variance)
  threshold <- bin_centers[optimal_idx]
  return(threshold)
}

flat_distances <- as.vector(acm_eucli_df)
hist_data <- hist(flat_distances, breaks = 256, plot = FALSE)
hist_data 
plot(hist_data)

seuil <- otsu_threshold(acm_eucli_df)
print(seuil)

# maladies_complexes_apres_selection <- acm_eucli_df <= seuil
# 
# matrice_TRUE <- rowSums(maladies_complexes_apres_selection)
# View(matrice_TRUE)
# 
# any(apply(maladies_complexes_apres_selection, 1, all) == FALSE)
# sum(apply(maladies_complexes_apres_selection, 1, function(row) all(row == FALSE)))
# 
# # Supprimer les lignes avec que des FALSE
# matrice_filtered <- maladies_complexes_apres_selection[!apply(maladies_complexes_apres_selection, 1, function(row) all(row == FALSE)),]
# matrice_filtered
# 
# matrice_filtered[matrice_filtered==TRUE] <- 1
# matrice_filtered[matrice_filtered==FALSE] <- 0

# Créer la heatmap
# library(pheatmap)
# matrice_filtered
# View(matrice_filtered)
# pheatmap(matrice_filtered)
# 
# # Si vous voulez personnaliser la heatmap :
# pheatmap(matrice_filtered,
#          cluster_rows = TRUE,
#          cluster_cols = TRUE,
#          show_rownames = TRUE,
#          show_colnames = TRUE,
#          color = colorRampPalette(c("white", "red"))(50))

############################################################# Approche individuelle sur toutes les lignes
# 
# # Pour chaque ligne, calculer le premier quantile et l'appliquer comme seuil
# seuils_lignes <- apply(acm_eucli_df, 1, function(row) {
#   quantiles <- quantile(row, probs = seq(0, 1, 0.1))
#   return(quantiles[2])  # Premier décile (10%)
# })
# maladies_complexes_apres_selection_individuelle <- t(apply(acm_eucli_df, 1,
#                                               function(row, seuil) row <= seuil,
#                                               seuil = seuils_lignes))
# 
# # Même traitement qu'avant pour la visualisation
# 
# any(apply(maladies_complexes_apres_selection_individuelle, 1, function(row) all(row == FALSE)))
# sum(apply(maladies_complexes_apres_selection_individuelle, 1, function(row) all(row == FALSE)))
# matrice_filtered_individuelle <- maladies_complexes_apres_selection_individuelle[!apply(maladies_complexes_apres_selection_individuelle,
#                                                                            1,
#                                                               function(row) all(row == FALSE)),]
# View(acm_eucli_df)
# View(matrice_filtered_individuelle)
# 
# matrice_TRUE_individuelle <- rowSums(matrice_filtered_individuelle)
# View(matrice_TRUE_individuelle)
# 
# lignes_filtrees <- maladies_complexes_apres_selection_individuelle[matrice_TRUE_individuelle < 10,]
# View(lignes_filtrees)
# 
# row_index <- which(rownames(matrice_filtered_individuelle) == "Cardiac congenital anomalies")
# colnames(matrice_filtered_individuelle)[matrice_filtered_individuelle[row_index,] == TRUE]
# 
# row_index <- which(rownames(matrice_filtered_individuelle) == "Diplopia and disorders of binocular vision")
# colnames(matrice_filtered_individuelle)[matrice_filtered_individuelle[row_index,] == TRUE]
# 
# matrice_filtered_individuelle[matrice_filtered_individuelle == TRUE] <- 1
# matrice_filtered_individuelle[matrice_filtered_individuelle == FALSE] <- 0
# 
# pheatmap(matrice_filtered_individuelle,
#          cluster_rows = TRUE,
#          cluster_cols = TRUE,
#          show_rownames = TRUE,
#          show_colnames = TRUE,
#          color = colorRampPalette(c("white", "red"))(50))
# 
# ################ Autres seuils individuels
# 
# seuils_lignes <- apply(acm_eucli_df, 1, function(row) {
#   quantiles <- quantile(row, probs = c(0.00001))
#   return(quantiles)
# })

############################################################################ Standardisation euclidienne
library(reshape2)
library(dplyr)
library(stringr)

# Transformation au format long
View(acm_eucli_df)
acm_eucli_df <- as.matrix(acm_eucli_df) #matrice des distances
df_eucli <- melt(acm_eucli_df)
View(df_eucli)
colnames(df_eucli ) <- c("mc", "ms", "distance")
df_eucli$ms <- sapply(strsplit(as.character(df_eucli $ms), " "), `[`, 1)
View(df_eucli)

# Obtention du tableau pour le seuil global
s_eucli <- 0.02 # 0.09 obtenu avec le code précédent avec la méthode d'Otsu  
df_filtre_eucli <- df_eucli %>% 
  group_by(mc) %>% 
  filter(distance<=s_eucli) %>%   
  summarise( mc = mc, liste = paste(ms,collapse = " ; "))%>%
  rename(ms = liste)

df_filtre_eucli <- unique(df_filtre_eucli)

comptage <- df_filtre_eucli %>%
  mutate(nb_elements = str_count(ms, ";") + 1,
         ms = gsub("\\s*\\([^)]*\\)", "", ms)) 

View(comptage)
View(df_filtre_eucli)


# Obtention du tableau pour le seuil individuel par ligne

q <- 0.01 # fixe le quantile qui servira de limite 

df_filtre_eucli_ligne <- df_eucli %>% 
  group_by(mc) %>%
  filter(distance <= quantile(distance,q)) %>%   
  summarise( mc = mc, liste = paste(ms,collapse = " ; "))  %>%
  rename(ms = liste)

df_filtre_eucli_ligne  <- unique(df_filtre_eucli_ligne )

comptage_individuel <- df_filtre_eucli_ligne %>%
  mutate(nb_elements = str_count(ms, ";") + 1,
         ms = gsub("\\s*\\([^)]*\\)", "", ms)) 
 
#####################################################################################Distance manhattan

acm_manhattan_df <- as.matrix(acm_manhattan_df)
############################################################ Approche globale
otsu_threshold <- function(distances) {
  flat_distances <- as.vector(distances)
  hist_data <- hist(flat_distances, breaks = 256, plot = FALSE)
  bin_centers <- (hist_data$breaks[-1] + hist_data$breaks[-length(hist_data$breaks)]) / 2
  hist_counts <- hist_data$counts / sum(hist_data$counts)
  weight1 <- cumsum(hist_counts)
  weight2 <- rev(cumsum(rev(hist_counts)))
  mean1 <- cumsum(hist_counts * bin_centers) / weight1
  mean2 <- rev(cumsum(rev(hist_counts * bin_centers)) / weight2)
  variance <- weight1[-length(weight1)] * weight2[-1] * 
    (mean1[-length(mean1)] - mean2[-1])^2
  optimal_idx <- which.max(variance)
  threshold <- bin_centers[optimal_idx]
  return(threshold)
}

flat_distances <- as.vector(acm_manhattan_df)
hist_data <- hist(flat_distances, breaks = 256, plot = FALSE)
hist_data 
plot(hist_data)
seuil_manhattan <- otsu_threshold(acm_manhattan_df)
print(seuil_manhattan)

############################################################################ Standardisation Manhattan
# Transformation au format long
df_manhattan <- melt(acm_manhattan_df)
colnames(df_manhattan) <- c("mc", "ms", "distance")
df_manhattan$ms <- sapply(strsplit(as.character(df_manhattan$ms), " "), `[`, 1)

# Obtention du tableau pour le seuil global
# seuil_manhattan = 0,17 définit précédemment
s_manhattan <- 0.04
df_filtre_manhattan <- df_manhattan %>% 
  group_by(mc) %>% 
  filter(distance<=s_manhattan) %>%   
  summarise(mc = mc, liste = paste(ms,collapse = " ; ")) %>%
  rename(ms = liste)
df_filtre_manhattan <- unique(df_filtre_manhattan)
comptage_manhattan <- df_filtre_manhattan %>%
  mutate(nb_elements = str_count(ms, ";") + 1) 
View(comptage_manhattan)
max(comptage_manhattan$nb_elements)

# Obtention du tableau pour le seuil individuel par ligne
q <- 0.01
df_filtre_manhattan_ligne <- df_manhattan %>% 
  group_by(mc) %>%
  filter(distance <= quantile(distance,q)) %>%   
  summarise(mc = mc, liste = paste(ms,collapse = " ; ")) %>%
  rename(ms = liste)
df_filtre_manhattan_ligne <- unique(df_filtre_manhattan_ligne)
comptage_individuel_manhattan <- df_filtre_manhattan_ligne %>%
  mutate(nb_elements = str_count(ms, ";") + 1)
is.data.frame(comptage_individuel_manhattan)

############################################################################ Heatmap ou graphique d'association

View(comptage)
View(comptage_individuel)
View(df_eucli)
View(df_filtre_eucli)
View(df_filtre_eucli_ligne)

View(comptage_manhattan)
View(comptage_individuel_manhattan)
View(df_manhattan)
View(df_filtre_manhattan)
View(df_filtre_manhattan_ligne)

# Manhattan individuel

df_filtre_manhattan_ligne$ms <- gsub("\\s*;\\s*", ";", df_filtre_manhattan_ligne$ms)
head(df_filtre_manhattan_ligne$ms)

library(tidyr)
data_long <- df_filtre_manhattan_ligne %>%
  separate_rows(ms, sep = ";")
View(data_long)
data_long
data_long <- unique(data_long)
association_matrix <- table(data_long$ms, data_long$mc)
association_matrix
association_matrix <- as.data.frame(association_matrix)
View(association_matrix)

unique(association_matrix$Freq)
lines_with_freq_2 <- association_matrix[association_matrix$Freq == 2, ]
View(lines_with_freq_2)

lignes_dupliquees <- df_filtre_manhattan_ligne[duplicated(df_filtre_manhattan_ligne) | duplicated(df_filtre_manhattan_ligne, fromLast = TRUE), ]
View(lignes_dupliquees)
lignes_dupliquees <- data_long[duplicated(data_long) | duplicated(data_long, fromLast = TRUE), ]
View(lignes_dupliquees)

 library(tidyr)
 association_wide <- association_matrix %>%
   pivot_wider(names_from = Var2, values_from = Freq)
 association_wide <- as.data.frame(association_wide)
 rownames(association_wide) <- association_wide$Var1
 association_wide <- as.matrix(association_wide[, -1]) # Enlever la colonne Var1
 View(association_wide)
 
table(association_wide)
 
 heatmap(
   association_wide,
   col = c("white", "black"), # Blanc pour 0, noir pour 1
   scale = "none",            # Pas de normalisation
   Rowv = NA,                 # Pas de clustering des lignes
   Colv = NA,                 # Pas de clustering des colonnes
   main = "Binary Heatmap (Native)"
 )
 
 # Pour Manhattan global
 df_filtre_manhattan$ms <- gsub("\\s*;\\s*", ";", df_filtre_manhattan$ms)
 
 data_long_manhattan <- df_filtre_manhattan %>%
   separate_rows(ms, sep = ";")
 data_long_manhattan <- unique(data_long_manhattan)
 
 association_matrix_manhattan <- table(data_long_manhattan$ms, data_long_manhattan$mc)
 association_matrix_manhattan <- as.data.frame(association_matrix_manhattan)
 
 association_wide_manhattan <- association_matrix_manhattan %>%
   pivot_wider(names_from = Var2, values_from = Freq)
 
 association_wide_manhattan <- as.data.frame(association_wide_manhattan)
 rownames(association_wide_manhattan) <- association_wide_manhattan$Var1
 association_wide_manhattan <- as.matrix(association_wide_manhattan[, -1])
 
 heatmap(
   association_wide_manhattan,
   col = c("white", "black"),
   scale = "none",
   Rowv = NA,
   Colv = NA,
   main = "Manhattan Global"
 )
 
 # Pour Euclidien global
 df_filtre_eucli$ms <- gsub("\\s*;\\s*", ";", df_filtre_eucli$ms)
 
 data_long_eucli <- df_filtre_eucli %>%
   separate_rows(ms, sep = ";") 
 data_long_eucli <- unique(data_long_eucli)
 
 association_matrix_eucli <- table(data_long_eucli$ms, data_long_eucli$mc)
 association_matrix_eucli <- as.data.frame(association_matrix_eucli)
 
 association_wide_eucli <- association_matrix_eucli %>%
   pivot_wider(names_from = Var2, values_from = Freq)
 
 association_wide_eucli <- as.data.frame(association_wide_eucli)
 rownames(association_wide_eucli) <- association_wide_eucli$Var1
 association_wide_eucli <- as.matrix(association_wide_eucli[, -1])
 
 heatmap(
   association_wide_eucli,
   col = c("white", "black"),
   scale = "none",
   Rowv = NA,
   Colv = NA,
   main = "Euclidien Global"
 )
 
 # Pour Euclidien individuel
 df_filtre_eucli_ligne$ms <- gsub("\\s*;\\s*", ";", df_filtre_eucli_ligne$ms)
 
 data_long_eucli_ligne <- df_filtre_eucli_ligne %>%
   separate_rows(ms, sep = ";") 
 data_long_eucli_ligne <- unique(data_long_eucli_ligne)
 
 association_matrix_eucli_ligne <- table(data_long_eucli_ligne$ms, data_long_eucli_ligne$mc)
 association_matrix_eucli_ligne <- as.data.frame(association_matrix_eucli_ligne)
 
 association_wide_eucli_ligne <- association_matrix_eucli_ligne %>%
   pivot_wider(names_from = Var2, values_from = Freq)
 
 association_wide_eucli_ligne <- as.data.frame(association_wide_eucli_ligne)
 rownames(association_wide_eucli_ligne) <- association_wide_eucli_ligne$Var1
 association_wide_eucli_ligne <- as.matrix(association_wide_eucli_ligne[, -1])
 
 heatmap(
   association_wide_eucli_ligne,
   col = c("white", "black"),
   scale = "none",
   Rowv = NA,
   Colv = NA,
   main = "Euclidien Individuel"
 )
 
View(association_wide_eucli)
association_wide
association_wide_eucli
association_wide_eucli_ligne
association_wide_manhattan

# Pour sauvegarder en RData
save(association_wide,
     file = "matrices_asso_manhattan_individuel.RData")

save(association_wide_eucli, 
     file = "matrices_asso_euclidien_global.RData")

save(association_wide_eucli_ligne,
     file = "matrices_asso_euclidien_individuel.RData")

save(association_wide_manhattan,
     file = "matrices_asso_manhattan_global.RData")

# Pour exporter en CSV
write.csv(association_wide, "asso_manhattan_individuel.csv", row.names = TRUE)
write.csv(association_wide_eucli, "asso_euclidien_global.csv", row.names = TRUE)
write.csv(association_wide_eucli_ligne, "asso_euclidien_individuel.csv", row.names = TRUE)
write.csv(association_wide_manhattan, "asso_manhattan_global.csv", row.names = TRUE)

library(ggplot2)

asso_acm_eucli <- data.frame(ifelse(acm_eucli_df<=0.02,1,0))
asso_acm_manhattan <- data.frame(ifelse(acm_manhattan_df<=0.04,1,0))
save(asso_acm_manhattan,
     file = "asso_acm_manhattan.RData")
save(asso_acm_eucli, file = "asso_acm_eucli.RData")

######################################################################################"
# Fonction pour standardiser une matrice de distance
standardize_robust <- function(matrix) {
  med <- median(matrix)
  iqr <- IQR(matrix)
  # Éviter la division par zéro si l'IQR est nul
  if (iqr == 0) {
    stop("L'IQR est nul, la standardisation robuste n'est pas possible.")
  }
  (matrix - med) / iqr
}

# Standardisation robuste pour chaque matrice
dist_acm_eucli_std <- standardize_robust(dist_acm_eucli_mx)
dist_acm_manhattan_std <- standardize_robust(dist_acm_manhattan_mx)
dist_or_ochiai_std <- standardize_robust(dist_or_ochiai_mx)
dist_or_sorensen_std <- standardize_robust(dist_or_sorensen_mx)

# Vérifier les statistiques après transformation
summary(as.vector(dist_acm_eucli_std))
summary(as.vector(dist_acm_manhattan_std))
summary(as.vector(dist_or_ochiai_std))
summary(as.vector(dist_or_sorensen_std))
