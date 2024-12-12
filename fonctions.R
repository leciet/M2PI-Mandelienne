# generate_permutations <- function(row, n_perm = 200) {
#   permutation_list <- vector("list", n_perm)
#   # Récupérer le rowname de la ligne d'entrée
#   original_rowname <- rownames(row)
#   for(i in 1:n_perm) {
#     # Permuter la ligne
#     permuted_row <- sample(row, length(row), replace = FALSE)
#     # Créer un nouveau dataframe avec la ligne permutée
#     new_df <- data.frame(matrix(permuted_row, nrow=1))
#     # Assigner le rowname original
#     rownames(new_df) <- original_rowname
#     # Stocker dans la liste
#     permutation_list[[i]] <- new_df
#   }
#   return(permutation_list)
# }
# permutations <- generate_permutations(or_df0[6127,])
# View(permutations[[2]])

rm(omim)
rm(phecode)

generate_permutations <- function(row, n_perm = 200) {
  permutation_list <- vector("list", n_perm)
  # Récupérer le rowname de la ligne d'entrée
  original_rowname <- rownames(row)
  for(i in 1:n_perm) {
    # Permuter la ligne
    permuted_row <- sample(row, length(row), replace = FALSE)
    # Créer un nouveau dataframe avec la ligne permutée
    new_df <- data.frame(matrix(permuted_row, nrow=1))
    # Assigner le rowname original
    rownames(new_df) <- original_rowname
    # Stocker dans la liste
    permutation_list[[i]] <- new_df
  }
  return(permutation_list)
}

debut <- proc.time()
permutations <- generate_permutations(or_df0[6127,])
fin <- proc.time()
temps_execution <- fin - debut
print(temps_execution)  # 341 secondes soit environ 5,6 minutes
View(permutations)

dataframe <- as.data.frame(permutations[[5]])
View(dataframe)
colonnes_indices <- which(dataframe[1, ] == 1)
colonnes_noms <- names(dataframe)[which(dataframe[1, ] == 1)]

############################################# validé jusque là


# Initialiser un compteur global
counter <- 0
compute_euclidean_distances <- function(permuted_df, or_df0) {
  # Calcule la distance euclidienne entre une ligne permutée et toutes les maladies simples
  distances <- apply(or_df0[1:6125,], 1, function(simple_disease) 
    {
    # Extraction des valeurs des listes
    v1 <- as.numeric(simple_disease)
    v2 <- unlist(permuted_df[1,]) # délistifie la ligne permutée
    v2 <- as.data.frame(v2)
    # Incrémenter le compteur et afficher la progression
    counter <<- counter + 1
    cat(sprintf("\rCalcul de la distance %d/6125", counter))
    flush.console()
    # Calcul de la distance
    sqrt(sum((v1 - v2)^2))
  })
  return(distances)
}

View(permutations[[1]])
View(unlist(permutations[[1]]))
liste <- permutations[[1]]
yo <- unlist(permutations[[1]])
yo2 <- as.data.frame(yo)

# Utilisation
counter <- 0  # Réinitialiser le compteur
debut <- proc.time()
deuxieme_fonction <- compute_euclidean_distances(permutations[[1]], or_df0)
fin <- proc.time()
temps_execution <- fin - debut
print(temps_execution)  # 478 sec donc 8 minutes

troisieme_essai_deuxieme_fonction <- compute_euclidean_distances(permutations[[1]], or_df0)
View(troisieme_essai_deuxieme_fonction)

############## validé mais bcp trop long

compute_euclidean_distances_v2 <- function(permuted_df, or_df0) {
  # Convertir la ligne permutée une seule fois
  v2 <- unlist(permuted_df[1,])
  #v2 <- as.data.frame(v2)
  # Convertir toutes les maladies simples en une seule matrice
  mat_maladies <- as.matrix(or_df0[1:6125,])
  # Calcul vectorisé des distances
  distances <- sqrt(rowSums((mat_maladies - v2)^2))
  # Afficher une seule fois à la fin
  cat("Calcul vectorisé terminé - 6125 distances calculées\n")
  return(distances)
}
# Test version vectorisée
debut2 <- proc.time()
dist2 <- compute_euclidean_distances_v2(permutations[[1]], or_df0)
fin2 <- proc.time()
temps2 <- fin2 - debut2
print(temps2)

# 8,88 secondes c'est imbattable 

# c'est validé

##########################################################
compute_all_distances_v2 <- function(permutation_list, or_df0) {
  lapply(permutation_list, function(perm) compute_euclidean_distances_v2(perm, or_df0))
}

# Exécution avec mesure du temps
debut <- proc.time()
resultats <- compute_all_distances_v2(permutations, or_df0)
fin <- proc.time()
temps_total <- fin - debut
print(temps_total)   # 339,42 secondes soit 6 minutes environ
View(resultats)

str(resultats)
resultats$

resultats2 <- unlist(resultats)
View(resultats2)
resultats2 <- as.data.frame(resultats2)

dataframe_final <- matrix(0, nrow = 200, ncol = 6125)
dataframe_final <- as.data.frame(dataframe_final)
View(dataframe_final)

for (j in 1:200) {
  for (i in 1:6125){
  l=1
  dataframe_final[j, i] <- resultats2[l, 1]
  l=l+1
  }
}
View(dataframe_final)

min_distances <- apply(dataframe_final, 1, min)
View(min_distances)
min_distances <- as.data.frame(min_distances)
min_of_mins <- min(min_distances)
max_of_mins <- max(min_distances)

View(permutations[[1]])
View(permutations[[1]][,2])
rownames(permutations[[1]])
nom_mc <- rownames(permutations[[1]])

which(rownames(or_df0)==nom_mc)
or_df0[which(rownames(or_df0)==nom_mc),]
View(or_df0[which(rownames(or_df0)==nom_mc),])
mc_originale <- or_df0[which(rownames(or_df0)==nom_mc),]

View(dist2)
View(mc_originale)
mat2_maladies <- as.matrix(or_df0[1:6125,])
mc_originale <- as.integer(mc_originale)
dim(mc_originale)
dim(mat2_maladies)
# a <- unlist(permutations[[1]][1,])
# View(a)
distances_mc_originale <- sqrt(rowSums((mat2_maladies - mc_originale)^2))
View(distances_mc_originale)

dataframe_graphique <- rbind(dataframe_final, distances_mc_originale)
View(dataframe_graphique)
str(dataframe_graphique)
summary(dataframe_graphique)
dim(dataframe_graphique)

dataframe_graphique2 <- as.matrix(dataframe_graphique)
transposee_dataframe_graphique <- t(dataframe_graphique2)
plot(dataframe_graphique)

library(DataExplorer)
plot_histogram(dataframe_graphique)
library(ggmcmc)
ggs_histogram(dataframe_graphique, family = "V", bins = 30, greek = FALSE)



# generate_permutations_and_distances(or_df0[6127,], or_df0)
# 
# # 1. Récupération des distances minimales
# min_distances <- sapply(resultats$distances, min)
# # 2. Préparation des données pour ggplot
# df_distances <- data.frame(
#   min_distance = min_distances,
#   permutation = 1:length(min_distances)
# )
# # 3. Calcul des valeurs pour les lignes verticales
# min_of_mins <- min(min_distances)
# max_of_mins <- max(min_distances)
# 
# # 4. Création du graphique
# library(ggplot2)
# ggplot(df_distances, aes(x = min_distance)) +
#   geom_histogram(binwidth = (max(min_distances) - min(min_distances))/30, 
#                  fill = "skyblue", color = "black", alpha = 0.7) +
#   geom_vline(xintercept = min_of_mins, color = "red", linetype = "dashed", size = 1) +
#   geom_vline(xintercept = max_of_mins, color = "blue", linetype = "dashed", size = 1) +
#   annotate("text", x = min_of_mins, y = Inf, label = sprintf("Min: %.3f", min_of_mins),
#            vjust = 2, hjust = -0.1, color = "red") +
#   annotate("text", x = max_of_mins, y = Inf, label = sprintf("Max: %.3f", max_of_mins),
#            vjust = 2, hjust = 1.1, color = "blue") +
#   labs(title = "Distribution des distances minimales pour les 200 permutations",
#        x = "Distance minimale",
#        y = "Fréquence") +
#   theme_minimal() +
#   theme(
#     plot.title = element_text(hjust = 0.5),
#     axis.text = element_text(size = 10),
#     axis.title = element_text(size = 12)
#   )
