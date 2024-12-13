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
permutations <- generate_permutations(or_df0[6129,])
fin <- proc.time()
temps_execution <- fin - debut
print(temps_execution)  # 341 secondes soit environ 5,6 minutes
View(permutations)

dataframe <- as.data.frame(permutations[[5]])
View(dataframe)
colonnes_indices <- which(dataframe[1, ] == 1)
colonnes_noms <- names(dataframe)[which(dataframe[1, ] == 1)]

dataframe_global <- do.call(rbind, permutations)
View(dataframe_global)
# Compare toutes les lignes avec la première
all(sapply(2:nrow(dataframe_global), function(i) 
  isTRUE(all.equal(dataframe_global[1,], dataframe_global[i,]))))

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
resultats

resultats2 <- unlist(resultats)
View(resultats2)
resultats2 <- as.data.frame(resultats2)

dataframe_final <- matrix(0, nrow = 200, ncol = 6125)
dataframe_final <- as.data.frame(dataframe_final)
View(dataframe_final)

# Convertir la liste de résultats en matrice directement
dataframe_final <- matrix(unlist(resultats), nrow = 200, byrow = TRUE)
dataframe_final <- as.data.frame(dataframe_final)
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

library(DataExplorer)
plot_histogram(dataframe_graphique)
library(ggmcmc)
ggs_histogram(dataframe_graphique, family = "V", bins = 30, greek = FALSE)


##########################################################  Sorensen

compute_sorensen_distances_v2 <- function(permuted_df, or_df0) {
  # Convertir la ligne permutée une seule fois
  v2 <- unlist(permuted_df[1,])
  
  # Convertir toutes les maladies simples en une seule matrice
  mat_maladies <- as.matrix(or_df0[1:6125,])
  
  # Calcul vectorisé des distances de Sorensen
  # Pour chaque ligne, calculer :
  # 1 - (2 * nombre d'éléments communs) / (somme des éléments dans chaque vecteur)
  common <- rowSums(mat_maladies & v2)  # éléments communs (1 et 1)
  total <- rowSums(mat_maladies) + sum(v2)  # somme des éléments dans chaque vecteur
  
  distances <- 1 - (2 * common / total)
  
  cat("Calcul vectorisé terminé - 6125 distances de Sorensen calculées\n")
  return(distances)
}

# Test
debut <- proc.time()
dist_sorensen <- compute_sorensen_distances_v2(permutations[[1]], or_df0)
fin <- proc.time()   #11 seconde c'est sympa aussi
temps <- fin - debut
print(temps)

compute_all_distances_sorensen <- function(permutation_list, or_df0) {
  lapply(permutation_list, function(perm) compute_sorensen_distances_v2(perm, or_df0))
}

# Exécution avec mesure du temps
debut <- proc.time()
resultats_sorensen <- compute_all_distances_sorensen(permutations, or_df0)
fin <- proc.time()  #316 secondes
temps_total <- fin - debut
print(temps_total)

# Manipulation des résultats
resultats2_sorensen <- unlist(resultats_sorensen)
resultats2_sorensen <- as.data.frame(resultats2_sorensen)

# Convertir la liste de résultats en matrice directement
dataframe_final_sorensen <- matrix(unlist(resultats_sorensen), nrow = 200, byrow = TRUE)
dataframe_final_sorensen <- as.data.frame(dataframe_final_sorensen)

# Calculer les minimums
min_distances_sorensen <- apply(dataframe_final_sorensen, 1, min)
min_distances_sorensen <- as.data.frame(min_distances_sorensen)
min_of_mins_sorensen <- min(min_distances_sorensen)
max_of_mins_sorensen <- max(min_distances_sorensen)

# Récupérer la maladie complexe originale
nom_mc <- rownames(permutations[[1]])
nom_mc
mc_originale <- or_df0[which(rownames(or_df0)==nom_mc),]
mat2_maladies <- as.matrix(or_df0[1:6125,])
mc_originale <- as.integer(mc_originale)

# Calculer les distances de Sorensen pour la maladie complexe originale
common_original <- rowSums(mat2_maladies & mc_originale)
total_original <- rowSums(mat2_maladies) + sum(mc_originale)
distances_mc_originale_sorensen <- 1 - (2 * common_original / total_original)

# Créer le dataframe final pour le graphique
dataframe_graphique_sorensen <- rbind(dataframe_final_sorensen, distances_mc_originale_sorensen)
View(dataframe_graphique_sorensen)

# Préparation des données pour les courbes
df_long <- as.data.frame(dataframe_graphique_sorensen) %>%
  mutate(row_id = row_number()) %>%
  gather(key = "column", value = "value", -row_id)

# 2. Plot des densités
p2 <- ggplot(df_long, aes(x = value, group = row_id)) +
  geom_density(aes(color = row_id == 201), alpha = 0.1) +
  scale_color_manual(values = c("black", "red")) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "Densités des distances pour chaque permutation",
       x = "Distance",
       y = "Densité")
p2

p2 <- ggplot(df_long, aes(x = value, group = row_id)) +
  geom_density(aes(color = row_id == 201), alpha = 0.1) +
  scale_color_manual(values = c("black", "red")) +
  # Limiter les axes x et y
  coord_cartesian(xlim = c(0.85, 1), ylim = c(0, 50)) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "Densités des distances pour chaque permutation",
       x = "Distance",
       y = "Densité")
p2

