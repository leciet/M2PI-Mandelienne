# Présence commune (1-1) plus informative que l'absence commune (0-0)
# sokal_michener peu adapté car donne autant de poids aux absences communes qu'aux présences communes
# Bray-Curtis plus adapté aux données de comptage
# Jaccard donne autant de poids aux présences communes qu'aux absences communes, et ne gère pas
# bien les données creuses
# sokal_sneath est normalement utilisé pour des valeurs continues, pas des valeurs binaires 

# Ochiai permet de gérer les maladies qui ont des différences importantes dans le nombre 
# de phénotypes associés (avec des 1)

# Distance cosinus est adapté pour les très grandes dimensions

# Sorensen et hamann donnent plus de poids aux présences commune

# MDS unfolding ----
# MDS unfolding pour les nouvelles matrices
coord_fact_mds_sokal_sneath <- execute_mds_unfolding(dist_or_sokal_sneath_mx, "sokal_sneath")
coord_fact_mds_sokal_michener <- execute_mds_unfolding(dist_or_sokal_michener_mx, "sokal_michener")
coord_fact_mds_rogers_tanimoto <- execute_mds_unfolding(dist_or_rogers_tanimoto_mx, "rogers_tanimoto")
coord_fact_mds_hamann <- execute_mds_unfolding(dist_or_hamann_mx, "hamann")
coord_fact_mds_cosine <- execute_mds_unfolding(dist_cosine_mx, "cosine")  

# Sauvegarder les résultats MDS
save(coord_fact_mds_sokal_sneath, file = "coord_fact_mds_sokal_sneath.RData")
save(coord_fact_mds_sokal_michener, file = "coord_fact_mds_sokal_michener.RData") 
save(coord_fact_mds_rogers_tanimoto, file = "coord_fact_mds_rogers_tanimoto.RData")
save(coord_fact_mds_hamann, file = "coord_fact_mds_hamann.RData")
save(coord_fact_mds_cosine, file = "coord_fact_mds_cosine.RData")

# Charger les données pour AFM
load("coord_fact_mds_sokal_sneath.RData")
load("coord_fact_mds_sokal_michener.RData")
load("coord_fact_mds_rogers_tanimoto.RData") 
load("coord_fact_mds_hamann.RData")
load("coord_fact_mds_cosine.RData")

# Extraire les coordonnées des maladies complexes pour chaque méthode
maladies_sokal_sneath <- as.data.frame(coord_fact_mds_sokal_sneath$conf.row)
maladies_sokal_michener <- as.data.frame(coord_fact_mds_sokal_michener$conf.row)
maladies_rogers_tanimoto <- as.data.frame(coord_fact_mds_rogers_tanimoto$conf.row)
maladies_hamann <- as.data.frame(coord_fact_mds_hamann$conf.row)
maladies_cosine <- as.data.frame(coord_fact_mds_cosine$conf.row)

# AFM sur maladies complexes ----
# Mettre à jour data_pour_afm en incluant toutes les méthodes
data_pour_afm.toutes.distances <- cbind(maladies_acm, 
                                        maladies_sorensen, 
                                        maladies_ochiai, 
                                        maladies_jaccard, 
                                        maladies_ot_jaccard,
                                        maladies_embedding,  # Données existantes
                       maladies_sokal_sneath,
                       maladies_sokal_michener, 
                       maladies_rogers_tanimoto,
                       maladies_hamann,
                       maladies_cosine)
save(data_pour_afm.toutes.distances, file = "data_pour_afm.toutes.distances.RData")

# Mettre à jour l'AFM avec toutes les méthodes
res.mfa.toutes.distances <- MFA(data_pour_afm.toutes.distances, 
               group = c(384, 384, 384, 384, 384, 384,  # Groupes existants
                         384, 384, 384, 384, 384),  # Nouveaux groupes
               type = c("c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c"), 
               ncp = 2,             
               name.group = c("ACM", "Sorensen", "Ochiai", "Jaccard", "Jaccard_OT", "Embedding",
                              "sokal_sneath", "sokal_michener", "rogers_tanimoto", "Hamann", "Cosine"))

# Sauvegarder l'AFM mise à jour
save(res.mfa.toutes.distances, file = "res.mfa.toutes.distances.RData")

# Création des datas pour les futurs classifs ----

# Classification hiérarchique pour les nouvelles méthodes
data_classif_sokal_sneath <- rbind(coord_fact_mds_sokal_sneath$conf.row, coord_fact_mds_sokal_sneath$conf.col)
data_classif_sokal_michener <- rbind(coord_fact_mds_sokal_michener$conf.row, coord_fact_mds_sokal_michener$conf.col)
data_classif_rogers_tanimoto <- rbind(coord_fact_mds_rogers_tanimoto$conf.row, coord_fact_mds_rogers_tanimoto$conf.col)
data_classif_hamann <- rbind(coord_fact_mds_hamann$conf.row, coord_fact_mds_hamann$conf.col)
data_classif_cosine <- rbind(coord_fact_mds_cosine$conf.row, coord_fact_mds_cosine$conf.col)
data_classif_cosine <- as.data.frame(data_classif_cosine)
data_classif_hamann <- as.data.frame(data_classif_hamann)
data_classif_sokal_michener <- as.data.frame(data_classif_sokal_michener)
data_classif_rogers_tanimoto <- as.data.frame(data_classif_rogers_tanimoto)
# maladies complexes en premier puis simples en second
rm(coord_fact_mds)
# Standardisation des noms de colonnes
names(data_classif_sokal_sneath) <- paste0("D", 1:ncol(data_classif_sokal_sneath))
names(data_classif_sokal_michener) <- paste0("D", 1:ncol(data_classif_sokal_michener))
names(data_classif_rogers_tanimoto) <- paste0("D", 1:ncol(data_classif_rogers_tanimoto))
names(data_classif_hamann) <- paste0("D", 1:ncol(data_classif_hamann))
names(data_classif_cosine) <- paste0("D", 1:ncol(data_classif_cosine))

save(data_classif_sokal_sneath, file="data_classif_sokal_sneath.RData")
save(data_classif_sokal_michener, file="data_classif_sokal_michener.RData")
save(data_classif_rogers_tanimoto, file="data_classif_rogers_tanimoto.RData")
save(data_classif_hamann, file="data_classif_hamann.RData")
save(data_classif_cosine, file="data_classif_cosine.RData")

# Permutations sur toutes les distances ----

# Fonctions vectorisées pour chaque type de distance
calc_distances <- function(matrix1, vector1, method = "sorensen") {
  # Points d'intersection et sommes (utilisés dans plusieurs métriques)
  intersections <- matrix1 %*% vector1
  sums1 <- rowSums(matrix1)
  sum2 <- sum(vector1)
  switch(method,
         "sorensen" = {
           1 - (2 * intersections) / (sums1 + sum2)
         },
         "ochiai" = {
           1 - intersections / sqrt(sums1 * sum2)
         },
         "jaccard" = {
           1 - intersections / (sums1 + sum2 - intersections)
         },
         "sokal_sneath" = {
           intersections / sqrt((sums1 + sum2 - intersections) * intersections)
         },
         "sokal_michener" = {
           rowSums(abs(matrix1 - vector1))
         },
         "rogers_tanimoto" = {
           1 - 0.5 * (intersections/sums1 + intersections/sum2)
         },
         "hamann" = {
             2 * intersections / (sums1 + sum2)
         },
         "cosinus" = {
           1 - intersections / (sqrt(sums1) * sqrt(sum2))
         }
  )
}

# Fonction de permutation pour une ligne
permut_all <- function(dta, n_perm, line_idx, distance_method = "sorensen") {
  dta <- as.matrix(dta)
  row <- dta[line_idx, ]
  original_rowname <- rownames(dta)[line_idx]
  pmin <- numeric(n_perm)
  other_rows <- matrix(as.numeric(dta[-line_idx, ]), 
                       nrow = nrow(dta) - 1,
                       byrow = FALSE)
  for(i in 1:n_perm) {
    permuted_row <- sample(row, length(row), replace = FALSE)
    distances <- calc_distances(other_rows, permuted_row, method = distance_method)
    pmin[i] <- min(distances)
  }
  
  return(c(rowname = original_rowname, threshold = min(pmin)))
}

# Fonction pour traiter plusieurs lignes
process_multiple_rows <- function(dta, n_perm, start_idx, end_idx, distance_method = "sorensen") {
  n_rows <- end_idx - start_idx + 1
  results <- matrix(nrow = n_rows, ncol = 2)
  for(i in 1:n_rows) {
    curr_idx <- start_idx + i - 1
    results[i, ] <- permut_all(dta, n_perm, curr_idx, distance_method)
  }
  results_df <- as.data.frame(results)
  colnames(results_df) <- c("rowname", "threshold")
  return(results_df)
}

# Fonction pour exécuter toutes les méthodes de distance
run_all_distances <- function(phenotype_matrix, n_perm = 10, start_idx = 6103, end_idx = 7064) {
  distance_methods <- c("sorensen", "ochiai", "jaccard", "sokal_sneath", 
                        "sokal_michener", "rogers_tanimoto", "hamann", "cosinus")
  results <- list()
  for(method in distance_methods) {
    cat("Processing", method, "distance...\n")
    results[[method]] <- process_multiple_rows(
      phenotype_matrix, 
      n_perm = n_perm, 
      start_idx = start_idx, 
      end_idx = end_idx,
      distance_method = method
    )
    results[[method]]$threshold <- as.numeric(results[[method]]$threshold)
  }
  return(results)
}
# rm(list=ls())
# Exécuter toutes les distances
results_all <- run_all_distances(phenotype_maladie_s_c)
# Ou exécuter une distance spécifique
# results_sorensen <- process_multiple_rows(
#   phenotype_maladie_s_c, 
#   n_perm = 10, 
#   start_idx = 6103, 
#   end_idx = 7064,
#   distance_method = "sorensen"
# )
# Sauvegarder les résultats
save(results_all, file = "results_all_distances.RData")

# Matrice d'assignation obtenus grâce aux seuils des permutations pour toutes les distances ----

# Fonction pour créer les matrices d'association pour toutes les distances
create_all_association_matrices <- function(results_all) {
  # Liste des matrices de distance à charger
  distance_files <- list(
    sorensen = "dist_or_sorensen_mx.RData",
    ochiai = "dist_or_ochiai_mx.RData",
    jaccard = "dist_or_jaccard_mx.RData",
    sokal_sneath = "dist_or_sokal_sneath_mx.RData",
    sokal_michener = "dist_or_sokal_michener_mx.RData",
    rogers_tanimoto = "dist_or_rogers_tanimoto_mx.RData",
    hamann = "dist_or_hamann_mx.RData",
    cosine = "dist_cosine_mx.RData"
  )
  # Stocker les résultats
  association_matrices <- list()
  # Pour chaque distance
  for(dist_name in names(distance_files)) {
    # Charger la matrice de distance
    load(distance_files[[dist_name]])
    # Obtenir le nom de l'objet chargé
    dist_matrix_name <- gsub("\\.RData$", "", distance_files[[dist_name]])
    dist_matrix <- get(dist_matrix_name)
    # Créer la matrice d'association
    cat("Processing", dist_name, "distance...\n")
    # Conversion en dataframe
    dist_df <- as.data.frame(dist_matrix)
    # Créer une matrice vide
    association_matrix <- matrix(0, nrow = nrow(dist_df), ncol = ncol(dist_df))
    # Pour chaque ligne
    for(i in 1:nrow(dist_df)) {
      # Récupérer le seuil correspondant
      current_threshold <- results_all[[dist_name]]$threshold[i]
      # Comparer les valeurs au seuil
      association_matrix[i,] <- ifelse(dist_df[i,] <= current_threshold, 1, 0)
    }
    # Convertir en dataframe
    association_matrix <- as.data.frame(association_matrix)
    rownames(association_matrix) <- rownames(dist_df)
    colnames(association_matrix) <- colnames(dist_df)
    # Stocker la matrice
    association_matrices[[dist_name]] <- association_matrix
    # Sauvegarder la matrice
    save_name <- paste0("association_matrix_", dist_name, ".RData")
    save(association_matrix, file = save_name)
  }
  return(association_matrices)
}

# Créer toutes les matrices d'association
association_matrices <- create_all_association_matrices(results_all)
View(rowSums(association_matrices$sorensen))
save(association_matrices, file="association_matrices.RData")

# Clustering ----

library(clustertend)
hopkins_stat <- hopkins(coord_fact_acm, n = 7063)
print(hopkins_stat)
VAT(coord_fact_acm)

# Implémentation du COP-KMEANS
cop_kmeans <- function(data, k, mc_inhamanns, max_iter = 100) {
  # Initialisation
  n <- nrow(data)
  centroids <- as.matrix(data[mc_inhamanns[1:k], ])  # Utiliser les MC comme centroids initiaux
  clusters <- rep(-1, n)
  clusters[mc_inhamanns] <- 1:k  # Fixer les MC à leurs clusters
  # Itérations
  for(iter in 1:max_iter) {
    old_clusters <- clusters
    # Pour tous les points sauf les MC
    for(i in setdiff(1:n, mc_inhamanns)) {
      # Calculer distances aux centroids
      distances <- apply(centroids, 1, function(cent) {
        sqrt(sum((data[i,] - cent)^2))
      })
      # Assigner au centroid le plus proche
      clusters[i] <- which.min(distances)
    }
    # Mettre à jour les centroids (sauf pour les clusters contenant uniquement une MC)
    for(j in 1:k) {
      cluster_points <- data[clusters == j & !seq_len(n) %in% mc_inhamanns, ]
      if(nrow(as.matrix(cluster_points)) > 0) {
        centroids[j,] <- colMeans(as.matrix(cluster_points))
      }
    }
    # Vérifier convergence
    if(all(clusters == old_clusters)) break
  }
  return(list(clusters = clusters, centroids = centroids))
}

mc_inhamanns=c(1:962)
  
# Application du COP-KMEANS sur chaque jeu de données
result_classif_hamann <- cop_kmeans(data = data_classif_hamann,
                                      k = length(mc_inhamanns),
                                      mc_inhamanns = mc_inhamanns)

result_classif_ochiai <- cop_kmeans(data = data_classif_ochiai,
                                    k = length(mc_inhamanns),
                                    mc_inhamanns = mc_inhamanns)

result_classif_jaccard <- cop_kmeans(data = data_classif_jaccard,
                                     k = length(mc_inhamanns),
                                     mc_inhamanns = mc_inhamanns)

result_classif_ot_jaccard <- cop_kmeans(data = data_classif_ot_jaccard,
                                        k = length(mc_inhamanns),
                                        mc_inhamanns = mc_inhamanns)

result_classif_embedding <- cop_kmeans(data = data_classif_embedding,
                                       k = length(mc_inhamanns),
                                       mc_inhamanns = mc_inhamanns)

result_classif_acm <- cop_kmeans(data = data_classif_acm,
                                 k = length(mc_inhamanns),
                                 mc_inhamanns = mc_inhamanns)

# Création d'une liste pour stocker tous les résultats
resultats_classification <- list(
  sorensen = result_classif_sorensen,
  ochiai = result_classif_ochiai,
  jaccard = result_classif_jaccard,
  ot_jaccard = result_classif_ot_jaccard,
  embedding = result_classif_embedding,
  acm = result_classif_acm
)

# Sauvegarder les résultats
save(resultats_classification, file = "resultats_classification.RData")

# # Bibliothèques nécessaires
# library(cluster)
# library(dplyr)
# library(Matrix)
# # Bibliothèques nécessaires
# library(stats)
# library(progress)  # Pour la barre de progression
# 
# library(progress)
# library(Matrix)

# library(progress)
# library(Matrix)
# 
# fuzzy_constrained_clustering <- function(
#     data,              # Jeu de données complet
#     mc_inhamanns,         # Inhamanns des maladies complexes
#     ms_inhamanns,         # Inhamanns des maladies simples
#     n_clusters,         # Nombre de clusters
#     m = 2,              # Paramètre de fuzzification (défaut: 2)
#     max_iter = 30,      # Nombre maximum d'itérations
#     epsilon = 0.05      # Seuil de convergence
# ) {
#   # Vérifications préliminaires
#   stopifnot(length(mc_inhamanns) > 0, length(ms_inhamanns) > 0)
#   stopifnot(max(c(mc_inhamanns, ms_inhamanns)) <= nrow(data))
#   # Conversion en matrice
#   data <- as.matrix(data)
#   # Initialisation des centres de clusters
#   set.seed(123)
#   cluster_centers <- data[sample(nrow(data), n_clusters), ]
#   # Matrice d'appartenance initiale
#   membership_matrix <- matrix(0, nrow = nrow(data), ncol = n_clusters)
#   # Fonction de calcul des distances
#   calculate_distances <- function(points, centers) {
#     # Calcul vectoriel des distances
#     distances_matrix <- matrix(0, nrow = nrow(points), ncol = nrow(centers))
#     for (i in 1:nrow(points)) {
#       distances_matrix[i, ] <- colSums((t(centers) - points[i, ])^2)
#     }
#     sqrt(distances_matrix)
#   }
#   # Contraintes pour les maladies complexes (appartenance unique)
#   mc_points <- data[mc_inhamanns, , drop = FALSE]
#   mc_distances <- calculate_distances(mc_points, cluster_centers)
#   mc_best_clusters <- apply(mc_distances, 1, which.min)
#   for (i in 1:length(mc_inhamanns)) {
#     membership_matrix[mc_inhamanns[i], ] <- 0
#     membership_matrix[mc_inhamanns[i], mc_best_clusters[i]] <- 1
#   }
#   
#   # Initialisation floue pour les maladies simples
#   ms_points <- data[ms_inhamanns, , drop = FALSE]
#   # Boucle principale de l'algorithme
#   for (iter in 1:max_iter) {
#     # Sauvegarde de l'ancienne matrice d'appartenance
#     old_membership <- membership_matrix
#     # Mise à jour des centres de clusters
#     membership_power <- membership_matrix^m
#     cluster_centers <- t(apply(membership_power, 2, function(cluster_weight) {
#       colSums(data * cluster_weight) / sum(cluster_weight)
#     }))
#     # Calcul des distances pour les maladies simples
#     ms_distances <- calculate_distances(ms_points, cluster_centers)
#     # Mise à jour des appartenances pour les maladies simples
#     for (i in 1:length(ms_inhamanns)) {
#       distances <- ms_distances[i, ]
#       # Calcul des appartenances floues
#       membership_values <- (1 / sqrt(distances))^(2/(m-1))
#       membership_values <- membership_values / sum(membership_values)
#       membership_matrix[ms_inhamanns[i], ] <- membership_values
#     }
#     # Vérification de la convergence
#     if (max(abs(membership_matrix - old_membership), na.rm = TRUE) < epsilon) {
#       break
#     }
#   }
#   # Retour simplifié
#   list(
#     membership_matrix = membership_matrix,
#     cluster_centers = cluster_centers
#   )
# }
# hc_constr_hamann <- fuzzy_constrained_clustering(data_classif_hamann, 
#                                                  mc_inhamanns = c(1:962), 
#                                                  ms_inhamanns = c(963:7064),
#                                                  n_clusters = 962)
# head(data_classif_hamann)
# View(rowSums(as.data.frame(hc_constr_hamann$membership_matrix)))
# View(hc_constr_hamann$membership_matrix)
# hc_constr_hamann <- as.data.frame(hc_constr_hamann)
# sum(is.na(data_classif_hamann))
# sum(is.nan(data_classif_hamann))
# sum(is.infinite(data_classif_hamann))

# hc_constr_sorensen <- fuzzy_constrained_clustering(
#   data = data_classif_sorensen, 
#   mc_inhamanns = c(1:962), 
#   ms_inhamanns = c(963:7064),
#   n_clusters = 962
# )
# 
# hc_constr_embedding <- fuzzy_constrained_clustering(
#   data = data_classif_embedding, 
#   mc_inhamanns = c(1:962), 
#   ms_inhamanns = c(963:7064),
#   n_clusters = 962
# )
# hc_constr_acm <- fuzzy_constrained_clustering(
#   data = data_classif_acm, 
#   mc_inhamanns = c(1:962), 
#   ms_inhamanns = c(963:7064),
#   n_clusters = 962
# )
# 
# save(hc_constr_sorensen, file = "hc_constr_sorensen.RData")
# save(hc_constr_hamann, file = "hc_constr_hamann.RData")
# save(hc_constr_embedding, file = "hc_constr_embedding,.RData")
# save(hc_constr_acm, file = "hc_constr_acm.RData")
# 
# 
# 
# hc_constr_ochiai <- fuzzy_constrained_clustering(
#   data = data_classif_ochiai, 
#   mc_inhamanns = c(1:962), 
#   ms_inhamanns = c(963:7064),
#   n_clusters = 962
# )
# 
# hc_constr_jaccard <- fuzzy_constrained_clustering(
#   data = data_classif_jaccard, 
#   mc_inhamanns = c(1:962), 
#   ms_inhamanns = c(963:7064),
#   n_clusters = 962
# )
# 
# hc_constr_sokal_sneath <- fuzzy_constrained_clustering(
#   data = data_classif_sokal_sneath, 
#   mc_inhamanns = c(1:962), 
#   ms_inhamanns = c(963:7064),
#   n_clusters = 962
# )
# 
# hc_constr_sokal_michener <- fuzzy_constrained_clustering(
#   data = data_classif_sokal_michener, 
#   mc_inhamanns = c(1:962), 
#   ms_inhamanns = c(963:7064),
#   n_clusters = 962
# )
# 
# hc_constr_rogers_tanimoto <- fuzzy_constrained_clustering(
#   data = data_classif_rogers_tanimoto, 
#   mc_inhamanns = c(1:962), 
#   ms_inhamanns = c(963:7064),
#   n_clusters = 962
# )
# 
# hc_constr_cosine <- fuzzy_constrained_clustering(
#   data = data_classif_cosine, 
#   mc_inhamanns = c(1:962), 
#   ms_inhamanns = c(963:7064),
#   n_clusters = 962
# )
# 
# # Sauvegarde des résultats
# save(hc_constr_sorensen, file = "hc_constr_sorensen.RData")
# save(hc_constr_ochiai, file = "hc_constr_ochiai.RData")
# save(hc_constr_jaccard, file = "hc_constr_jaccard.RData")
# save(hc_constr_sokal_sneath, file = "hc_constr_sokal_sneath.RData")
# save(hc_constr_sokal_michener, file = "hc_constr_sokal_michener.RData")
# save(hc_constr_rogers_tanimoto, file = "hc_constr_rogers_tanimoto.RData")
# save(hc_constr_hamann, file = "hc_constr_hamann.RData")
# save(hc_constr_cosine, file = "hc_constr_cosine.RData")

# Créatin des matrices d'assignation à partir du clustering ----

# Création des matrices d'assignation pour les nouvelles méthodes
assignment_matrix_sokal_sneath <- create_assignment_matrix(hc_sokal_sneath$data.clust)
assignment_matrix_sokal_michener <- create_assignment_matrix(hc_sokal_michener$data.clust)
assignment_matrix_rogers_tanimoto <- create_assignment_matrix(hc_rogers_tanimoto$data.clust)
assignment_matrix_hamann <- create_assignment_matrix(hc_hamann$data.clust)
assignment_matrix_cosine <- create_assignment_matrix(hc_cosine$data.clust)
assignment_matrix_bray <- create_assignment_matrix(hc_bray$data.clust)

# Convertir en dataframe
assignment_matrix_sokal_sneath <- as.data.frame(assignment_matrix_sokal_sneath)
assignment_matrix_sokal_michener <- as.data.frame(assignment_matrix_sokal_michener)
assignment_matrix_rogers_tanimoto <- as.data.frame(assignment_matrix_rogers_tanimoto)
assignment_matrix_hamann <- as.data.frame(assignment_matrix_hamann)
assignment_matrix_cosine <- as.data.frame(assignment_matrix_cosine)
assignment_matrix_bray <- as.data.frame(assignment_matrix_bray)

# Mise à jour de la liste des matrices d'assignation
matrices_assignation_classif$assignment_matrix_classif_sokal_sneath <- assignment_matrix_sokal_sneath
matrices_assignation_classif$assignment_matrix_classif_sokal_michener <- assignment_matrix_sokal_michener
matrices_assignation_classif$assignment_matrix_classif_rogers_tanimoto <- assignment_matrix_rogers_tanimoto
matrices_assignation_classif$assignment_matrix_classif_hamann <- assignment_matrix_hamann
matrices_assignation_classif$assignment_matrix_classif_cosine <- assignment_matrix_cosine
matrices_assignation_classif$assignment_matrix_classif_bray <- assignment_matrix_bray

# Sauvegarder la liste mise à jour
save(matrices_assignation_classif, file = "matrices_assignation.RData")