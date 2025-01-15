#Library ----
library(tidyverse)
library(FactoMineR)
library(Factoshiny)
library(ade4)
library(vegan)
library(reshape)
library(DataExplorer)
library(smacof)

#Import initial ----
rm(list=ls())
load("phenotype_maladie_s_c.RData")

######################################################### 1) Matrices de distance
# Import Jaccard distance ----
# matrice_distance_Jules <- read.csv("distance_matrix_dist_Jules.csv")
# library(dplyr)
# matrice_distance_Jules_filtered <- matrice_distance_Jules %>%
#   select(-c(Acute.bronchospasm, 
#             Nonspecific.abnormal.findings.on.radiological.and.other.examination.of.other.intrathoracic.organs..echocardiogram..etc.))
# matrice_distance_Jules_filtered <- t(matrice_distance_Jules_filtered)
# matrice_distance_Jules_filtered  <- as.data.frame(matrice_distance_Jules_filtered )
# colnames(matrice_distance_Jules_filtered) <- matrice_distance_Jules_filtered[1, ]
# matrice_distance_Jules_filtered <- matrice_distance_Jules_filtered[-1, ]
# dist_or_jaccard <- matrice_distance_Jules_filtered 
# save(dist_or_jaccard ,  file='dist_or_jaccard.RData')

load("dist_or_jaccard.RData")
dist_or_jaccard_mx <- as.matrix(dist_or_jaccard)

# Calculate Sorensen distance ----
dist_or_sorensen <- dist.binary(phenotype_maladie_s_c, method = 5)
dist_or_sorensen <- as.data.frame(as.matrix(dist_or_sorensen))
dist_or_sorensen <- dist_or_sorensen[6126:7089,1:6125]
dist_or_sorensen_mx <- as.matrix(dist_or_sorensen)

# Calculate Ochiai distance ----

dist_or_ochiai <- dist.binary(phenotype_maladie_s_c, method = 7)
dist_or_ochiai <- as.data.frame(as.matrix(dist_or_ochiai))
dist_or_ochiai <- dist_or_ochiai[6126:7089,1:6125]
dist_or_ochiai_mx <- as.matrix(dist_or_ochiai)

# Save distance matrices ----

save(dist_or_sorensen_mx, file="dist_or_sorensen_mx.RData")
save(dist_or_ochiai_mx, file="dist_or_ochiai_mx.RData")
load("dist_or_ochiai_mx.RData")
load("dist_or_sorensen_mx.RData")

######################################################### 2) Obtention des coordonnées factorielles

# Embedding ----
load("coord_fact_embedding.RData")

# Matrice de Jaccard avec transport optimal ----

# matrice_similarite <- read.csv("agro_to.csv") 
# matrice_similarite_t <- matrice_similarite %>%
#   select(-c(Acute.bronchospasm, 
#             Nonspecific.abnormal.findings.on.radiological.and.other.examination.of.other.intrathoracic.organs..echocardiogram..etc.))
# matrice_similarite_t <- t(matrice_similarite_t)
# matrice_similarite_t  <- as.data.frame(matrice_similarite_t)
# colnames(matrice_similarite_t ) <- matrice_similarite_t [1, ]
# matrice_similarite_t  <- matrice_similarite_t[-1, ]
# matrice_similarite_t_numeric <- matrice_similarite_t %>%
#   mutate(across(everything(), ~ as.numeric(as.character(.))))
# matrice_distance_apres_OT <- as.data.frame(1 - matrice_similarite_t_numeric)
# dist_ot_jaccard <- matrice_distance_apres_OT
# save(dist_ot_jaccard, file = "dist_ot_jaccard.RData") 

load("dist_ot_jaccard.RData")
dist_ot_jaccard_mx <- as.matrix(dist_ot_jaccard)

# ACM ----

 dta_acm <- data.frame(lapply(phenotype_maladie_s_c, as.factor))
 res.acm<- MCA(dta_acm, ind.sup = c(6126:7089), graph = TRUE, ncp=384)
 save(res.acm, file = "res.acm.RData")
 col_acm <- res.acm$ind$coord
 row_acm <- res.acm$ind.sup$coord
 coord_fact_acm2 <- rbind(col_acm, row_acm)
 coord_fact_acm2 <- as.data.frame(coord_fact_acm2)
 rm(coord_fact_acm)
 save(coord_fact_acm2, file = "coord_fact_acm2.RData")
 load("res.acm.RData")
 load("coord_fact_acm2.RData")
 
# ACM_5_dim <- readRDS("ACM_5_dim.rds")
# col_acm <- ACM_5_dim$ind$coord
# row_acm <- ACM_5_dim$ind.sup$coord
# coord_fact_acm <- rbind(col_acm, row_acm)
# coord_fact_acm <- as.data.frame(coord_fact_acm)
# save(coord_fact_acm, file = "coord_fact_acm.RData")
 load("coord_fact_acm.RData")
 dim(mandale$ind$coord)
 

# MDS unfolding ----

options(timeout = 6000)      # Augmenter le timeout à 10 minutes
execute_mds_unfolding <- function(distance_matrix, name) {
  coord_fact_mds <- smacofRect(delta = distance_matrix, 
                           ndim = 384, 
                           type = "ratio", 
                           itmax = 30, 
                           eps = 1e-6, 
                           verbose = TRUE)
  save(coord_fact_mds, file = paste0("coord_fact_mds_", name, ".RData"))
  return(coord_fact_mds)
}

rm(result_sorensen)
load("mds_result_sorensen_384.RData")
head(result_sorensen$conf.row[1:2, 1:2])
head(coord_fact_mds_sorensen$conf.row[1:2, 1:2])

coord_fact_mds_sorensen <- execute_mds_unfolding(dist_or_sorensen_mx, "sorensen")
coord_fact_mds_ochiai <- execute_mds_unfolding(dist_or_ochiai_mx, "ochiai")
coord_fact_mds_ot_jaccard <- execute_mds_unfolding(dist_ot_jaccard_mx, "ot_jaccard")

head(dist_or_ochiai_mx)
head(dist_or_jaccard_mx)
dist_or_jaccard_mx <- as.matrix(dist_or_jaccard)
dist_or_jaccard_mx <- apply(dist_or_jaccard_mx, 2, as.numeric)
coord_fact_mds_jaccard <- execute_mds_unfolding(dist_or_jaccard_mx, "jaccard")

save(coord_fact_mds_sorensen, "coord_fact_mds_sorensen.RData")
save(coord_fact_mds_ochiai, "coord_fact_mds_ochiai.RData")
save(coord_fact_mds_ot_jaccard, "coord_fact_mds_ot_jaccard.RData")
save(coord_fact_mds_jaccard, "coord_fact_mds_jaccard.RData")

rm(coord_fact_mds_jaccard)
rm(coord_fact_mds)
load("coord_fact_mds_sorensen.RData")
coord_fact_mds_sorensen <- coord_fact_mds
load("coord_fact_mds_ochiai.RData")
coord_fact_mds_ochiai <- coord_fact_mds
load("coord_fact_mds_jaccard.RData")
coord_fact_mds_jaccard <- coord_fact_mds
load("coord_fact_mds_ot_jaccard.RData")
coord_fact_mds_ot_jaccard <- coord_fact_mds
rm(coord_fact_mds)
# AFM sur maladies complexes ----

maladies_acm <- res.acm$ind.sup$coord
maladies_sorensen <- as.data.frame(coord_fact_mds_sorensen$conf.row)
maladies_ochiai <- as.data.frame(coord_fact_mds_ochiai$conf.row)
maladies_jaccard <- as.data.frame(coord_fact_mds_jaccard$conf.row)
maladies_ot_jaccard <- as.data.frame(coord_fact_mds_ot_jaccard$conf.row)
head(coord_fact_embedding)
maladies_embedding <- as.data.frame(coord_fact_embedding[6126:7089, ])

data_pour_afm <- cbind(maladies_acm, 
                       maladies_sorensen, 
                       maladies_ochiai, 
                       maladies_jaccard, 
                       maladies_ot_jaccard,
                       maladies_embedding)

# Est ce que l'embedding sur 964 lignes va donner la même chose que sur 7089 lignes ?

data_pour_afm <- as.data.frame(data_pour_afm)
save(data_pour_afm, file = "data_pour_afm.RData")
res.mfa <- MFA(data_pour_afm, 
               group = c(384, 384, 384, 384, 384, 384),
               type = c("c", "c", "c" , "c", "c", "c"),  # toutes les variables sont continues
               ncp = 2,             # nombre de dimensions à retenir
               name.group = c("ACM", "Sorensen", "Ochiai", "Jaccard", "Jaccard_OT","Embedding"))
save(res.mfa, file = "res.mfa.RData")

# Permutations ----
# Optimized permutation analysis function
# Optimized Sorensen distance function using vectorization
sorensen_distance <- function(x, y) {
  intersection <- sum(x & y)
  union <- sum(x) + sum(y)
  1 - (2 * intersection) / union
}

# Vectorized distance calculation for matrix
calc_distances <- function(matrix1, vector1) {
  intersections <- matrix1 %*% vector1  # Calcule tous les points d'intersection d'un coup
  sums1 <- rowSums(matrix1)  # Somme des lignes
  sum2 <- sum(vector1)  # Somme du vecteur
  1 - (2 * intersections) / (sums1 + sum2)
}

permut_all <- function(dta, n_perm, line_idx) {
  # Convert data to matrix format for faster operations
  dta <- as.matrix(dta)
  # Extract and store the row to permute
  row <- dta[line_idx, ]
  original_rowname <- rownames(dta)[line_idx]
  # Pre-allocate vector for minimum distances
  pmin <- numeric(n_perm)
  # Convert other rows to numeric matrix once
  other_rows <- matrix(as.numeric(dta[-line_idx, ]), 
                       nrow = nrow(dta) - 1,
                       byrow = FALSE)
  # Perform permutations efficiently
  for(i in 1:n_perm) {
    permuted_row <- sample(row, length(row), replace = FALSE)
    distances <- calc_distances(other_rows, permuted_row)
    pmin[i] <- min(distances)
  }
  return(c(rowname = original_rowname, threshold = min(pmin)))
}

# Optimized multiple rows processing
process_multiple_rows <- function(dta, n_perm, start_idx, end_idx) {
  n_rows <- end_idx - start_idx + 1
  results <- matrix(nrow = n_rows, ncol = 2)
  # Process in batches for better performance
  for(i in 1:n_rows) {
    curr_idx <- start_idx + i - 1
    results[i, ] <- permut_all(dta, n_perm, curr_idx)
  }
  results_df <- as.data.frame(results)
  colnames(results_df) <- c("rowname", "threshold")
  return(results_df)
}

tableau_permutation <- process_multiple_rows(phenotype_maladie_s_c, n_perm = 1, start_idx = 6126, end_idx = 7089)
save(tableau_permutation, file = "tableau_permutation.RData")

load ("tableau_permutation.RData")
# ensuite une fois que j'ai comparé mes méthodes entre elle j'obtiens une matrice d'assignation 
# avec les permutations qui font que j'associe un grand nombre de maladies simples à chaque 
# maladie complexe

# comment restreindre mon pool de maladies simples associées à chaque maladie complexe ? et 
# comment analyser les différents résultats obtenus des acm, embedding et autres distances

# voir si le clustering offre les mêmes qualités, les mêmes sorties que la méthode des permutations
# avec un pool restreint de maladies simples

# l'acm met les maladies complexes dans l'espace des maladies simples
# si on a des résultats différents c'est qu'on a des relations, des structures entre les maladies
# complexes qui sont captées dans les autres méthodes qu'on a pas dans l'acm

# donc potentiellement on a des maladies complexes qui sont associées entre elles et qui ne sont
# pas forcément liées dans l'acm


# Matrice d'assignation sorensen ----

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

association_df <- create_association_matrix(dist_or_sorensen_mx, tableau_permutation)
rowSums(association_df)
matrice_assignation_sorensen <- association_df 
save(matrice_assignation_sorensen, file="matrice_assignation_sorensen.RData")
rm(association_df)
rm(matrice_assignation_sorensen)

load("matrice_assignation_sorensen.RData")
rowSums(matrice_assignation_sorensen)
colSums(matrice_assignation_sorensen)

# Clustering ----

row_sorensen = coord_fact_mds_sorensen$conf.row
row_sorensen <- as.data.frame(row_sorensen)
col_sorensen = coord_fact_mds_sorensen$conf.col
col_sorensen <- as.data.frame(col_sorensen)
data_classif_sorensen <- rbind(row_sorensen, col_sorensen)
save(data_classif_sorensen, file="data_classif_sorensen.RData")

row_ochiai = coord_fact_mds_ochiai$conf.row
row_ochiai <- as.data.frame(row_ochiai)
col_ochiai = coord_fact_mds_ochiai$conf.col
col_ochiai  <- as.data.frame(col_ochiai )
data_classif_ochiai <- rbind(row_ochiai, col_ochiai)
save(data_classif_ochiai, file="data_classif_ochiai.RData")

row_jaccard = coord_fact_mds_jaccard$conf.row
row_jaccard <- as.data.frame(row_jaccard)
col_jaccard = coord_fact_mds_jaccard$conf.col
col_jaccard <- as.data.frame(col_jaccard)
data_classif_jaccard <- rbind(row_jaccard, col_jaccard)
save(data_classif_jaccard, file="data_classif_jaccard.RData")

row_ot_jaccard = coord_fact_mds_ot_jaccard$conf.row
row_ot_jaccard <- as.data.frame(row_ot_jaccard)
col_ot_jaccard = coord_fact_mds_ot_jaccard$conf.col
col_ot_jaccard <- as.data.frame(col_ot_jaccard)
data_classif_ot_jaccard <- rbind(row_ot_jaccard, col_ot_jaccard)
save(data_classif_ot_jaccard, file="data_classif_ot_jaccard.RData")
save(data_classif_embedding, file="data_classif_embedding.RData")
save(data_classif_acm, file="data_classif_acm.RData")

load("data_classif_sorensen.RData")
load("data_classif_ochiai.RData")
load("data_classif_jaccard.RData")
load("data_classif_ot_jaccard.RData")
load("coord_fact_embedding.RData")
load("coord_fact_acm2.RData")
data_classif_embedding <- as.data.frame(coord_fact_embedding)
data_classif_acm <- as.data.frame(coord_fact_acm2)

names(data_classif_sorensen) <- paste0("D", 1:ncol(data_classif_sorensen))
names(data_classif_ochiai) <- paste0("D", 1:ncol(data_classif_ochiai))
names(data_classif_jaccard) <- paste0("D", 1:ncol(data_classif_jaccard))
names(data_classif_ot_jaccard) <- paste0("D", 1:ncol(data_classif_ot_jaccard))
names(data_classif_embedding) <- paste0("D", 1:ncol(data_classif_embedding))
names(data_classif_acm) <- paste0("D", 1:ncol(data_classif_acm))

mc_indices <- 6126:7089

library(cluster)
library(factoextra)
silhouette_sorensen  <- fviz_nbclust(data_classif_sorensen, 
                                     kmeans, 
                                     method = "silhouette", 
                                     k.max = 400)
silhouette_sorensen
hc_sorensen <- HCPC(data_classif_sorensen,
           nb.clust = -1)
hc_sorensen <- hc
save(hc_sorensen, file="hc_sorensen.RData")
rm(hc)

silhouette_ochiai  <- fviz_nbclust(data_classif_ochiai, 
                                   kmeans, 
                                   method = "silhouette", 
                                   k.max = 400)
silhouette_ochiai
hc_ochiai <- HCPC(data_classif_ochiai,
                  nb.clust = -1)


silhouette_jaccard <- fviz_nbclust(data_classif_jaccard, 
                                   kmeans, 
                                   method = "silhouette", 
                                   k.max = 400)
silhouette_jaccard
hc_jaccard <- HCPC(data_classif_jaccard,
                   nb.clust = -1)

silhouette_ot_jaccard<- fviz_nbclust(data_classif_ot_jaccard, 
                                     kmeans, 
                                     method = "silhouette", 
                                     k.max = 400)
silhouette_ot_jaccard
hc_ot_jaccard <- HCPC(data_classif_ot_jaccard,
                      nb.clust = -1)

silhouette_embedding- fviz_nbclust(data_classif_embedding, 
                                   kmeans, 
                                   method = "silhouette", 
                                   k.max = 400)
silhouette_embedding
hc_embedding <- HCPC(data_classif_embedding,
                     nb.clust = -1)


silhouette_acm <- fviz_nbclust(data_classif_acm, 
                               kmeans, 
                               method = "silhouette", 
                               k.max = 400)
silhouette_acm
hc_acm <- HCPC(data_classif_acm,
               nb.clust = -1)

# Création d'une liste pour stocker tous les résultats
resultats_classification <- list(
  hc_sorensen = hc_sorensen,
  hc_ochiai = hc_ochiai,
  hc_jaccard = hc_jaccard,
  hc_ot_jaccard = hc_ot_jaccard,
  hc_embedding = hc_embedding,
  hc_acm = hc_acm
)

# Sauvegarder les résultats
save(resultats_classification, file = "resultats_classification.RData")

# Création des nouvelles matrices d'assignation à partir du clustering ----

hc_sorensen$desc.var
hc_sorensen$desc.axes
hc_sorensen$desc.ind
cluster_sorensen <- hc_sorensen$data.clust
str(cluster_sorensen)
summary(cluster_sorensen)

cluster_ochiai <- hc_ochiai$data.clust
cluster_jaccard <- hc_jaccard$data.clust
cluster_ot_jaccard <- hc_ot_jaccard$data.clust
cluster_embedding <- hc_embedding$data.clust
cluster_acm <- hc_acm$data.clust

# Convert clustering results to assignment matrix
create_assignment_matrix <- function(cluster_data) {
  # Extract cluster assignments and row names
  cluster_assignments <- cluster_data[,385]
  all_names <- rownames(cluster_data)
  # Extract names of simple diseases (from index 965 to end)
  ms_names <- all_names[965:length(all_names)]
  mc_names <- all_names[1:964]
  # Create empty matrix
  n_mc <- 964
  n_ms <- 7089 - 964
  assignment_matrix <- matrix(0, nrow=n_mc, ncol=n_ms)
  # Fill matrix - 1 if MC and MS are in same cluster
  for(i in 1:n_mc) {
    for(j in 1:n_ms) {
      if(cluster_assignments[i] == cluster_assignments[j+964]) {
        assignment_matrix[i,j] <- 1
      }
    }
  }
  # Add row and column names
  rownames(assignment_matrix) <- mc_names
  colnames(assignment_matrix) <- ms_names
  return(assignment_matrix)
}
# Use function
assignment_matrix <- create_assignment_matrix(cluster_sorensen)
assignment_matrix <- as.data.frame(assignment_matrix)

assignment_matrix_ochiai <- create_assignment_matrix(cluster_ochiai)
assignment_matrix_jaccard <- create_assignment_matrix(cluster_jaccard)
assignment_matrix_ot_jaccard <- create_assignment_matrix(cluster_ot_jaccard)
assignment_matrix_embedding <- create_assignment_matrix(cluster_embedding)
assignment_matrix_acm <- create_assignment_matrix(cluster_acm)

assignment_matrix_ochiai <- as.data.frame(assignment_matrix_ochiai)
assignment_matrix_jaccard <- as.data.frame(assignment_matrix_jaccard)
assignment_matrix_ot_jaccard <- as.data.frame(assignment_matrix_ot_jaccard)
assignment_matrix_embedding <- as.data.frame(assignment_matrix_embedding)
assignment_matrix_acm <- as.data.frame(assignment_matrix_acm)

# Création d'une liste pour stocker toutes les matrices d'assignation
matrices_assignation_classif <- list(
  assignment_matrix_classif_sorensen = assignment_matrix_sorensen,
  assignment_matrix_classif_ochiai = assignment_matrix_ochiai,
  assignment_matrix_classif_jaccard = assignment_matrix_jaccard,
  assignment_matrix_classif_ot_jaccard = assignment_matrix_ot_jaccard,
  assignment_matrix_classif_embedding = assignment_matrix_embedding,
  assignment_matrix_classif_acm = assignment_matrix_acm
)

# Sauvegarde des matrices d'assignation
save(matrices_assignation_classif, file = "matrices_assignation.RData")

# Comparaison des matrices d'assignation de sorensen entre elles ----

# Fonction pour compter les associations par maladie complexe
matrices_list <- list(
  sorensen_permut = matrice_assignation_sorensen,
  sorensen_mds_classif = assignment_matrix
)
count_associations <- function(row) {
  sum(row == 1)
}
# Initialiser une matrice pour stocker les résultats
result <- data.frame(row.names = rownames(matrice_assignation_sorensen))
# Fusionner toutes les matrices pour la comparaison
merged_matrix <- Reduce(`+`, matrices_list)
# Nombre total de matrices
n <- length(matrices_list)
# Remplir les colonnes pour chaque niveau d'association
for (i in n:0) {
  col_name <- ifelse(i == n, 
                     "Association_commune", 
                     ifelse(i == 0, 
                            "Aucune_association_commune", 
                            paste0("Association_commune_", i, "_matrices")))
  # Comparer et remplir les colonnes selon le nombre d'associations
  result[[col_name]] <- apply(merged_matrix, 1, function(x) sum(x == i))
}
# Vérifier les résultats
head(result)

rowSums(assignment_matrix[1:10,])
rowSums(matrice_assignation_sorensen[1,])

#2283 -40 = 2243

# Reste de la comparaison des matrices d'assignation ----

# Comparaison des matrices d'assignation entre elles
# Liste des matrices à comparer
matrices_list <- list(
  sorensen = assignment_matrix,
  ochiai = assignment_matrix_ochiai,
  jaccard = assignment_matrix_jaccard,
  ot_jaccard = assignment_matrix_ot_jaccard,
  embedding = assignment_matrix_embedding,
  acm = assignment_matrix_acm
)
# Fonction pour compter les associations par maladie complexe
count_associations <- function(row) {
  sum(row == 1)
}
# Initialiser une matrice pour stocker les résultats
result <- data.frame(row.names = rownames(assignment_matrix))
# Fusionner toutes les matrices pour la comparaison
merged_matrix <- Reduce(`+`, matrices_list)

# Nombre total de matrices
n <- length(matrices_list)
# Remplir les colonnes pour chaque niveau d'association
for (i in n:0) {
  col_name <- case_when(
    i == n ~ "Association_commune",
    i == 0 ~ "Aucune_association_commune",
    TRUE ~ paste0("Association_commune_", i, "_matrices")
  )
  # Comparer et remplir les colonnes selon le nombre d'associations
  result[[col_name]] <- apply(merged_matrix, 1, function(x) sum(x == i))
}
# Sauvegarder les résultats
save(result, file = "resultats_comparaison_matrices.RData")
# Examiner quelques résultats
head(result)
summary(result)
# Vérifier les nombres d'associations pour quelques maladies
print("Nombre d'associations par méthode pour les 5 premières MC:")
for(name in names(matrices_list)) {
  cat("\n", name, ":", rowSums(matrices_list[[name]][1:5,]))
}


