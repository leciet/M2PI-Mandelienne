# Création du dataframe initial ----

# setwd("~/S9_ACO/Projet 2 mois/M2PI-Mendelienne/Matrices Jules sans hiérarchie")
# Phe_OMIM <- read.csv('Profils_Base_OMIM.csv')
# Phe_Phecode <- read.csv('Profils_Base_Phecodes.csv')
# 
# rownames(Phe_OMIM) <- Phe_OMIM[,1]
# rownames(Phe_Phecode) <- Phe_Phecode[,1]
# Phe_OMIM <- Phe_OMIM[,-1]
# Phe_Phecode <- Phe_Phecode[,-1]
# 
# colonnes_communes <- intersect(colnames(Phe_OMIM), colnames(Phe_Phecode))
# Phe_OMIM_communes <- Phe_OMIM[, colonnes_communes, drop = FALSE]
# Phe_Phecode_communes <- Phe_Phecode[, colonnes_communes, drop = FALSE]
# mc_ms_communes <- rbind(Phe_OMIM_communes, Phe_Phecode_communes)
# 
# library(dplyr)
# mc_ms_communes_filtre_row <- mc_ms_communes %>% 
#   filter(rowSums(mc_ms_communes)!=0)
# mc_ms_communes_filtre_col_row <- mc_ms_communes_filtre_row %>%
#   select(where(~sum(.) != 0))
# mc_ms_communes_filtre_col_row <- mc_ms_communes_filtre_row[, colSums(mc_ms_communes_filtre_row) != 0]
# 
# phenotype_maladie_s_c <- mc_ms_communes_filtre_col_row
# phenotype_maladie_s_c <-  as.data.frame(phenotype_maladie_s_c)

# save(phenotype_maladie_s_c, file="phenotype_maladie_s_c.RData")

rm(list=ls())
load("phenotype_maladie_s_c.RData")
phenotype_maladie_s_c_binary <- as.matrix(phenotype_maladie_s_c == 1)

# Library ----
library(tidyverse)
library(FactoMineR)
library(Factoshiny)
library(ade4)
library(vegan)
library(reshape)
library(DataExplorer)
library(smacof)

######################################################### 1) Matrices de distance

# Autres distances ----

dist_or_sokal_sneath <- dist.binary(phenotype_maladie_s_c, method = 3) 
dist_or_sokal_sneath <- as.data.frame(as.matrix(dist_or_sokal_sneath))
dist_or_sokal_sneath <- dist_or_sokal_sneath[6103:7064, 1:6102]
dist_or_sokal_sneath_mx <- as.matrix(dist_or_sokal_sneath)

dist_or_sokal_michener <- dist.binary(phenotype_maladie_s_c, method = 2)
dist_or_sokal_michener <- as.data.frame(as.matrix(dist_or_sokal_michener))
dist_or_sokal_michener <- dist_or_sokal_michener[6103:7064, 1:6102]
dist_or_sokal_michener_mx <- as.matrix(dist_or_sokal_michener)

dist_or_rogers_tanimoto  <- dist.binary(phenotype_maladie_s_c, method = 4)
dist_or_rogers_tanimoto <- as.data.frame(as.matrix(dist_or_rogers_tanimoto))
dist_or_rogers_tanimoto <- dist_or_rogers_tanimoto[6103:7064, 1:6102]
dist_or_rogers_tanimoto_mx <- as.matrix(dist_or_rogers_tanimoto)

dist_or_hamann <- dist.binary(phenotype_maladie_s_c, method = 6) 
dist_or_hamann <- as.data.frame(as.matrix(dist_or_hamann))
dist_or_hamann <- dist_or_hamann[6103:7064, 1:6102]
dist_or_hamann_mx <- as.matrix(dist_or_hamann)


dist_cosine <- dist(phenotype_maladie_s_c, method = "cosine")
dist_cosine <- as.data.frame(as.matrix(dist_cosine))
dist_cosine <- dist_cosine[6103:7064, 1:6102]
dist_or_cosine_mx <- as.matrix(dist_cosine)

save(dist_or_cosine_mx, file="dist_or_cosine_mx.RData")
save(dist_or_sokal_sneath_mx, file="dist_or_sokal_sneath_mx.RData")
save(dist_or_sokal_michener_mx, file="dist_or_sokal_michener_mx.RData")
save(dist_or_rogers_tanimoto_mx, file="dist_or_rogers_tanimoto_mx.RData")
save(dist_or_hamann_mx, file="dist_or_hamann_mx.RData")

rm(dist_or_tanimoto_mx)
load("dist_or_sokal_sneath_mx.RData")
dist_or_sokal_sneath_mx<- dist_or_tanimoto_mx
rm(dist_or_tanimoto_mx)
load("dist_or_hamann_mx.RData")
dist_or_hamann_mx <- dist_or_dice_mx
rm(dist_or_dice_mx)
load("dist_or_rogers_tanimoto_mx.RData")
dist_or_rogers_tanimoto_mx <- dist_or_kulczynski_mx
rm(dist_or_kulczynski_mx)
load("dist_or_sokal_michener_mx.RData")
dist_or_sokal_michener_mx <- dist_or_hamming_mx
rm(dist_or_hamming_mx)
load("dist_or_cosine_mx.RData")

# Calculate Jaccard distance ----


index_ASCVD <- which(rownames(phenotype_maladie_s_c) == "ASCVD")
rownames(phenotype_maladie_s_c[6103,])
rownames(phenotype_maladie_s_c[6102,])

dist_or_jaccard <- dist.binary(phenotype_maladie_s_c, method = 1)
dist_or_jaccard <- as.data.frame(as.matrix(dist_or_jaccard))
dist_or_jaccard <- dist_or_jaccard[6103:7064, 1:6102]  # Extraction du sous-ensemble MC-MS
dist_or_jaccard_mx <- as.matrix(dist_or_jaccard)

# Calculate Sorensen distance ----

dist_or_sorensen <- dist.binary(phenotype_maladie_s_c, method = 5)
dist_or_sorensen <- as.data.frame(as.matrix(dist_or_sorensen))
dist_or_sorensen <- dist_or_sorensen[6103:7064, 1:6102]
dist_or_sorensen_mx <- as.matrix(dist_or_sorensen)

# Calculate Ochiai distance ----

dist_or_ochiai <- dist.binary(phenotype_maladie_s_c, method = 7)
dist_or_ochiai <- as.data.frame(as.matrix(dist_or_ochiai))
dist_or_ochiai <- dist_or_ochiai[6103:7064, 1:6102]
dist_or_ochiai_mx <- as.matrix(dist_or_ochiai)

# Save distance matrices ----

save(dist_or_jaccard_mx, file="dist_or_jaccard_mx.RData")
save(dist_or_sorensen_mx, file="dist_or_sorensen_mx.RData")
save(dist_or_ochiai_mx, file="dist_or_ochiai_mx.RData")

# Load distance matrices ----
load("dist_or_jaccard_mx.RData")
load("dist_or_ochiai_mx.RData")
load("dist_or_sorensen_mx.RData")
load("dist_or_sorensen_mx.RData")

######################################################### 2) Obtention des coordonnées factorielles

# Embedding ----

# Load necessary library
# if (!require("ontologyIndex")) install.packages("ontologyIndex", repos = "http://cran.us.r-project.org")
# library(ontologyIndex)

# Load the HPO file
# hpo_url <- "https://raw.githubusercontent.com/obophenotype/human-phenotype-ontology/master/hp.obo"
# hpo <- get_ontology(hpo_url, extract_tags = 'everything')
# name_ph <- hpo$name
# formatted_colnames <- gsub("\\.", ":", colnames(phenotype_maladie_s_c))
# colnames(phenotype_maladie_s_c) <- name_ph[formatted_colnames]

# transformation au format long 
# df_long <- melt(as.matrix(phenotype_maladie_s_c)) 
# colnames(df_long) <- c('maladie','phenotype','presence')
# liste_phenotype <- df_long %>% 
#   filter(presence == 1) %>% 
#   group_by(maladie) %>% 
#   summarise( liste_phenotype = paste('We are a group of individuals suffering from the same disease, but it manifests itself differently in each of us. Here is a list of the symptoms we may experience 
# ',paste(phenotype,collapse = ", "),sep = ' : '),.groups = 'drop')
# liste_phenotype <- unique(liste_phenotype)[,2]
# write_csv(x = liste_phenotype,file = 'Data/embedding_0hierarchie.csv')

# embedding par Massamba ici

# coord_fact_embedding <- read.csv("embeddings_Pheno_SBERT_sans_hierarchie.csv")
# save(coord_fact_embedding, file = "coord_fact_embedding.RData")

load("coord_fact_embedding.RData")

# Matrice de Jaccard avec transport optimal ----

sinkhorn_transport <- function(a, b, M, reg, numItermax=1000, stopThr=1e-9) {
  # a, b: distributions marginales
  # M: matrice de coût
  # reg: paramètre de régularisation
  # numItermax: nombre maximum d'itérations
  # stopThr: seuil de convergence
  n <- length(a)
  m <- length(b)
  # Initialisation
  u <- rep(1/n, n)
  v <- rep(1/m, m)
  K <- exp(-M/reg)
  
  # Boucle de Sinkhorn
  for(i in 1:numItermax) {
    u_prev <- u
    # Mise à jour de u et v
    u <- a/(K %*% v)
    v <- b/(t(K) %*% u)
    # Vérification de la convergence
    err <- sum(abs(u - u_prev))
    if(err < stopThr) break
  }
  # Calcul de la matrice de transport
  P <- diag(as.vector(u)) %*% K %*% diag(as.vector(v))
  return(P)
}

# Utilisation sur votre matrice Jaccard
# Normalisation des vecteurs marginaux
a <- rep(1/nrow(dist_or_jaccard_mx), nrow(dist_or_jaccard_mx))
b <- rep(1/ncol(dist_or_jaccard_mx), ncol(dist_or_jaccard_mx))

# Application de l'algorithme de Sinkhorn
transport_matrix <- sinkhorn_transport(
  a = a,
  b = b,
  M = dist_or_jaccard_mx,
  reg = 0.1  # paramètre de régularisation à ajuster
)

# Calcul de la matrice de distance avec transport optimal
ot_matrix <- dist_or_jaccard_mx * transport_matrix
ot_matrix <- as.matrix(ot_matrix)
ot_matrix <- as.data.frame(ot_matrix)

# ot_matrix_numeric <- ot_matrix %>%
#    mutate(across(everything(), ~ as.numeric(as.character(.))))
dist_ot_jaccard  <- as.data.frame(1 - ot_matrix)
dist_ot_jaccard_mx <- as.matrix(dist_ot_jaccard)
save(dist_ot_jaccard_mx, file="dist_ot_jaccard_mx.RData")

load("dist_ot_jaccard_mx.RData")

colonnes_communes <- intersect(colnames(Phe_OMIM), colnames(Phe_Phecode))
Phe_OMIM_communes <- Phe_OMIM[, colonnes_communes, drop = FALSE]

# ACM ----

dta_acm <- data.frame(lapply(phenotype_maladie_s_c, as.factor))
res.acm<- MCA(dta_acm, ind.sup = c(6103:7064), graph = TRUE, ncp=384)
save(res.acm, file = "res.acm.RData")
col_acm <- res.acm$ind$coord
row_acm <- res.acm$ind.sup$coord
coord_fact_acm <- rbind(col_acm, row_acm)
coord_fact_acm <- as.data.frame(coord_fact_acm)
save(coord_fact_acm, file = "coord_fact_acm.RData")
load("res.acm.RData")
load("coord_fact_acm.RData")

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

coord_fact_mds_sorensen <- execute_mds_unfolding(dist_or_sorensen_mx, "sorensen")
coord_fact_mds_ochiai <- execute_mds_unfolding(dist_or_ochiai_mx, "ochiai")
coord_fact_mds_jaccard <- execute_mds_unfolding(dist_or_jaccard_mx, "jaccard")
coord_fact_mds_ot_jaccard <- execute_mds_unfolding(dist_ot_jaccard_mx, "ot_jaccard")

save(coord_fact_mds_sorensen, "coord_fact_mds_sorensen.RData")
save(coord_fact_mds_ochiai, "coord_fact_mds_ochiai.RData")
save(coord_fact_mds_jaccard, "coord_fact_mds_jaccard.RData")
save(coord_fact_mds_ot_jaccard, "coord_fact_mds_ot_jaccard.RData")

rm(list=ls())
load("coord_fact_mds_sorensen.RData")
coord_fact_mds_sorensen <- coord_fact_mds
load("coord_fact_mds_ochiai.RData")
coord_fact_mds_ochiai <- coord_fact_mds
load("coord_fact_mds_jaccard.RData")
coord_fact_mds_jaccard <- coord_fact_mds
load("coord_fact_mds_ot_jaccard.RData")
coord_fact_mds_ot_jaccard <- coord_fact_mds
rm(coord_fact_mds)

load("coord_fact_mds_hamann.RData")
coord_fact_mds_hamann <- coord_fact_mds_dice
rm(coord_fact_mds_dice)
load("coord_fact_mds_rogers_tanimoto.RData")
coord_fact_mds_rogers_tanimoto <- coord_fact_mds_kulczynski
rm(coord_fact_mds_kulczynski)
load("coord_fact_mds_sokal_michener.RData")
coord_fact_mds_sokal_michener <- coord_fact_mds_hamming
rm(coord_fact_mds_hamming)
load("coord_fact_mds_sokal_sneath.RData")
coord_fact_mds_sokal_sneath <- coord_fact_mds_tanimoto
rm(coord_fact_mds_tanimoto)

# AFM sur maladies complexes ----

maladies_acm <- res.acm$ind.sup$coord
maladies_sorensen <- as.data.frame(coord_fact_mds_sorensen$conf.row)
maladies_ochiai <- as.data.frame(coord_fact_mds_ochiai$conf.row)
maladies_jaccard <- as.data.frame(coord_fact_mds_jaccard$conf.row)
maladies_ot_jaccard <- as.data.frame(coord_fact_mds_ot_jaccard$conf.row)
head(coord_fact_embedding)
maladies_embedding <- as.data.frame(coord_fact_embedding[6103:7064, ])

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

# Visualisation des densités de distance avec ou sans permutation ----
phenotype_maladie_s_c2 <- as.matrix(phenotype_maladie_s_c)
row <- phenotype_maladie_s_c2[6103, ]
other_rows <- matrix(as.numeric(phenotype_maladie_s_c2[-6103, ]), 
                     nrow = nrow(phenotype_maladie_s_c2) - 1,
                     byrow = FALSE)
calc_distances <- function(matrix1, vector1) {
  intersections <- matrix1 %*% vector1  # Calcule tous les points d'intersection d'un coup
  sums1 <- rowSums(matrix1)  # Somme des lignes
  sum2 <- sum(vector1)  # Somme du vecteur
  1 - (2 * intersections) / (sums1 + sum2)
}
distance_originale <- calc_distances(other_rows, row)
View(distance_originale)

permuted_row1 <- sample(row, length(row), replace = FALSE)
distance_permut_1 <- calc_distances(other_rows, permuted_row1)

permuted_row2 <- sample(row, length(row), replace = FALSE)
distance_permut_2 <- calc_distances(other_rows, permuted_row2)

permuted_row3 <- sample(row, length(row), replace = FALSE)
distance_permut_3 <- calc_distances(other_rows, permuted_row3)

permuted_row4 <- sample(row, length(row), replace = FALSE)
distance_permut_4 <- calc_distances(other_rows, permuted_row4)

permuted_row5 <- sample(row, length(row), replace = FALSE)
distance_permut_5 <- calc_distances(other_rows, permuted_row5)

library(ggplot2)

# Regrouper les distances dans un data.frame
distance_data <- data.frame(
  Distance = c(distance_originale, 
               distance_permut_1, 
               distance_permut_2, 
               distance_permut_3, 
               distance_permut_4, 
               distance_permut_5),
  Type = factor(c(
    rep("Originale", length(distance_originale)),
    rep("Permutation 1", length(distance_permut_1)),
    rep("Permutation 2", length(distance_permut_2)),
    rep("Permutation 3", length(distance_permut_3)),
    rep("Permutation 4", length(distance_permut_4)),
    rep("Permutation 5", length(distance_permut_5))
  ))
)

# Visualisation avec ggplot2
ggplot(distance_data, aes(x = Distance, color = Type, fill = Type)) +
  geom_density(alpha = 0.3, size = 1) +  # Courbes de densité avec transparence
  scale_color_manual(values = c("Originale" = "red", 
                                "Permutation 1" = "black", 
                                "Permutation 2" = "blue", 
                                "Permutation 3" = "green", 
                                "Permutation 4" = "purple", 
                                "Permutation 5" = "orange")) +
  scale_fill_manual(values = c("Originale" = "red", 
                               "Permutation 1" = "black", 
                               "Permutation 2" = "blue", 
                               "Permutation 3" = "green", 
                               "Permutation 4" = "purple", 
                               "Permutation 5" = "orange")) +
  labs(title = "Densité des distances (Originale vs 5 permutations)",
       x = "Distance (Sørensen-Dice)",
       y = "Densité") +
  theme_minimal() +
  theme(legend.title = element_blank())  # Retire le titre de la légende

# Permutations ----

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

# dta3 <- permut_all(phenotype_maladie_s_c, n_perm=5, line_idx=6103)
# dta3
# dta <- permut_all(phenotype_maladie_s_c, n_perm=1000, line_idx=6104)
# dta
# dta2 <- permut_all(phenotype_maladie_s_c, n_perm=5, line_idx=6104)
# dta2
# dt4 <- permut_all(phenotype_maladie_s_c, n_perm=2, line_idx=6105)
# dt4
# compte_par_ligne <- apply(phenotype_maladie_s_c, 1, function(ligne) sum(ligne != 1))
# print(compte_par_ligne)

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

tableau_permutation <- process_multiple_rows(phenotype_maladie_s_c, n_perm = 5, start_idx = 6103, end_idx = 7064)
tableau_permutation$threshold <- as.numeric(tableau_permutation$threshold)
head(tableau_permutation)
save(tableau_permutation, file="tableau_permutation.RData")

load ("tableau_permutation.RData")

# head(dist_or_sorensen[1:10, 1:10])
# 
# dist_or_sorensen <- as.data.frame(dist_or_sorensen_mx)
# compte_par_ligne <- apply(dist_or_sorensen, 1, function(ligne) sum(ligne != 1))
# print(compte_par_ligne)
# 
# all(rownames(dist_or_sorensen_mx) == tableau_permutation$rowname)
# dim(dist_or_sorensen_mx)
# dim(tableau_permutation)


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

matrice_assignation_sorensen_permutation <- create_association_matrix(dist_or_sorensen_mx, 
                                                                      tableau_permutation)
rowSums(matrice_assignation_sorensen_permutation)
save(matrice_assignation_sorensen_permutation, file="matrice_assignation_sorensen.RData")

matrice_filtrée <- matrice_assignation_sorensen_permutation[rowSums(matrice_assignation_sorensen_permutation) != 0, ]
View(matrice_filtrée)

load("matrice_assignation_sorensen_permutation.RData")
rowSums(matrice_assignation_sorensen_permutation)
colSums(matrice_assignation_sorensen_permutation)

# Clustering ----

row_sorensen = coord_fact_mds_sorensen$conf.row
row_sorensen <- as.data.frame(row_sorensen)
col_sorensen = coord_fact_mds_sorensen$conf.col
col_sorensen <- as.data.frame(col_sorensen)
data_classif_sorensen <- rbind(row_sorensen, col_sorensen)

row_ochiai = coord_fact_mds_ochiai$conf.row
row_ochiai <- as.data.frame(row_ochiai)
col_ochiai = coord_fact_mds_ochiai$conf.col
col_ochiai  <- as.data.frame(col_ochiai )
data_classif_ochiai <- rbind(row_ochiai, col_ochiai)

row_jaccard = coord_fact_mds_jaccard$conf.row
row_jaccard <- as.data.frame(row_jaccard)
col_jaccard = coord_fact_mds_jaccard$conf.col
col_jaccard <- as.data.frame(col_jaccard)

row_ot_jaccard = coord_fact_mds_ot_jaccard$conf.row
row_ot_jaccard <- as.data.frame(row_ot_jaccard)
col_ot_jaccard = coord_fact_mds_ot_jaccard$conf.col
col_ot_jaccard <- as.data.frame(col_ot_jaccard)
data_classif_ot_jaccard <- rbind(row_ot_jaccard, col_ot_jaccard)

data_classif_embedding <- as.data.frame(coord_fact_embedding)
section1 <- data_classif_embedding[1:6102, ]
section2 <- data_classif_embedding[6103:7064, ]
data_classif_embedding <- rbind(section2, section1)
rownames(data_classif_embedding) <- 1:nrow(data_classif_embedding)
save(data_classif_embedding, file = "data_classif_embedding.RData")

data_classif_acm <- as.data.frame(coord_fact_acm)
section1 <- data_classif_acm[1:6102, ]
section2 <- data_classif_acm[6103:7064, ]
data_classif_acm <- rbind(section2, section1)
rownames(data_classif_acm) <- 1:nrow(data_classif_acm)
save(data_classif_acm, file = "data_classif_acm.RData")

names(data_classif_sorensen) <- paste0("D", 1:ncol(data_classif_sorensen))
names(data_classif_ochiai) <- paste0("D", 1:ncol(data_classif_ochiai))
names(data_classif_jaccard) <- paste0("D", 1:ncol(data_classif_jaccard))
names(data_classif_ot_jaccard) <- paste0("D", 1:ncol(data_classif_ot_jaccard))
names(data_classif_embedding) <- paste0("D", 1:ncol(data_classif_embedding))
names(data_classif_acm) <- paste0("D", 1:ncol(data_classif_acm))

save(data_classif_sorensen, file="data_classif_sorensen.RData")
save(data_classif_sorensen, file="data_classif_sorensen.RData")
save(data_classif_jaccard, file="data_classif_jaccard.RData")
save(data_classif_ot_jaccard, file="data_classif_ot_jaccard.RData")
save(data_classif_embedding, file="data_classif_embedding.RData")
save(data_classif_acm, file="data_classif_acm.RData")

load("data_classif_sorensen.RData")
load("data_classif_ochiai.RData")
load("data_classif_jaccard.RData")
load("data_classif_ot_jaccard.RData")
load("data_classif_embedding.RData")
load("data_classif_acm.RData")

mc_indices <- 6103:7064

library(cluster)
library(factoextra)
silhouette_sorensen  <- fviz_nbclust(data_classif_sorensen, 
                                     kmeans, 
                                     method = "silhouette", 
                                     k.max = 700)
silhouette_sorensen
hc_sorensen <- HCPC(data_classif_sorensen,
                    nb.clust = -1)
save(hc_sorensen, file="hc_sorensen.RData")

silhouette_ochiai  <- fviz_nbclust(data_classif_ochiai, 
                                   kmeans, 
                                   method = "silhouette", 
                                   k.max = 700)
silhouette_ochiai
hc_ochiai <- HCPC(data_classif_ochiai,
                  nb.clust = -1)
save(hc_ochiai, file="hc_ochiai.RData")


silhouette_jaccard <- fviz_nbclust(data_classif_jaccard, 
                                   kmeans, 
                                   method = "silhouette", 
                                   k.max = 700)
silhouette_jaccard
hc_jaccard <- HCPC(data_classif_jaccard,
                   nb.clust = -1)
save(hc_jaccard, file="hc_jaccard.RData")

silhouette_ot_jaccard<- fviz_nbclust(data_classif_ot_jaccard, 
                                     kmeans, 
                                     method = "silhouette", 
                                     k.max = 700)
silhouette_ot_jaccard
hc_ot_jaccard <- HCPC(data_classif_ot_jaccard,
                      nb.clust = -1)
save(hc_ot_jaccard, file="hc_ot_jaccard.RData")

silhouette_embedding- fviz_nbclust(data_classif_embedding, 
                                   kmeans, 
                                   method = "silhouette", 
                                   k.max = 700)
silhouette_embedding
hc_embedding <- HCPC(data_classif_embedding,
                     nb.clust = -1)
save(hc_embedding, file="hc_embedding.RData")

silhouette_acm <- fviz_nbclust(data_classif_acm, 
                               kmeans, 
                               method = "silhouette", 
                               k.max = 700)
silhouette_acm
hc_acm <- HCPC(data_classif_acm,
               nb.clust = -1)
save(hc_acm, file="hc_acm.RData")

resultats_classification <- list(
  hc_sorensen = hc_sorensen,
  hc_ochiai = hc_ochiai,
  hc_jaccard = hc_jaccard,
  hc_ot_jaccard = hc_ot_jaccard,
  hc_embedding = hc_embedding,
  hc_acm = hc_acm
)

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
assignment_matrix_sorensen <- create_assignment_matrix(cluster_sorensen)
assignment_matrix_sorensen <- as.data.frame(assignment_matrix_sorensen)

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


