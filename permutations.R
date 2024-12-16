rm(list = ls())


# Calcul de la distance de Sørensen
# Fonction de distance personnalisée
sorensen_distance <- function(x, y) {
  1 - (2 * sum(x & y)) / (sum(x) + sum(y))
}


load(file = 'Data/data_clean0.RData')
load(file = 'Data/distance_origin_sorensen.RData')
n_perm <- 100
row <- or_df0[6127,]
permutation_list <- vector("list", n_perm)
distances_list <- vector("list", n_perm)
original_rowname <- rownames(row)
original_colnames <- colnames(or_df0)


for(i in 1:n_perm) {
  permuted_row <- sample(row, length(row), replace = FALSE)
  new_df <- permuted_row # créé matrice 1 seule ligne avec ligne permutée 
  colnames(new_df) <- original_colnames
  rownames(new_df) <- original_rowname
  permutation_list[[i]] <- new_df   # ajoute ma ligne permutée en dernier indice de ma liste de permutation
  temp_df <- rbind(or_df0[1:6125,], new_df) # créé temporairement une matrice les 6125 maladies 
  # simples et une maladie complexe 
  distances <- apply(temp_df[1:6125,], 1, function(row) sorensen_distance(as.numeric(row), as.numeric(new_df)))
  distances_list[[i]] <- distances # je stocke la matrice de distances à 1 ligne pour ma maladie simple
  # dans une liste de distances
}



test <- t(as.matrix(distances_list[[1]]))


plot(density(as.matrix(dist_or_sorensen_mx[2,])),col='red')
for(i in 1:10){
  lines(density(t(as.matrix(distances_list[[i]]))))
}



# Distribution des mc selon leur nombre de phénotypes associés

somPh <- as.matrix(apply(or_df0[6126:7091,],1,sum))

library(tidyverse)

as.data.frame(somPh) %>% 
  ggplot()+
  aes(x=V1)+
  geom_density(stat = 'density',linewidth = 1)+
  geom_vline(xintercept=330,col='blue')+
  geom_vline(xintercept=108,col='blue')+
  geom_vline(xintercept=67,col='blue')+
  geom_vline(xintercept=7,col='blue')+
  labs(title ='Nombre de phénotypes associés aux maladies complexes' )+
  xlab('Nombre de phénotypes')+
  ylab('Densité')+
  theme_classic()



# Définir les sommes cibles
target_sums <- c(330, 108, 67, 7)

# Initialiser un vecteur pour stocker les noms correspondants
matching_rows <- c()

# Parcourir chaque somme cible
for (target in target_sums) {
  # Trouver la première ligne correspondant à la somme cible
  matching_row <- rownames(or_df0[6126:7091, ])[which(apply(or_df0[6126:7091, ], 1, sum) == target)[1]]
  # Ajouter le nom trouvé au vecteur
  matching_rows <- c(matching_rows, matching_row)
}

# Afficher les noms de lignes correspondant aux sommes cibles
matching_rows


# test sur ces 4 maladies
# [1] "Congenital anomaly of fingers/toes"  330 -----------------------------------


n_perm <- 10
row <- or_df0[c("Congenital anomaly of fingers/toes"),]
permutation_list <- vector("list", n_perm)
distances_list <- vector("list", n_perm)
original_rowname <- rownames(row)
original_colnames <- colnames(or_df0)


for(i in 1:n_perm) {
  permuted_row <- sample(row, length(row), replace = FALSE)
  new_df <- permuted_row # créé matrice 1 seule ligne avec ligne permutée 
  colnames(new_df) <- original_colnames
  rownames(new_df) <- original_rowname
  permutation_list[[i]] <- new_df   # ajoute ma ligne permutée en dernier indice de ma liste de permutation
  temp_df <- rbind(or_df0[1:6125,], new_df) # créé temporairement une matrice les 6125 maladies 
  # simples et une maladie complexe 
  distances <- apply(temp_df[1:6125,], 1, function(row) sorensen_distance(as.numeric(row), as.numeric(new_df)))
  distances_list[[i]] <- distances # je stocke la matrice de distances à 1 ligne pour ma maladie simple
  # dans une liste de distances
}


plot(density(as.matrix(dist_or_sorensen_mx[c("Congenital anomaly of fingers/toes"),])),col='red',main = 'Congenital anomaly of fingers/toes 330 phenotypes')
pmin <- c()
for(i in 1:10){
  pmin <- append(pmin,min(distances_list[[i]]))
  lines(density(t(as.matrix(distances_list[[i]]))))
}

seuil <- min(pmin)

abd <- as.data.frame(dist_or_sorensen_mx[c("Congenital anomaly of fingers/toes"),])
# Sélectionner toutes les distances inférieures au seuil
selected_distances <- abd %>% 
  filter(`dist_or_sorensen_mx[c("Congenital anomaly of fingers/toes"), ]`<=seuil)
# on a 153 gènes sélectionné


# [2] "Other specified congenital anomalies of nervous system" 108 -----------------------------------


n_perm <- 10
row <- or_df0[c("Other specified congenital anomalies of nervous system"),]
permutation_list <- vector("list", n_perm)
distances_list <- vector("list", n_perm)
original_rowname <- rownames(row)
original_colnames <- colnames(or_df0)


for(i in 1:n_perm) {
  permuted_row <- sample(row, length(row), replace = FALSE)
  new_df <- permuted_row # créé matrice 1 seule ligne avec ligne permutée 
  colnames(new_df) <- original_colnames
  rownames(new_df) <- original_rowname
  permutation_list[[i]] <- new_df   # ajoute ma ligne permutée en dernier indice de ma liste de permutation
  temp_df <- rbind(or_df0[1:6125,], new_df) # créé temporairement une matrice les 6125 maladies 
  # simples et une maladie complexe 
  distances <- apply(temp_df[1:6125,], 1, function(row) sorensen_distance(as.numeric(row), as.numeric(new_df)))
  distances_list[[i]] <- distances # je stocke la matrice de distances à 1 ligne pour ma maladie simple
  # dans une liste de distances
}

plot(density(as.matrix(dist_or_sorensen_mx[c("Other specified congenital anomalies of nervous system"),])),col='red',main = 'Other specified congenital anomalies of nervous system : 108 phenotypes', sub = '1083 simple diseases kept')
pmin <- c()
for(i in 1:10){
  pmin <- append(pmin,min(distances_list[[i]]))
  lines(density(t(as.matrix(distances_list[[i]]))))
}

seuil <- min(pmin)

library(tidyverse)

abd <- as.data.frame(dist_or_sorensen_mx[c("Other specified congenital anomalies of nervous system"),])
# Sélectionner toutes les distances inférieures au seuil
selected_distances <- abd %>% 
  filter(`dist_or_sorensen_mx[c("Other specified congenital anomalies of nervous system"), ]`<=seuil)
# on a 1041 gènes sélectionné



# regarder les phénotypes
library(reshape)

liste_assign <- c(rownames(selected_distances))

phenotype <- or_df0[liste_assign,]

# Combiner les deux dataframes pour comparaison
comparison <- phenotype %>%
  mutate(across(everything(), ~ ifelse(. == 1 & new_df[1, cur_column()]==1, 1, 0)))  # Réattacher les noms des maladies


nb1 <- as.data.frame(apply(comparison,1,sum))

nb <- nb1 %>% filter(`apply(comparison, 1, sum)`!=0)
# Résultat
comparison

# Comparer et compter les phénotypes
result <- phenotype %>%
  rowwise() %>% # Permet de traiter chaque ligne individuellement
  mutate(
    Commun = sum(c_across() == 1 & new_df[1, ]==1),        # Compter les phénotypes communs
    Autre = sum(c_across() != new_df[1, ])      # Compter les phénotypes différents
  ) %>%
  ungroup() %>%
  select( rownames(phenotype),Commun, Autre) # Garder uniquement les colonnes nécessaires

# Afficher le résultat
result


# [3] "Congenital osteodystrophies" 67 -----------------------------------


n_perm <- 10
row <- or_df0[c("Congenital osteodystrophies"),]
permutation_list <- vector("list", n_perm)
distances_list <- vector("list", n_perm)
original_rowname <- rownames(row)
original_colnames <- colnames(or_df0)


for(i in 1:n_perm) {
  permuted_row <- sample(row, length(row), replace = FALSE)
  new_df <- permuted_row # créé matrice 1 seule ligne avec ligne permutée 
  colnames(new_df) <- original_colnames
  rownames(new_df) <- original_rowname
  permutation_list[[i]] <- new_df   # ajoute ma ligne permutée en dernier indice de ma liste de permutation
  temp_df <- rbind(or_df0[1:6125,], new_df) # créé temporairement une matrice les 6125 maladies 
  # simples et une maladie complexe 
  distances <- apply(temp_df[1:6125,], 1, function(row) sorensen_distance(as.numeric(row), as.numeric(new_df)))
  distances_list[[i]] <- distances # je stocke la matrice de distances à 1 ligne pour ma maladie simple
  # dans une liste de distances
}

plot(density(as.matrix(dist_or_sorensen_mx[c('Congenital osteodystrophies'),])),col='red',main = 'Congenital osteodystrophies 67 phenotypes')
pmin <- c()
for(i in 1:10){
  pmin <- append(pmin,min(distances_list[[i]]))
  lines(density(t(as.matrix(distances_list[[i]]))))
}
    
seuil <- min(pmin)

abd <- as.data.frame(dist_or_sorensen_mx[c('Congenital osteodystrophies'),])
# Sélectionner toutes les distances inférieures au seuil
selected_distances <- abd %>% 
  filter(`dist_or_sorensen_mx[c("Congenital osteodystrophies"), ]`<=seuil)
# on a 530 gènes sélectionné


# [4] "Abdominal pain"  7 -----------------------------------


n_perm <- 10
row <- or_df0[c("Abdominal pain"),]
permutation_list <- vector("list", n_perm)
distances_list <- vector("list", n_perm)
original_rowname <- rownames(row)
original_colnames <- colnames(or_df0)


for(i in 1:n_perm) {
  permuted_row <- sample(row, length(row), replace = FALSE)
  new_df <- permuted_row # créé matrice 1 seule ligne avec ligne permutée 
  colnames(new_df) <- original_colnames
  rownames(new_df) <- original_rowname
  permutation_list[[i]] <- new_df   # ajoute ma ligne permutée en dernier indice de ma liste de permutation
  temp_df <- rbind(or_df0[1:6125,], new_df) # créé temporairement une matrice les 6125 maladies 
  # simples et une maladie complexe 
  distances <- apply(temp_df[1:6125,], 1, function(row) sorensen_distance(as.numeric(row), as.numeric(new_df)))
  distances_list[[i]] <- distances # je stocke la matrice de distances à 1 ligne pour ma maladie simple
  # dans une liste de distances
}

plot(density(as.matrix(dist_or_sorensen_mx[c("Abdominal pain"),])),col='red',main = 'Abdominal pain : 7 phenotypes', sub = '10 simple diseases kept')
pmin <- c()
for(i in 1:10){
  pmin <- append(pmin,min(distances_list[[i]]))
  lines(density(t(as.matrix(distances_list[[i]]))))
}

seuil <- min(pmin)

abd <- as.data.frame(dist_or_sorensen_mx[c("Abdominal pain"),])
# Sélectionner toutes les distances inférieures au seuil
selected_distances <- abd %>% 
  filter(`dist_or_sorensen_mx[c("Abdominal pain"), ]`<=seuil)
# on a 4 gènes sélectionné








