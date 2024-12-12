rm(list = ls())


# Calcul de la distance de Sørensen
# Fonction de distance personnalisée
sorensen_distance <- function(x, y) {
  1 - (2 * sum(x & y)) / (sum(x) + sum(y))
}


load(file = 'Data/data_clean0.RData')
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

plot(density(as.matrix(dist_or_sorensen_mx[2,])),col='red',main = 'Congenital anomaly of fingers/toes 330 phenotypes')
for(i in 1:10){
  lines(density(t(as.matrix(distances_list[[i]]))))
}


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

plot(density(as.matrix(dist_or_sorensen_mx[2,])),col='red',main = 'Other specified congenital anomalies of nervous system 108 phenotypes')
for(i in 1:10){
  lines(density(t(as.matrix(distances_list[[i]]))))
}



