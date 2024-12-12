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

















