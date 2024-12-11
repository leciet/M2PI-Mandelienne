library(ade4)

generate_permutations_and_distances <- function(row, or_df0, n_perm = 200) {
  permutation_list <- vector("list", n_perm)
  distances_list <- vector("list", n_perm)
  original_rowname <- rownames(row)
  original_colnames <- colnames(or_df0)
  
  for(i in 1:n_perm) {
    permuted_row <- sample(row, length(row), replace = FALSE)
    new_df <- data.frame(matrix(permuted_row, nrow=1)) # créé matrice 1 seule ligne avec ligne permutée 
    colnames(new_df) <- original_colnames
    rownames(new_df) <- original_rowname
    permutation_list[[i]] <- new_df   # ajoute ma ligne permutée en dernier indice de ma liste de permutation
    temp_df <- rbind(or_df0[1:6125,], permutation_list[[i]]) # créé temporairement une matrice les 6125 maladies 
    # simples et une maladie complexe 
    distances_sorensen <- as.matrix(dist.binary(temp_df,method = 5))[6126,1:6125]  # créé une matrice distance 
    distances_list[[i]] <- distances_sorensen # je stocke la matrice de distances à 1 ligne pour ma maladie simple
    # dans une liste de distances
  }
  return(list(
    permutations = permutation_list,
    distances = distances_list
  ))
}

test <- generate_permutations_and_distances(or_df0[6127,], or_df0,n_perm = 1)

n_perm <- 1
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
  permutation_list[[2]] <- new_df   # ajoute ma ligne permutée en dernier indice de ma liste de permutation
  temp_df <- rbind(or_df0[1:6125,], new_df) # créé temporairement une matrice les 6125 maladies 
  # simples et une maladie complexe 
  distances_sorensen <- as.matrix(dist.binary(temp_df,method = 5))[6126,1:6125]  # créé une matrice distance 
  distances_list[[2]] <- distances_sorensen # je stocke la matrice de distances à 1 ligne pour ma maladie simple
  # dans une liste de distances
}


test <- as.data.frame(distances_list[[1]])
distances_sorensen <- t(as.data.frame(distances_sorensen))

plot(density(as.matrix(o_sorensen[2,])))
lines(density(distances_sorensen),col='red')


permuted_row <- sample(row, length(row), replace = FALSE)
new_df <- permuted_row # créé matrice 1 seule ligne avec ligne permutée 
colnames(new_df) <- original_colnames
rownames(new_df) <- original_rowname
permutation_list[[2]] <- new_df   # ajoute ma ligne permutée en dernier indice de ma liste de permutation
temp_df <- rbind(or_df0[1:6125,], new_df) # créé temporairement une matrice les 6125 maladies 
# simples et une maladie complexe 
distances <- as.matrix(dist.binary(temp_df[1:2,],method = 5))




distances <- apply(temp_df, 1, FUN = function(row){
  as.matrix(dist.binary(temp_df[c(row,6126),],method = 5))[2,1]
  }
  )

