
# Application de l'algorithme de permutation sur l'ensemble des maladies complexes
# Uniquement avec la distance de Sorensen pour le moment 




#------------------------- Seuil -----------------------------------------------

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

# Optimized permutation analysis function
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

# Example usage:
results <- process_multiple_rows(or_df0, n_perm = 1, start_idx = 6126, end_idx = 7091)

#save(results, file = "permuations_10_all.RData")

#------------------------ Nb MS ----------------------------------------------
#Ajout temporaire noms
or_df0$rowname <- rownames(or_df0)

# Fusionner sur la colonne "rowname"
merged_results <- merge(or_df0, results, by = "rowname", all.x = TRUE)

# Optionnel : Rétablir les noms de lignes dans le résultat final
rownames(merged_results) <- merged_results$rowname
merged_results$rowname <- NULL

# Remplacer les NA par 0 dans les colonnes venant de results 
merged_results[is.na(merged_results)] <- 0

merged_results$Count <- apply(merged_results, 1, function(row) sum(as.numeric(row), na.rm = TRUE))
df_plot <- merged_results[,5384:5385]
