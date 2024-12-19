# ------------------------------------------------------------------------------
# 18/12/2024
# Test hiérarchie des phénotypes
# ------------------------------------------------------------------------------

rm(list = ls())


# Load necessary library
if (!require("ontologyIndex")) install.packages("ontologyIndex", repos = "http://cran.us.r-project.org")
library(ontologyIndex)

# Load the HPO file
hpo_url <- "https://raw.githubusercontent.com/obophenotype/human-phenotype-ontology/master/hp.obo"
hpo <- get_ontology(hpo_url, extract_tags = 'everything')

parents <- hpo$parents



# Fonction récursive pour trouver les chemins vers la racine
find_paths_to_root <- function(term, parents, path = NULL) {
  path <- c(term, path)
  if (length(parents[[term]]) == 0) {
    return(list(path)) # At root, return the path
  }
  paths <- list()
  for (parent in parents[[term]]) {
    paths <- c(paths, find_paths_to_root(parent, parents, path))
  }
  return(paths)
}

# Récupérer tous les chemins
all_paths <- list()
for (term in names(parents)) {
  paths <- find_paths_to_root(term, parents)
  for (path in paths) {
    all_paths <- c(all_paths, list(path))
  }
}

# Trouver la profondeur maximale pour construire les colonnes
max_depth <- max(sapply(all_paths, length))

# Construire un dataframe avec des colonnes pour chaque niveau
all_paths_named <- lapply(all_paths, function(path) {
  c(path, rep(NA, max_depth - length(path))) # Compléter les chemins plus courts avec NA
})
df <- do.call(rbind, all_paths_named)
colnames(df) <- paste0("Level_", seq_len(max_depth))

df <- as.data.frame(df)
df <- df %>% filter(Level_2 == 'HP:0000118')


# Retirer les lignes ou la hiérarchie n'est pas complète 
new_df <- df
for (i in 1:dim(new_df)[1]) {
  for (j in 1:(dim(new_df)[2]-1)) {
    if (length(children[[new_df[dim(new_df)[1]+1-i,j]]])!=0 & is.na(new_df[dim(new_df)[1]+1-i,j+1])) {
      new_df <- new_df[-(dim(new_df)[1]+1-i),]
    }
    
  }
  
}

#  comparaison des colonnes

# Fonction pour comparer deux colonnes avec correspondances exactes
compare_columns_exact <- function(col1, col2) {
  common <- intersect(col1, col2) # Correspondances exactes
  return(data.frame(value = common))
}

# Comparer toutes les colonnes deux à deux dans le dataframe
compare_all_columns <- function(df) {
  result <- list()
  col_names <- colnames(df)
  
  for (i in 1:(ncol(df) - 1)) {
    for (j in (i + 1):ncol(df)) {
      col1 <- df[[i]]
      col2 <- df[[j]]
      result[[paste0(col_names[i], "_vs_", col_names[j])]] <- compare_columns_exact(col1, col2)
    }
    }
  return(result)
}

# Comparaison avec correspondances exactes
exact_matches <- compare_all_columns(new_df)
print("Correspondances exactes :")
print(exact_matches)


new_df <- new_df %>% 
  mutate_if(is.character, as.factor)

sum <- 0
for (i in 1:136) {
  sum <- sum+length(exact_matches[[i]][[1]])
}


df2 <- as.data.frame(sapply(new_df,gsub,pattern=":",replacement="."))

# Filtrer les lignes avec au moins un match
filtered_df <- df2[apply(df2, 1, function(row) any(row %in% colnames(or_df0))), ]



# Fonction pour trouver la dernière valeur non-NA dans une ligne
last_non_na <- function(row) {
  row <- as.character(row) # S'assurer que les valeurs sont des chaînes
  row <- row[!is.na(row)]  # Filtrer les valeurs non-NA
  if (length(row) > 0) {
    return(tail(row, 1))  # Retourner la dernière valeur
  } else {
    return(NA)  # Retourner NA si la ligne est entièrement NA
  }
}

# Filtrer les lignes où la dernière valeur non-NA est dans les noms des colonnes de or_df0
filtered_df <- df2[apply(df2, 1, function(row) last_non_na(row) %in% colnames(or_df0)), ]

# Comparaison avec correspondances exactes
exact_matches_filtered <- compare_all_columns(filtered_df)



