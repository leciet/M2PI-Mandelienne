library(tidyverse)

analyser_patterns_maladies <- function(data) {
  # Convertir en matrice si ce n'est pas déjà le cas
  data <- as.matrix(data)
  
  # Créer un identifiant unique pour chaque pattern de 1
  patterns <- apply(data, 1, paste, collapse = "")
  
  # Identifier les patterns qui apparaissent plus d'une fois
  patterns_multiples <- patterns[duplicated(patterns) | duplicated(patterns, fromLast = TRUE)]
  patterns_uniques <- unique(patterns_multiples)
  
  # Créer le dataframe pour les patterns qui se répètent
  resultats <- data.frame(
    Pattern_ID = integer(0),
    Maladies = character(0),
    Nombre_1 = integer(0),
    Positions_1 = character(0),
    Nombre_Apparitions = integer(0),
    stringsAsFactors = FALSE
  )
  
  # Analyser chaque pattern qui se répète
  for (i in seq_along(patterns_uniques)) {
    pattern <- patterns_uniques[i]
    maladies_concernees <- rownames(data)[patterns == pattern]
    nombre_1 <- sum(data[maladies_concernees[1],] == 1)
    positions_1 <- paste(which(data[maladies_concernees[1],] == 1), collapse = ", ")
    
    resultats <- rbind(resultats, data.frame(
      Pattern_ID = i,
      Maladies = paste(maladies_concernees, collapse = ", "),
      Nombre_1 = nombre_1,
      Positions_1 = positions_1,
      Nombre_Apparitions = length(maladies_concernees)
    ))
  }
  
  return(resultats)
}

