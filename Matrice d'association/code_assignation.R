# rm(list = ls())

# 03/12/2024
# Leslie Cieters
# Assignation des maladies 
# -----------------------------------------------------------------------------


library(tidyverse)
library(reshape)


### On travaille pour le moment sur la matrice originel avec la distance de 
### Sorensen

load('distances.RData')

matrice_sorensen <- as.data.frame(as.matrix(origin_sorensen)) #matrice complète


### On récupère la partie d'intérêt dans la matrice complète

o_sorensen <- matrice_sorensen[6126:7089,1:6125] #dataframe de la matrice des distances

o_sorensen_mat <- as.matrix(o_sorensen) #matrice des distances


# Création de fonctions

#' AssignGene
#' @description
#' A short description...
#'
#' @param dist 
#' @param method 
#' @param s 
#' @param q
#' @param graph 
#'
#' @return
#' @export
#'
#' @examples
AssignGene <- function(dist, method='seuil' , s = 0.5 , q = 0.25 , graph = TRUE){
  # Conversion de la matrice en format long
  flong <- melt(dist)
  colnames(flong) <- c("mc", "ms", "distance")
  
  if(method=='seuil'){
    flong_filtre <- flong %>% 
      filter(distance<=s) %>% 
      group_by(mc) %>% 
      summarise( mc = mc, ms=ms, distance=distance, liste_genes = paste(ms,collapse = " ; "),.groups = 'drop')
  } else if(method=='quantile'){
    flong_filtre <- flong %>% 
      filter(distance<=quantile(distance,q))%>% 
      group_by(mc) %>% 
      summarise( mc = mc, ms=ms, distance=distance, liste_genes = paste(ms,collapse = " ; "),.groups = 'drop')
  } else{
    stop('méthode non reconnue')
  }
  assign_liste <- unique(flong_filtre[,c(1,4)])
  assign_matrice <-  flong_filtre[,1:3] %>% 
    mutate(distance=ifelse(distance!=0,1,0)) %>% 
    pivot_wider(names_from = ms,values_from = distance) 
  
  if(graph==TRUE){
    plot <- flong_filtre %>% 
      mutate(distance=ifelse(distance!=0,1,0)) %>% 
      ggplot( aes(x = ms, y = mc, fill = distance)) +
        geom_tile() +
        scale_fill_gradient(low = "pink2", high = "firebrick") +
        labs(title = "Matrice d'Association", x = "MS", y = "MC", fill = "Distance")+
        theme(axis.text.x = element_text(angle = 30,size = 2),
            axis.text.y = element_text(size=3))
    print(plot)
  }
  
  return(list(assign_liste,assign_matrice))
  
}

test <- AssignGene(o_sorensen_mat)

test1 <- as.data.frame(test[[1]])
test2 <- as.data.frame(test[[2]])





















