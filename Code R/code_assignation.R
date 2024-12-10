# rm(list = ls())

# 03/12/2024
# Leslie Cieters
# Assignation des maladies 
# -----------------------------------------------------------------------------


library(tidyverse)
library(reshape)




### On travaille pour le moment sur la matrice originel avec la distance de 
### Sorensen


load("Data/distance_origin_ochiai.RData")
load("Data/distance_origin_sorensen.RData")


# Mise en forme
dist_or_ochiai_df <- as.data.frame(dist_or_ochiai_mx)

# Appliquer le test de Shapiro-Wilk à chaque ligne
norm <- apply(dist_or_sorensen_mx, 1, function(row) {
  ks.test(row, "pnorm", mean = mean(row), sd = sd(row))$p.value
})
norm <- as.data.frame(norm)


# On n'a pas de distribution normale 
# Pour passer outre on peut s'intérésser à la p value empirique

compute_empirical_pvalues <- function(row) {
  sapply(row, function(d) {
    sum(row <= d) / length(row)  # Proportion de distances <= d
  })
}

# Calculer les p-values empiriques pour chaque ligne
empirical_pvalues <- t(apply(dist_or_sorensen_mx, 1, compute_empirical_pvalues))


# Sélection des gènes avec un seuil de p-value
alpha <- 0.005
selected_genes <- data.frame(ifelse(empirical_pvalues<=alpha,1,0))

selected_genes_long <- melt(as.matrix(selected_genes))
colnames(selected_genes_long) <- c("mc", "ms", "assignation")

selected_genes_long <- selected_genes_long %>% 
  filter(assignation == 1) %>% 
  group_by(mc) %>% 
  summarise( mc = mc, liste_genes = paste(ms,collapse = " ; "),.groups = 'drop')

selected_genes_liste <- as.data.frame(unique(selected_genes_long))


print("Empirical p-values:")
print(empirical_pvalues)
print("Selected genes:")
print(selected_genes)






# Création de fonctions

#' AssignGene
#' @description
#' Assigne une liste de gènes à chaque maladie complexe, la matrice de distances a en colonne les ms et en ligne les mc
#'
#' @param dist objet  matrix correspondant à la matrice de distances
#' @param method 'seuil' ou 'quantile' pour le critère d'association
#' @param s paramètre pour la méthode seuil
#' @param q paramètre pour la méthode quantile
#' @param graph si TRUE affiche la matrice d'association
#'
#' @return une liste contenant 2 objets correspondant à l'assignation d'une liste de gènes par maladie et la matrice d'assignation
#' @export
#'
#' @examples
AssignGene <- function(dist, method='seuil' , s = 0.5 , q = 0.25 , graph = TRUE){
  # Conversion de la matrice en format long
  flong <- melt(dist)
  colnames(flong) <- c("mc", "ms", "distance")
  
  if(method=='seuil'){
    if (s<0 | s>1) {
      stop("Le seuil doit être compris entre 0 et 1")
    }
    flong_filtre <- flong %>% 
      filter(distance<=s) %>% 
      group_by(mc) %>% 
      summarise( mc = mc, ms=ms, distance=distance, liste_genes = paste(ms,collapse = " ; "),.groups = 'drop')
  } else if(method=='quantile'){
    if (q<0 | q>1) {
      stop("Le quantile doit être compris entre 0 et 1")
    }
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


### On récupère la partie d'intérêt dans la matrice complète


matrice_sorensen <- as.data.frame(as.matrix(origin_sorensen)) #matrice complète

o_sorensen <- matrice_sorensen[6126:7089,1:6125] #dataframe de la matrice des distances

o_sorensen_mat <- as.matrix(o_sorensen) #matrice des distances

test <- AssignGene(dist_or_sorensen_mx,method = 'seuil',s = 0.84)

test1 <- as.data.frame(test[[1]])
test2 <- as.data.frame(test[[2]])



matrice_ochiai <- as.data.frame(as.matrix(origin_ochiai)) #matrice complète

o_ochiai <- matrice_ochiai[6126:7089,1:6125] #dataframe de la matrice des distances

o_ochiai_mat <- as.matrix(o_ochiai) #matrice des distances

asso_origin_ochiai <- AssignGene(o_ochiai_mat,method = 'seuil',s = 0.84)[[2]]




asso_origin_ochiai <- o_ochiai
asso_origin_ochiai <- data.frame(ifelse(asso_origin_ochiai<=0.84,1,0))

asso_origin_sorensen <- o_sorensen
asso_origin_sorensen <- data.frame(ifelse(dist_or_sorensen_mx<=0.84,1,0))

save(asso_origin_ochiai,asso_origin_sorensen,file="Matrice d'association/asso_origin.RData")








