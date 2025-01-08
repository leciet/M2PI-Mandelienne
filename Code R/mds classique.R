rm(list=ls())

actif <- mandale$ind$coord
sup <- mandale$ind.sup$coord
dta_dist <- rbind(actif, sup)

library(factoextra)
dist_acm_eucli <- get_dist(dta_dist, method = "euclidean")
dist_acm_eucli <- as.data.frame(as.matrix(dist_acm_eucli))

dist_acm_manhattan<- get_dist(dta_dist, method = "manhattan")
dist_acm_manhattan <- as.data.frame(as.matrix(dist_acm_manhattan))

or_filter_df0 <- or_df0 %>% 
  filter(rowSums(or_df0)!=0)

dist_or_sorensen <- dist.binary(or_filter_df0,method = 5) #Long 
dist_or_sorensen <- as.data.frame(as.matrix(dist_or_sorensen))

dist_or_ochiai<- dist.binary(or_filter_df0,method = 7)
dist_or_ochiai <- as.data.frame(as.matrix(dist_or_ochiai))

dist_or_sorensen_mx<- as.matrix(dist_or_sorensen)
dist_or_ochiai_mx  <- as.matrix(dist_or_ochiai)
dist_acm_manhattan_mx <- as.matrix(dist_acm_manhattan)
dist_acm_eucli_mx <- as.matrix(dist_acm_eucli)

library(smacof)
data_pour_mds <- cbind(dist_or_sorensen_mx,
                       dist_or_ochiai_mx,
                       dist_acm_eucli_mx,
                       dist_acm_manhattan_mx
)

# Si vous avez plusieurs matrices de distance
list_matrices <- list(
  matrice1 = as.dist(dist_or_sorensen_mx),
  matrice2 = as.dist(dist_or_ochiai_mx)
#  matrice3 = as.dist(dist_acm_eucli_mx),
#  matrice4 = as.dist(dist_acm_manhattan_mx)
)

save(dist_or_sorensen_mx, file = "dist_or_sorensen_mx.RData")
save(dist_or_ochiai_mx, file = "dist_or_ochiai_mx.RData")
save(dist_acm_eucli_mx, file = "dist_acm_eucli_mx.RData")
save(dist_acm_manhattan_mx, file = "dist_acm_manhattan_mx.RData")

# Pour n matrices de distance dans une liste
mds_indscal <- smacofIndDiff(delta = list_matrices,  
                             ndim = 2,       # dimensions
                             type = "ratio", 
                             itmax = 15,
                             verbose = TRUE)

# Visualisation
plot(mds_indscal)

mds_result <- smacofSym(as.dist(dist_or_sorensen_mx), 
                        ndim = 100000,
                        type = "ratio",
                        itmax = 1000,
                        verbose = TRUE)

mds_result_acm_eucli <- smacofSym(as.dist(dist_acm_eucli_mx), 
                                  ndim = 100000,
                                  type = "ratio",
                                  itmax = 1000,
                                  verbose = TRUE)

mds_result_acm_manhattan <- smacofSym(as.dist(dist_acm_manhattan_mx), 
                                      ndim = 100000,
                                      type = "ratio",
                                      itmax = 1000,
                                      verbose = TRUE)

mds_result_avant_OT <- smacofRect(delta = matrice_distance_Jules_filtered,  # votre matrice
                                  ndim = 100000,        # même nombre de dimensions
                                  type = "ratio",   # même type d'échelle
                                  itmax = 250,    # nombre max d'itérations
                                  eps = 1e-6,      # critère de convergence
                                  verbose = TRUE)   # affichage des étapes

mds_result_apres_OT <- smacofRect(delta = matrice_distance_apres_OT,  # votre matrice
                                  ndim = 100000,        # même nombre de dimensions
                                  type = "ratio",   # même type d'échelle
                                  itmax = 250,    # nombre max d'itérations
                                  eps = 1e-6,      # critère de convergence
                                  verbose = TRUE)   # affichage des étapes

mds_result_ochiai <- smacofSym(as.dist(dist_or_ochiai_mx), 
                               ndim = 2,
                               type = "ratio",
                               itmax = 1000,
                               verbose = TRUE)

