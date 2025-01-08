library(MASS)
library(smacof)

mds_result_sorensen <- smacofRect(delta = dist_or_sorensen_mx,  # votre matrice
                                  ndim = 100000,        # même nombre de dimensions
                                  type = "ratio",   # même type d'échelle
                                  itmax = 250,    # nombre max d'itérations
                                  eps = 1e-6,      # critère de convergence
                                  verbose = TRUE)   # affichage des étapes

mds_result_ochiai <- smacofRect(delta = dist_or_ochiai_mx,  # votre matrice
                                  ndim = 100000,        # même nombre de dimensions
                                  type = "ratio",   # même type d'échelle
                                  itmax = 250,    # nombre max d'itérations
                                  eps = 1e-6,      # critère de convergence
                                  verbose = TRUE)   # affichage des étapes

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

mds_result_acm_eucli <- smacofRect(delta = dist_acm_eucli_mx,  # votre matrice
                                  ndim = 100000,        # même nombre de dimensions
                                  type = "ratio",   # même type d'échelle
                                  itmax = 250,    # nombre max d'itérations
                                  eps = 1e-6,      # critère de convergence
                                  verbose = TRUE)   # affichage des étapes

mds_result_acm_manhattan <- smacofRect(delta = dist_acm_manhattan_mx,  # votre matrice
                                  ndim = 100000,        # même nombre de dimensions
                                  type = "ratio",   # même type d'échelle
                                  itmax = 250,    # nombre max d'itérations
                                  eps = 1e-6,      # critère de convergence
                                  verbose = TRUE)   # affichage des étapes

mds_result_sorensen <- unfolding(dist_or_sorensen_mx, 
                        ndim = 2,        # nombre de dimensions
                        type = "ratio",   # type d'échelle
                        verbose = TRUE)  # pour éviter trop de messages

mds_result_ochiai <- unfolding(dist_or_ochiai_mx, 
                                 ndim = 2,        # nombre de dimensions
                                 type = "ratio",   # type d'échelle
                                 verbose = TRUE)  # pour éviter trop de messages

mds_result_avant_OT <- unfolding(matrice_distance_Jules_filtered, 
                               ndim = 2,        # nombre de dimensions
                               type = "ratio",   # type d'échelle
                               verbose = TRUE)  # pour éviter trop de messages

mds_result_apres_OT <- unfolding(matrice_distance_apres_OT, 
                               ndim = 2,        # nombre de dimensions
                               type = "ratio",   # type d'échelle
                               verbose = TRUE)  # pour éviter trop de messages

mds_result_acm_eucli <- unfolding(dist_acm_eucli_mx, 
                                 ndim = 2,        # nombre de dimensions
                                 type = "ratio",   # type d'échelle
                                 verbose = TRUE)  # pour éviter trop de messages

mds_result_acm_manhattan <- unfolding(dist_acm_manhattan_mx, 
                                 ndim = 2,        # nombre de dimensions
                                 type = "ratio",   # type d'échelle
                                 verbose = TRUE)  # pour éviter trop de messages