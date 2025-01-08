rm(dist_acm_eucli_mx)
rm(dist_acm_manhattan_mx)
rm(dist_or_ochiai_mx)
rm(dist_or_sorensen_mx)
rm(or_df0)
rm(phecode)
rm(mandale)
rm(matrice_distance)
rm(matrice_distance_ou_similarite)
rm(nom_col)
rm(colonnes_manquantes)

matrice_distance_Jules <- read.csv("distance_matrix_dist_Jules.csv")
dist_or_sorensen_mx <- as.data.frame(dist_or_sorensen_mx)
dist_or_ochiai_mx <- as.data.frame(dist_or_ochiai_mx )
dist_acm_eucli_mx <- as.data.frame(dist_acm_eucli_mx)
dist_acm_manhattan_mx <- as.data.frame(dist_acm_manhattan_mx )

# nom_ligne
# nom_ligne <- rownames(dist_acm_manhattan_mx)
# nom_ligne
# 
# nom_col <- colnames(matrice_distance_Jules)
# nom_col <- gsub("\\.", " ", nom_col)
# nom_col

library(dplyr)
matrice_distance_Jules_filtered <- matrice_distance_Jules %>%
  select(-c(Acute.bronchospasm, 
            Nonspecific.abnormal.findings.on.radiological.and.other.examination.of.other.intrathoracic.organs..echocardiogram..etc.))
matrice_distance_Jules_filtered <- t(matrice_distance_Jules_filtered)
matrice_distance_Jules_filtered  <- as.data.frame(matrice_distance_Jules_filtered )
colnames(matrice_distance_Jules_filtered) <- matrice_distance_Jules_filtered[1, ]
matrice_distance_Jules_filtered <- matrice_distance_Jules_filtered[-1, ]
# View(matrice_distance_Jules_filtered)

matrice_similarite <- read.csv("agro_to.csv")
matrice_similarite_t <- matrice_similarite %>%
  select(-c(Acute.bronchospasm, 
            Nonspecific.abnormal.findings.on.radiological.and.other.examination.of.other.intrathoracic.organs..echocardiogram..etc.))
matrice_similarite_t <- t(matrice_similarite_t)
matrice_similarite_t  <- as.data.frame(matrice_similarite_t)
colnames(matrice_similarite_t ) <- matrice_similarite_t [1, ]
matrice_similarite_t  <- matrice_similarite_t[-1, ]

matrice_similarite_t_numeric <- matrice_similarite_t %>%
  mutate(across(everything(), ~ as.numeric(as.character(.))))
matrice_distance_apres_OT <- as.data.frame(1 - matrice_similarite_t_numeric)
# View(matrice_distance_apres_OT)


max(dist_or_sorensen_mx)
min(dist_or_sorensen_mx)
max(dist_or_ochiai_mx)
min(dist_or_ochiai_mx)
max(dist_acm_eucli_mx)
min(dist_acm_eucli_mx)
max(dist_acm_manhattan_mx)
min(dist_acm_manhattan_mx)
matrice_distance_Jules_filtered <- matrice_distance_Jules_filtered %>%
  mutate(across(everything(), ~ as.numeric(as.character(.))))
max(matrice_distance_Jules_filtered)
min(matrice_distance_Jules_filtered)
max(matrice_distance_apres_OT)
min(matrice_distance_apres_OT)

save(matrice_distance_Jules_filtered, file = "matrice_avant_OT.RData")
save(matrice_distance_apres_OT, file = "matrice_apres_OT.RData")

data_pour_afc <- cbind(dist_or_sorensen_mx,
                       dist_or_ochiai_mx,
                       dist_acm_eucli_mx,
                       dist_acm_manhattan_mx,
                       matrice_distance_Jules_filtered,
                       matrice_distance_apres_OT
)

data_pour_afc <- as.data.frame(lapply(data_pour_afc, as.numeric))

summary(data_pour_afc)
str(data_pour_afc)
max(data_pour_afc)
min(data_pour_afc)

# data_pour_afc <- data.frame(data_pour_afc, check.names = FALSE)

# colnames(data_pour_afc)[1:6125] <- "Sorensen"
# colnames(data_pour_afc)[6126:12250] <- "Ochiai"
# colnames(data_pour_afc)[12251:18375] <- "Eucli"
# colnames(data_pour_afc)[18376:24500] <- "Manhattan"
# colnames(data_pour_afc)[24501:30625] <- "Avant_OT"
# colnames(data_pour_afc)[30626:36750] <- "Après_OT"
# 
# colnames(data_pour_afc) <- sub("\\..*", "", colnames(data_pour_afc))
# colnames(data_pour_afc) <- sapply(strsplit(colnames(data_pour_afc), "\\."), `[`, 1)
# 
# head(colnames(data_pour_afc))
# colnames(data_pour_afc)

data_pour_afc_simplifie <- data_pour_afc[1:10, ]
head(colnames(data_pour_afc_simplifie[24501:30625]))
data_pour_afc_simplifie <- as.data.frame(lapply(data_pour_afc_simplifie, as.numeric))

library(FactoMineR)
library(factoextra)

afc <- CA(data_pour_afc_simplifie)

res.mfa <- MFA(data_pour_afc_simplifie, 
               group = c(6125, 6125, 6125, 6125, 6125, 6125),
               type = c("c", "c", "c" , "c", "c", "c"),  # toutes les variables sont continues
               ncp = 5,             # nombre de dimensions à retenir
               name.group = c("Sorensen", "Ochiai", "Eucli", "Manhattan", "Avant_OT", "Apres_OT"))

##########################################################

data_pour_afc <- cbind(dist_or_sorensen_mx,
                       dist_or_ochiai_mx,
                       dist_acm_eucli_mx,
                       dist_acm_manhattan_mx,
                       matrice_distance_Jules_filtered,
                       matrice_distance_apres_OT
)

data_pour_afc <- as.data.frame(lapply(data_pour_afc, as.numeric))
save(data_pour_afc, file = "data_pour_afc.RData")
write.csv(data_pour_afc, "data_pour_afc.csv")

data_pour_afc_simplifie <- data_pour_afc[1:10, ]
data_pour_afc_simplifie <- as.data.frame(lapply(data_pour_afc_simplifie, as.numeric))

res.mfa <- MFA(data_pour_afc_simplifie, 
               group = c(6125, 6125, 6125, 6125, 6125, 6125),
               type = c("c", "c", "c" , "c", "c", "c"),  # toutes les variables sont continues
               ncp = 5,             # nombre de dimensions à retenir
               name.group = c("Sorensen", "Ochiai", "Eucli", "Manhattan", "Avant_OT", "Apres_OT"))

res.mfa$group$RV
# Fortes similarités (RV > 0.9)
# Similarités moyennes (RV entre 0.5 et 0.9)
# Faibles similarités (RV < 0.5)
plot(res.mfa, choix = "axes", habillage = "group")
# Comment chaque métrique de distance caractérise différemment les mêmes maladies

cor(res.mfa$separate.analyses$Sorensen$ind$coord,
    res.mfa$separate.analyses$Ochiai$ind$coord)
# Première ligne (Dim.1 de Sorensen) et colonne (Dim.1 de Ochiai)

fviz_mfa_ind(res.mfa, repel = TRUE)

fviz_mfa_var(res.mfa, repel = TRUE)

fviz_mfa_axes(res.mfa, repel = TRUE)

fviz_mfa(res.mfa, repel = TRUE)
