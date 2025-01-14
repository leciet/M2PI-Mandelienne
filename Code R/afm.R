col_acm <- mandale$ind$coord
row_acm <- mandale$ind.sup$coord
row_acm <- as.data.frame(row_acm)
dta_coord <- rbind(actif, sup)

save(dta_coord, file = "dta_coord_acm.RData")
head(dta_coord)

# Pour Sorensen
row_sorensen = result_sorensen$conf.row
row_sorensen <- as.data.frame(row_sorensen)

row_ochiai = result_ochiai$conf.row
row_ochiai <- as.data.frame(row_ochiai)

row_avant_OT = result_avant_OT$conf.row
row_avant_OT <- as.data.frame(row_avant_OT)

row_apres_OT = result_apres_OT$conf.row
row_apres_OT <- as.data.frame(row_apres_OT)

data_pour_afm <- cbind(row_sorensen, row_ochiai, row_avant_OT, row_apres_OT, row_acm)
data_pour_afm <- as.data.frame(data_pour_afm)
save(data_pour_afm, file = "data_pour_afm.RData")
head(data_pour_afm)
library(FactoMineR)


dac

res.mfa <- MFA(data_pour_afm, 
               group = c(5, 5, 5, 5, 5),
               type = c("c", "c", "c" , "c", "c"),  # toutes les variables sont continues
               ncp = 2,             # nombre de dimensions à retenir
               name.group = c("Sorensen", "Ochiai", "Jaccard", "Jacccard_OT", "ACM"))

res.mfa$summary.quanti
res.mfa$summary.quali
res.mfa$eig
res.mfa$group
res.mfa$inertia.ratio
res.mfa$ind$coord
res.mfa$quanti.var
res.mfa$partial.axes
res.mfa$call

# mtn let's go mettre mes maladies simples et mes maladies complexes dans le même plan

col_acm <- mandale$ind$coord
col_acm <- as.data.frame(col_acm)

col_sorensen = result_sorensen$conf.col
col_sorensen <- as.data.frame(col_sorensen)

col_ochiai = result_ochiai$conf.col
col_ochiai  <- as.data.frame(col_ochiai )

col_avant_OT = result_avant_OT$conf.col
col_avant_OT <- as.data.frame(col_avant_OT)

col_apres_OT = result_apres_OT$conf.col
col_apres_OT <- as.data.frame(col_apres_OT)

names(row_sorensen)
names(row_ochiai)
names(row_avant_OT)
names(row_apres_OT)
names(row_acm)
names(col_acm)
names(col_sorensen)
names(col_ochiai)
names(col_avant_OT)
names(col_apres_OT)

# Rename the columns of row_acm and col_acm to match the others
names(row_acm) <- c("D1", "D2", "D3", "D4", "D5")
names(col_acm) <- c("D1", "D2", "D3", "D4", "D5")

# Now try the rbind again
data_pour_afm_supp <- rbind(row_sorensen, 
                            row_ochiai, 
                            row_avant_OT, 
                            row_apres_OT, 
                            row_acm,
                            col_acm, 
                            col_sorensen,
                            col_ochiai, 
                            col_avant_OT, 
                            col_apres_OT)

data_pour_afm_supp <- as.data.frame(data_pour_afm_supp)
# data_pour_afm_supp_2 <- t(data_pour_afm_supp)
# data_pour_afm_supp_2 <- as.data.frame(data_pour_afm_supp_2)
# head(data_pour_afm_supp)

res.pca <- PCA(data_pour_afm_supp, 
               scale.unit = TRUE,  # Important pour standardiser
               ncp = 5)           # Nombre de dimensions à conserver

# Créer un vecteur de couleurs
colors <- c(rep("red", 4820), rep("blue", 30625))
plot(res.pca, choix = "ind", col.ind = colors)

# CAH avec la méthode Ward
res.hcpc <- HCPC(res.pca,             # Résultat de l'ACP
                 nb.clust = -1,        # Nombre de clusters (-1 = automatique)
                 consol = TRUE)        # Consolidation
plot(res.hcpc, choice = "3D.map")     # Carte 3D des clusters
plot(res.hcpc, choice = "tree")       # Dendrogramme


# K-means sur les coordonnées
coord_ind <- res.pca$ind$coord
res.kmeans <- kmeans(coord_ind,        # Coordonnées des individus
                     centers = 10,       # Nombre de clusters souhaité
                     nstart = 25)       # Nombre de départs aléatoires
data_pour_afm_supp$cluster <- as.factor(res.kmeans$cluster) 

library(cluster)
sil <- silhouette(res.kmeans$cluster, dist(coord_ind))
plot(sil)  
rm(sil)
rm(plot_data)
rm(res.kmeans)

save(data_pour_afm, file = "data_pour_afm.RData")
save(data_pour_afm_supp, file = "data_pour_afm_supp.RData")
rm(list=ls())
data_pour_afm_supp <- data_pour_afm_supp[,-6]
memory.size()
gc()
memory.size()
memory.limit(size = 16000)

library(ggplot2)
plot_data <- data.frame(
  PC1 = coord_ind[,1],
  PC2 = coord_ind[,2],
  Cluster = as.factor(res.kmeans$cluster)
)
ggplot(plot_data, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(title = "K-means Clustering in PCA Space")
