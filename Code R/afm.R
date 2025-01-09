actif <- mandale$ind$coord
row_acm <- mandale$ind.sup$coord
dta_coord <- rbind(actif, sup)

save(dta_coord, file = "dta_coord_acm.RData")
head(dta_coord)

# Pour Sorensen
row_sorensen = result_sorensen$conf.row
row_ochiai = result_ochiai$conf.row
row_avant_OT = result_avant_OT$conf.row
row_apres_OT = result_apres_OT$conf.row

data_pour_afm <- cbind(row_sorensen, row_ochiai, row_avant_OT, row_apres_OT, row_acm)
data_pour_afm <- as.data.frame(data_pour_afm)

library(FactoMineR)
res.mfa <- MFA(data_pour_afm, 
               group = c(5, 5, 5, 5, 5),
               type = c("c", "c", "c" , "c", "c"),  # toutes les variables sont continues
               ncp = 2,             # nombre de dimensions à retenir
               name.group = c("Sorensen", "Ochiai", "Avant_OT", "Apres_OT", "ACM"))

rm(dist_or_ochiai_mx)
rm(dist_or_sorensen_mx)
dist_or_ochiai_mx <- as.data.frame(dist_or_ochiai_mx)
dist_or_sorensen_mx <- as.data.frame(dist_or_sorensen_mx)
