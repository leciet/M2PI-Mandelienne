# ---------- 
# Étude préliminaire du jeu de données sans hiérarchie
# 15/01/2025
# Leslie Cieters
# ----------

# rm(list = ls())
# Regarder la structure des données sans hiérarchie

load("Matrices Jules sans hiérarchie/phenotype_maladie_s_c (1).RData")



min(rowSums(phenotype_maladie_s_c)) #il y a minimum 1 phénotype par maladie 
max(rowSums(phenotype_maladie_s_c)) #il y a maximum 203 phénotype par maladie

min(colSums(phenotype_maladie_s_c))
max(colSums(phenotype_maladie_s_c))


# on a 6102 maladies simples et 962 maladies complexes  
min(colSums(phenotype_maladie_s_c[1:6102,]))
max(colSums(phenotype_maladie_s_c[6103:7064,]))





