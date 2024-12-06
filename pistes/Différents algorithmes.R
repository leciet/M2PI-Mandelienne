# install.packages(IBFS)  pas sûr
# install.packages(rTephra) pas sûr
# install.packages(ReacTran) pas sûr
# install.packages(POT)

install.packages("transport")

install.packages("causalOT")
install.packages("fdWasserstein")
install.packages("T4transport")
install.packages("gridOT")
install.packages("OTrecod")
install.packages("otinference")
install.packages("TransP")
install.packages("lpSolve")
install.packages()


# Chargement des packages avec bcp de biblio
library(transport)
library(fdWasserstein)
library(T4transport)
library(OTrecod)
library(lpSolve)

# Chargement des packages avec peu de biblio
library(causalOT)
library(gridOT)

# Chargement des packages avec pas de biblio
library(TransP)
library(otinference)
#library(POT)

# Création des données de test
set.seed(123)
n <- 5
# Matrice de coût aléatoire
cost_matrix <- matrix(runif(n*n), nrow=n, ncol=n)
# Distributions de masse (uniformes pour l'exemple)
a <- rep(1/n, n)
b <- rep(1/n, n)

# 1. Package 'transport' : donne des matrices de similarité
# ----------------------
# Méthode du simplexe (version réseau)
network_result <- transport::transport(
  a, b, 
  costm = cost_matrix, 
  method = "networkflow"
)

# Méthode du simplexe court
shortsimplex_result <- transport::transport(
  a, b, 
  costm = cost_matrix, 
  method = "shortsimplex"
)

# Méthode du simplexe révisé
revsimplex_result <- transport::transport(
  a, b, 
  costm = cost_matrix, 
  method = "revsimplex"
)

# Méthode primale-duale
primaldual_result <- transport::transport(
  a, b, 
  costm = cost_matrix, 
  method = "primaldual"
)

# 2. Package causalOT
# ----------------------

# Méthode de sinkhorn

# Calcul de Sinkhorn avec régularisation

x2 <- cost_matrix

sinkhorn_result <- ot_distance(
  x1 = cost_matrix,  # Matrice des coûts
  x2 = x2,           # Matrice des covariables cible
  a = a,             # Poids source
  b = b,             # Poids cible
  penalty = 100,     # Paramètre de régularisation (épsilon)
  p = 2,             # Distance L2
  debias = TRUE,     # Option pour Sinkhorn débiassé
  niter = 25,      # Nombre d'itérations
  tol = 1e-07        # Tolérance de convergence
)

print("Divergences de Sinkhorn:")
print(sinkhorn_result)
summary(sinkhorn_result)


# Méthode du point intérieur
result_point_interieur <- causalOT(a, b, cost_matrix, method = "IPFP")

print("Matrice de transport optimale (Point Intérieur - causalOT):")
print(result_point_interieur$plan)





# Méthode des enchères avec T4transport
auction_result <- T4transport::auction(
  costm = cost_matrix,
  mass.1 = a,
  mass.2 = b
)

# Méthode stochastique avec causalOT
stochastic_result <- causalOT::wasserstein(
  x = matrix(1:n, ncol = 1),
  y = matrix(1:n, ncol = 1),
  w_x = a,
  w_y = b,
  p = 2,
  method = sgd,            # Méthode stochastique
  num_iterations = 1000
)

# Méthode multi-échelle avec gridOT
# Création d'une grille pour l'exemple
grid_x <- seq(0, 1, length.out = n)
grid_y <- seq(0, 1, length.out = n)
multi_scale_result <- gridOT::gridOT(
  source = matrix(a, nrow = n),
  target = matrix(b, nrow = n),
  grid_x = grid_x,
  grid_y = grid_y
)

# Affichage des résultats
cat(\nRésultats Sinkhorn:\n)
print(sinkhorn_result)

cat(\nRésultats Auction:\n)
print(auction_result)

cat(\nRésultats Stochastic:\n)
print(stochastic_result)

cat(\nRésultats Multi-échelle:\n)
print(multi_scale_result)

# Comparaison des temps d'exécution
library(microbenchmark)

benchmark_results <- microbenchmark(
  Networkflow = transport::transport(a, b, costm = cost_matrix, method = networkflow),
  Shortsimplex = transport::transport(a, b, costm = cost_matrix, method = shortsimplex),
  Revsimplex = transport::transport(a, b, costm = cost_matrix, method = revsimplex),
  Primaldual = transport::transport(a, b, costm = cost_matrix, method = primaldual),
  Sinkhorn = causalOT::wasserstein(matrix(1:n, ncol = 1), matrix(1:n, ncol = 1), 
                                   a, b, p = 2, method = sinkhorn, epsilon = 0.01),
  Auction = T4transport::auction(costm = cost_matrix, mass.1 = a, mass.2 = b),
  Stochastic = causalOT::wasserstein(matrix(1:n, ncol = 1), matrix(1:n, ncol = 1), 
                                     a, b, p = 2, method = sgd, num_iterations = 1000),
  times = 5
)

print(benchmark_results)

# Vérification de la conservation des masses pour chaque méthode
check_mass_conservation <- function(transport_plan, a, b) {
  row_sums <- rowSums(transport_plan)
  col_sums <- colSums(transport_plan)
  max_error_row <- max(abs(row_sums - a))
  max_error_col <- max(abs(col_sums - b))
  return(list(row_error = max_error_row, col_error = max_error_col))
}




