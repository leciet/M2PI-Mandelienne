# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 15:22:48 2024

@author: nemod
"""
#%%    
import numpy as np
#%%


# Définition des données initiales
D = np.array([
    [0.5, 2.1, 3.2],
    [1.8, 1.2, 2.5]
])
#%%

# a et b sont les distributions marginales, c'est-à-dire :
# a représente la quantité de "masse" disponible à chaque point de départ
# b représente la quantité de "masse" désirée à chaque point d'arrivée
# 2 entrepôts (sources) avec respectivement 100 et 50 unités de produits
# 3 magasins (destinations) qui ont besoin de 60, 50 et 40 unités
# Distributions marginales cibles
a = np.array([100, 50])  # somme = 150
b = np.array([60, 50, 40])  # somme = 150
# L'algorithme de Sinkhorn converge vers une solution unique qui respecte les
# contraintes marginales définies par a et b.
# C'est a et b qui définissent les contraintes marginales que le plan de
# transport final doit respecter
# C'est comme si a et b étaient les "cibles" que l'algorithme doit atteindre, 
# peu importe son point de départ.
#%%

# Paramètres
epsilon = 0.1  # paramètre de régularisation
n_iter = 5  # nombre d'itérations
#%%

# Calcul de la matrice de kernel
P = np.exp(-D/epsilon)
#%%

# Initialisation des vecteurs d'échelle
u = np.ones(len(a))
v = np.ones(len(b))
# On initialise u et v avec des 1 car :
# C'est une initialisation neutre et standard pour l'algorithme de Sinkhorn
# Ces valeurs vont être mises à jour itérativement pendant l'algorithme
# Le choix exact des valeurs initiales n'est pas crucial car l'algorithme 
# convergera vers la même solution
# Les mises à jour itératives de u et v servent justement à transformer la
#  matrice P en un plan de transport T qui respecte ces contraintes a et b
#%%

print("Matrice P (kernel) :")
print(P)
print("\nItérations de l'algorithme de Sinkhorn :")
print(f"Voici a : {a} et voici u : {u}")
print(f"Voici b : {b} et voici v : {v}")
#%%

for i in range(n_iter):
    # Mise à jour de u
    u = a / (P @ v)   # @ est le multiplicateur matricielle
    
    # Mise à jour de v
    v = b / (P.T @ u)  # transposée de la matrice P 
    
    # P1m = a
    # P⊤1n = b
    # elle satisfait les contraintes marginales
    # Calcul du plan de transport courant
    T = np.diag(u) @ P @ np.diag(v)
    
    print(f"\nItération {i+1}")
    print(f"u = {u}")
    print(f"v = {v}")
    print("Plan de transport T:")
    print(T)
    print(f"Marges lignes: {T.sum(axis=1)}")
    print(f"Marges colonnes: {T.sum(axis=0)}")
    #%%
print("P1m = ", T.sum(axis=1))  # devrait converger vers a
print("P⊤1n = ", T.sum(axis=0))  # devrait converger vers b
#%%
n_iter2 = 25
for i in range(n_iter2):
    # Mise à jour de u
    u = a / (P @ v)   # @ est le multiplicateur matricielle
    
    # Mise à jour de v
    v = b / (P.T @ u)  # transposée de la matrice P 
    
    # P1m = a
    # P⊤1n = b
    # elle satisfait les contraintes marginales
    # Calcul du plan de transport courant
    T = np.diag(u) @ P @ np.diag(v)
    
    print(f"\nItération {i+1}")
    print(f"u = {u}")
    print(f"v = {v}")
    print("Plan de transport T:")
    print(T)
    print(f"Marges lignes: {T.sum(axis=1)}")
    print(f"Marges colonnes: {T.sum(axis=0)}")
    #%%
print("P1m = ", T.sum(axis=1))  # devrait converger vers a
print("P⊤1n = ", T.sum(axis=0))  # devrait converger vers b
print(f"T = {T}")
#%%
epsilon = 0.4  # paramètre de régularisation

P = np.exp(-D/epsilon)
u = np.ones(len(a))
v = np.ones(len(b))
for i in range(n_iter):
    # Mise à jour de u
    u = a / (P @ v)   # @ est le multiplicateur matricielle
    
    # Mise à jour de v
    v = b / (P.T @ u)  # transposée de la matrice P 
    
    # P1m = a
    # P⊤1n = b
    # elle satisfait les contraintes marginales
    # Calcul du plan de transport courant
    T = np.diag(u) @ P @ np.diag(v)
    
    print(f"\nItération {i+1}")
    print(f"u = {u}")
    print(f"v = {v}")
    print("Plan de transport T:")
    print(T)
    print(f"Marges lignes: {T.sum(axis=1)}")
    print(f"Marges colonnes: {T.sum(axis=0)}")
    print("P1m = ", T.sum(axis=1))  # devrait converger vers a
    print("P⊤1n = ", T.sum(axis=0))  # devrait converger vers b
    
    # m représente le nombre de points de départ (la dimension de a)
    # n représente le nombre de points d'arrivée (la dimension de b)
    
    # plus le paramètre de régularisation est élevée plus j'obtiens rapidement
    # ma convergence vers a et b
    
    # il faut un nombre d'itération suffisante pour obtenir les distributions
    # marginales souhaitées
    
    # Transport = sinkhorn(a, b, M, reg, method='sinkhorn', numItermax=1000, 
    # stopThr=1e-9,verbose=False, log=False, warn=True, warmstart=None)
    
    # M :matrice des coûts/distances entre chaque paire de points
    # warmstart=None : permet de fournir une initialisation spécifique pour les
    # vecteurs u et v (None signifie initialisation par défaut avec des 1)
    # log=False : si True, renvoie des informations supplémentaires sur le
    # déroulement de l'algorithme
    
    # stopThr=1e-9 : seuil de convergence (l'algorithme s'arrête quand la 
    # différence entre deux itérations est inférieure à ce seuil)
    # Le seuil de convergence (stopThr) représente un critère d'arrêt qui mesure 
    # la stabilité de la solution entre deux itérations successives.
    # ||Tk+1 - Tk|| / ||Tk|| < stopThr
    # max|Pk1 - a| < stopThr et max|Pk^T1 - b| < stopThr
    # ||uk+1 - uk|| / ||uk|| < stopThr et ||vk+1 - vk|| / ||vk|| < stopThr
    
    # reg : paramètre de régularisation (epsilon)
    # Plus il est petit, plus la solution est proche du transport optimal non 
    # régularisé
    # La régularisation dans le transport optimal ajoute un terme d'entropie qui rend le problème plus "lisse" et plus facile à résoudre numériquement. Voici ce que cela implique :
    # Avec un petit epsilon (faible régularisation) :

    # ✓ Solution plus proche du transport optimal "pur"
    # ✓ Coût de transport plus faible
    # ✗ Convergence plus lente
    # ✗ Risque d'instabilités numériques
 
    # Avec un grand epsilon (forte régularisation) :

    # ✓ Convergence plus rapide
    # ✓ Plus stable numériquement
    # ✗ Solution plus "diffuse"
    # ✗ Coût de transport plus élevé
    
    # Avec epsilon = 0.001
    # P1 = exp(-5/0.001) ≈ 0.0000....(beaucoup de zéros)
    # L'ordinateur pourrait arrondir à 0, et alors 1/P1 = inf !

    # Avec epsilon = 0.1
    # P2 = exp(-5/0.1) ≈ 0.0067
    # Nombre plus "raisonnable" pour les calculs
 #%%
epsilon = 0.004  # paramètre de régularisation
# pour tester l'instabilité dans mon exemple on ne peut pas aller en dessous
# de 0.004, à 0.003 cela plante

P = np.exp(-D/epsilon)
u = np.ones(len(a))
v = np.ones(len(b))
for i in range(n_iter):
    # Mise à jour de u
    u = a / (P @ v)   # @ est le multiplicateur matricielle
    
    # Mise à jour de v
    v = b / (P.T @ u)  # transposée de la matrice P 
    
    # P1m = a
    # P⊤1n = b
    # elle satisfait les contraintes marginales
    # Calcul du plan de transport courant
    T = np.diag(u) @ P @ np.diag(v)
    
    print(f"\nItération {i+1}")
    print(f"u = {u}")
    print(f"v = {v}")
    print("Plan de transport T:")
    print(T)
    print(f"Marges lignes: {T.sum(axis=1)}")
    print(f"Marges colonnes: {T.sum(axis=0)}")
    print("P1m = ", T.sum(axis=1))  # devrait converger vers a
    print("P⊤1n = ", T.sum(axis=0))  # devrait converger vers b





   
    # Une approche courante est de :
    # Commencer avec un grand epsilon pour obtenir une première solution rapidement
    # Puis diminuer progressivement epsilon pour raffiner la solution
    # S'arrêter quand on atteint un bon compromis entre précision et stabilité
    
