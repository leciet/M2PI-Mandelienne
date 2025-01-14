# ------------------------------------------------------------------------------
# LLM via ollama
# 10/01/2025
# Leslie Cieters
# ------------------------------------------------------------------------------

# Nettoyage environnement
# rm(list=ls())

## On souhaite utiliser des LLM pour attribuer un score allant de 0 à 10 à    ##
## chaque association gène-maladie (depuis la matrice d'assignation issu du   ##
## seuil).                                                                    ##

# Packages ----
if (!require("ollamar")) install.packages("ollamar")

library(tidyverse)
library(ollamar)

# ollamar ----

# Vérification de la connexion
test_connection()

# meditron ---- 

## https://ollama.com/library/meditron                                        ##
## Open-source medical large language model adapted from Llama 2 to the       ##
## medical domain.                                                            ##

# Téléchargement du modèle
ollamar::pull("meditron")

prompt <-  "Give me the confidence level ( from 0 to 10) of the association between Adrenal hypofunction and the gene CACNA1H on a scale from 1 to 10. The confidence level should be based on peer-reviewed scientific publications, the reproducibility of results, and the biological plausibility of the association. \n
Evaluation criteria:\n
0-2: No or very weak scientific evidence.\n
3-5: Some studies but contradictory results.\n
6-8: Multiple concordant studies.\n
9-10: Strong scientific consensus with validated mechanisms.\n
Additionally, even if there are no explicit studies directly linking these two, please assess whether the underlying biological mechanisms involving CACNA1H are plausible in the context of adrenal function. Consider whether the biological role of CACNA1H in calcium signaling in adrenal cells might logically relate to adrenal hypofunction"
res <-  generate("llama3", prompt=prompt, output="df")
for (i in 2:10) {
  res[i,] <- generate("llama3", prompt=prompt, output="df")
}

 res

# 2ème test

prompt2 <- "Question: Give me a confidence level on a scale from 0 to 10 for the association between Adrenal hypofunction and the gene CACNA1H without explanation.
0-2: No or very weak scientific evidence.
3-5: Some studies but contradictory results.
6-8: Multiple concordant studies.
9-10: Strong scientific consensus with validated mechanisms.
Additionally, even if there are no explicit studies directly linking these two, 
please assess whether the underlying biological mechanisms involving CACNA1H are plausible in the context of adrenal function. Consider whether the biological role of CACNA1H in calcium signaling in adrenal cells might logically relate to adrenal hypofunction
"
res2 <- generate("cniongolo/biomistral", prompt=prompt2, output="df")
for (i in 2:10) {
  res2[i,] <- generate("cniongolo/biomistral", prompt=prompt2, output="df")
}

res
ollamar::list_models()


# Test de prompt modifié avec un autre vocabulaire

prompt3 <- "We suppose a relationship exist between CACNA1H variation and adrenal hypofunction. Give me a percentage to describe the confidence we have on this supposed causal relationship between CACNA1H variation and Adrenal hypofunction."

res3 <- generate("cniongolo/biomistral", prompt=prompt3, output="df")
for (i in 2:10) {
  res3[i,] <- generate("cniongolo/biomistral", prompt=prompt3, output="df")
}

# Test de prompt modifié pour les sources

prompt4 <- "We suppose a relationship exist between CACNA1H variation and adrenal hypofunction. Give me references (title, author and date)."

res4 <- generate("cniongolo/biomistral", prompt=prompt4, output="df")
for (i in 2:10) {
  res4[i,] <- generate("cniongolo/biomistral", prompt=prompt4, output="df")
}

# test avec llama3

pull('llama3')
