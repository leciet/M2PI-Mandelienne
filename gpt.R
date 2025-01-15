###### LLM ollamar

install.packages("ollamar")
library(ollamar)

pull("cniongolo/biomistral") 
pull("meditron") 
pull("llama3")

#Prompt  1 : 

disease <- "Adrenal hypofunction"
gene <- "CACNA1H"
prompt <- sprintf(
  "We've conducted an analysis and found an association between %s and %s. Please provide a score of confidence, from 1 to 10, for the suggested causal relationship between %s and gene %s",
  disease, gene, disease, gene
)

#Prompt  2 : 

disease <- "Adrenal hypofunction"
gene <- "CACNA1H"
prompt <- sprintf(
  "Give me the confidence score (1-10) for the association between %s and gene %s.
  Score the confidence level based on the following criteria:
  - 0-2: No or very weak evidence.
  - 3-5: Some studies exist, but results are contradictory.
  - 6-8: Multiple supporting studies, evidence is consistent.
  - 9-10: Strong scientific consensus, supported by validated mechanisms.
  Answer:",
  disease, gene
)

#Prompt 3 : 

disease <- "Adrenal hypofunction"
gene <- "CACNA1H"
prompt <- sprintf(
  "Question: Quantify on a scale from 0 to 10 the potential causal relationship between %s variations and %s, considering both direct evidence and mechanistic plausibility.

0-2: No or very weak scientific evidence and low mechanistic plausibility 
3-5: Limited direct evidence but significant mechanistic plausibility based on known pathways 
6-8: Strong mechanistic evidence and/or some direct clinical evidence 
9-10: Well-established causal relationship with validated mechanisms

Consider:
Direct evidence 
Molecular pathway analysis
Physiological relevance
Comparative analysis: similar genes/pathways 
  
RESPOND EXACTLY LIKE THIS TEMPLATE:
SCORE: [enter single number 0-10] 
REASON: [enter one short sentence]
",gene, disease, gene, disease
)

#Prompt Sebastien 
prompt <- "Based on the function of the CACNA1H gene, could this gene be involved in the development of adrenal hypofunction? Please provide a confidence score from 0 to 5, where 0 means no involvement at all and 5 means definite involvement.  
"
#Final 
library(stringr)
responses <- vector("character", 50)
disease <- "Adrenal hypofunction"
gene <- "CACNA1H"

prompt <- sprintf("Question: Based on the function of the %s gene, could this gene be involved in the development of %s? Quantify on a scale from 0 to 10 the potential causal relationship between %s variations and %s, considering both direct evidence and mechanistic plausibility.

0-2: No or very weak scientific evidence and low mechanistic plausibility 
3-5: Limited direct evidence but significant mechanistic plausibility based on known pathways 
6-8: Strong mechanistic evidence and/or some direct clinical evidence 
9-10: Well-established causal relationship with validated mechanisms

Consider:
Direct evidence 
Molecular pathway analysis
Physiological relevance
Comparative analysis: similar genes/pathways 
  
RESPOND EXACTLY LIKE THIS TEMPLATE:
SCORE: [enter single number 0-10] 
REASON: [enter one short sentence]
",gene, disease, gene, disease
)

for (i in 1:50){
  responses[i] <- generate("llama3", prompt = prompt, output = "text")
}

#Fonction automatisation

# Fonction pour collecter les scores
collect_scores <- function(prompt, gene, disease, n_scores = 5, max_attempts = 50, conf_level = 0.95) {
  scores <- numeric()
  attempts <- 0
  
  while (length(scores) < n_scores && attempts < max_attempts) {
    attempts <- attempts + 1
    tryCatch({
      response <- generate("llama3", prompt = prompt, output = "text")
      score <- as.numeric(sub(".*SCORE:\\s*(\\d+).*", "\\1", response))

      if (!is.na(score)) {
        scores <- c(scores, score)
      }
    }, error = function(e) {
      warning(sprintf("Tentative %d échouée : %s", attempts, e$message))
    })
    
    Sys.sleep(1)
  }
  
  if (length(scores) < n_scores) {
    stop(sprintf("Impossible d'obtenir %d scores valides après %d tentatives", 
                 n_scores, attempts))
  }
  
  # Calcul des statistiques
  mean_score <- mean(scores)
  se <- sd(scores) / sqrt(length(scores))
  t_value <- qt(1 - (1 - conf_level) / 2, df = length(scores) - 1)
  ci_lower <- mean_score - t_value * se 
  ci_upper <- mean_score + t_value * se
  
  # Retourne un dataframe avec une ligne
  return(data.frame(
    gene = gene,
    disease = disease,
    mean_score = mean_score,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    conf_level = conf_level,
    stringsAsFactors = FALSE
  ))
}

# Fonction pour analyser un pair gène-maladie
analyze_gene_disease <- function(gene, disease, n_scores = 5, conf_level = 0.95) {
  # Création du prompt avec le format original
  prompt <- sprintf("Question: Based on the function of the %s gene, could this gene be involved in the development of %s? Quantify on a scale from 0 to 10 the potential causal relationship between %s variations and %s, considering both direct evidence and mechanistic plausibility.
0-2: No or very weak scientific evidence and low mechanistic plausibility 
3-5: Limited direct evidence but significant mechanistic plausibility based on known pathways 
6-8: Strong mechanistic evidence and/or some direct clinical evidence 
9-10: Well-established causal relationship with validated mechanisms

Consider:
1. Direct evidence 
2. Molecular pathway analysis
3. Physiological relevance
4. Comparative analysis: similar genes/pathways 
  
RESPOND EXACTLY LIKE THIS TEMPLATE:
SCORE: [enter single number 0-10] 
REASON: [enter one short sentence]", 
                    gene, disease, gene, disease)
  
  # Collecte des scores et retourne directement le dataframe d'une ligne
  return(collect_scores(prompt, gene, disease, n_scores, conf_level = conf_level))
}

# Application sur le jeu de test
test <- df_tidy[sample(nrow(df_tidy), 1), ]
system.time(results <- do.call(rbind, apply(test, 1, function(row) {
  analyze_gene_disease(
    gene = row["gene"], 
    disease = row["disease"]
  )
}))
)

