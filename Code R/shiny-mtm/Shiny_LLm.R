library(shiny)
library(dplyr)
library(ollamar)
library(stringr)
library(sys)

# Vérification ollama en local
check_ollama <- function() {
  tryCatch({
    system2("ollama", "list", stdout = TRUE, stderr = TRUE)
    return(TRUE)
  }, error = function(e) {
    return(FALSE)
  })
}

pull("llama3")

# Fonction pour collecter les scores et les raisons
collect_scores <- function(prompt, gene, disease, n_scores = 5, max_attempts = 50, conf_level = 0.95) {
  scores <- numeric()
  reasons <- character()
  attempts <- 0
  
  while (length(scores) < n_scores && attempts < max_attempts) {
    attempts <- attempts + 1
    tryCatch({
      response <- generate("llama3", prompt = prompt, output = "text")
      score <- as.numeric(sub(".*SCORE:\\s*(\\d+).*", "\\1", response))
      reason <- sub(".*REASON:\\s*(.*)$", "\\1", response)
      
      if (!is.na(score)) {
        scores <- c(scores, score)
        reasons <- c(reasons, sprintf("Iteration %d : %s", length(scores), reason))
      }
    }, error = function(e) {
      warning(sprintf("Tentative %d échouée : %s", attempts, e$message))
    })
    
    Sys.sleep(1)
  }
  
  if (length(scores) < n_scores) {
    stop(sprintf("Impossible d'obtenir %d scores valides après %d tentatives", n_scores, attempts))
  }
  
  mean_score <- mean(scores)
  se <- sd(scores) / sqrt(length(scores))
  t_value <- qt(1 - (1 - conf_level) / 2, df = length(scores) - 1)
  ci_lower <- mean_score - t_value * se 
  ci_upper <- mean_score + t_value * se
  
  return(list(
    results = data.frame(
      gene = gene,
      disease = disease,
      mean_score = mean_score,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      conf_level = conf_level,
      stringsAsFactors = FALSE
    ),
    reasons = reasons
  ))
}

# Fonction pour analyser les paires gène-maladie
analyze_gene_disease <- function(gene, disease, n_scores = 5, conf_level = 0.95) {
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
  
  return(collect_scores(prompt, gene, disease, n_scores, conf_level = conf_level))
}

# Interface utilisateur
ui <- fluidPage(
  titlePanel("Analyse des associations gène-maladie"),
  
  # Message d'avertissement Ollama
  conditionalPanel(
    condition = "!output.ollama_installed",
    div(
      style = "background-color: #f8d7da; color: #721c24; padding: 10px; margin: 10px 0; border-radius: 5px;",
      "Veuillez télécharger Ollama sur ",
      tags$a(href = "https://ollama.com/", "https://ollama.com/", target = "_blank")
    )
  ),
  
  div(
    style = "background-color: #e2e3e5; padding: 10px; margin: 10px 0; border-radius: 5px;",
    "Note : Cette application nécessite d'avoir installé l'application Ollama en local pour fonctionner."
  ),
  
  tabsetPanel(
    # Onglet 1 : Recherche par maladie
    tabPanel("Recherche par Maladie",
             sidebarLayout(
               sidebarPanel(
                 h4("Sélectionnez une maladie et au moins un gène puis cliquez sur 'Lancer l'analyse'"),
                 selectInput("disease1", 
                             "Sélectionner une maladie:",
                             choices = NULL),
                 
                 selectizeInput("genes1",
                                "Sélectionner un ou plusieurs gènes:",
                                choices = NULL,
                                multiple = TRUE),
                 
                 actionButton("analyze1", "Lancer l'analyse",
                              class = "btn-primary"),
                 
                 numericInput("n_scores1", 
                              "Nombre de scores par association:",
                              value = 5,
                              min = 1,
                              max = 10),
                 
                 sliderInput("conf_level1",
                             "Niveau de confiance:",
                             min = 0.8,
                             max = 0.99,
                             value = 0.95)
               ),
               
               mainPanel(
                 h3("Résultats de l'analyse"),
                 tableOutput("results_table1"),
                 hr(),
                 h4("Raisons détaillées :"),
                 verbatimTextOutput("reasons1"),
                 textOutput("progress1")
               )
             )
    ),
    
    # Onglet 2 : Recherche par gène
    tabPanel("Recherche par Gène",
             sidebarLayout(
               sidebarPanel(
                 h4("Sélectionnez un gène et au moins une maladie puis cliquez sur 'Lancer l'analyse'"),
                 selectInput("gene2", 
                             "Sélectionner un gène:",
                             choices = NULL),
                 
                 selectizeInput("diseases2",
                                "Sélectionner une ou plusieurs maladies:",
                                choices = NULL,
                                multiple = TRUE),
                 
                 actionButton("analyze2", "Lancer l'analyse",
                              class = "btn-primary"),
                 
                 numericInput("n_scores2", 
                              "Nombre de scores par association:",
                              value = 5,
                              min = 1,
                              max = 10),
                 
                 sliderInput("conf_level2",
                             "Niveau de confiance:",
                             min = 0.8,
                             max = 0.99,
                             value = 0.95)
               ),
               
               mainPanel(
                 h3("Résultats de l'analyse"),
                 tableOutput("results_table2"),
                 hr(),
                 h4("Raisons détaillées :"),
                 verbatimTextOutput("reasons2"),
                 textOutput("progress2")
               )
             )
    )
  )
)

# Serveur
server <- function(input, output, session) {
  # Vérifier si Ollama est installé
  output$ollama_installed <- reactive({
    check_ollama()
  })
  outputOptions(output, "ollama_installed", suspendWhenHidden = FALSE)
  
  # Charger les données
  df_tidy <- reactive({
    read.csv("df_tidy.csv", stringsAsFactors = FALSE)
  })
  
  # Variables réactives pour les résultats et les raisons
  results1 <- reactiveVal(NULL)
  results2 <- reactiveVal(NULL)
  reasons1 <- reactiveVal(NULL)
  reasons2 <- reactiveVal(NULL)
  
  # Variables réactives pour suivre l'état de progression
  progress_msg1 <- reactiveVal(NULL)
  progress_msg2 <- reactiveVal(NULL)
  
  # Mise à jour des choix pour l'onglet 1
  observe({
    updateSelectInput(session, "disease1",
                      choices = sort(unique(df_tidy()$disease)))
  })
  
  observe({
    req(input$disease1)
    genes <- df_tidy()[df_tidy()$disease == input$disease1, "gene"]
    updateSelectizeInput(session, "genes1",
                         choices = sort(unique(genes)))
  })
  
  # Mise à jour des choix pour l'onglet 2
  observe({
    updateSelectInput(session, "gene2",
                      choices = sort(unique(df_tidy()$gene)))
  })
  
  observe({
    req(input$gene2)
    diseases <- df_tidy()[df_tidy()$gene == input$gene2, "disease"]
    updateSelectizeInput(session, "diseases2",
                         choices = sort(unique(diseases)))
  })
  
  # Gestion de l'analyse pour l'onglet 1
  observeEvent(input$analyze1, {
    req(input$disease1, input$genes1)
    
    if (!check_ollama()) {
      showNotification("Ollama n'est pas installé. Veuillez l'installer avant de continuer.", 
                       type = "error")
      return()
    }
    
    progress_msg1("Analyse en cours...")  # Message au début de l'analyse
    
    withProgress(message = 'Analyse en cours...', value = 0, {
      all_results <- lapply(input$genes1, function(g) {
        incProgress(1/length(input$genes1))
        analyze_gene_disease(
          gene = g,
          disease = input$disease1,
          n_scores = input$n_scores1,
          conf_level = input$conf_level1
        )
      })
      
      # Séparer les résultats et les raisons
      results1(do.call(rbind, lapply(all_results, function(x) x$results)))
      reasons1(lapply(seq_along(all_results), function(i) {
        sprintf("Pour %s - %s :\n%s\n",
                all_results[[i]]$results$gene,
                all_results[[i]]$results$disease,
                paste(all_results[[i]]$reasons, collapse = "\n"))
      }))
      
      progress_msg1("Analyse terminée")  # Mettre à jour le message une fois terminé
    })
  })
  
  # Gestion de l'analyse pour l'onglet 2
  observeEvent(input$analyze2, {
    req(input$gene2, input$diseases2)
    
    if (!check_ollama()) {
      showNotification("Ollama n'est pas installé. Veuillez l'installer avant de continuer.", 
                       type = "error")
      return()
    }
    
    progress_msg2("Analyse en cours...")  # Message au début de l'analyse
    
    withProgress(message = 'Analyse en cours...', value = 0, {
      all_results <- lapply(input$diseases2, function(d) {
        incProgress(1/length(input$diseases2))
        analyze_gene_disease(
          gene = input$gene2,
          disease = d,
          n_scores = input$n_scores2,
          conf_level = input$conf_level2
        )
      })
      
      # Séparer les résultats et les raisons
      results2(do.call(rbind, lapply(all_results, function(x) x$results)))
      reasons2(lapply(seq_along(all_results), function(i) {
        sprintf("Pour %s - %s :\n%s\n",
                all_results[[i]]$results$gene,
                all_results[[i]]$results$disease,
                paste(all_results[[i]]$reasons, collapse = "\n"))
      }))
      
      progress_msg2("Analyse terminée")  # Mettre à jour le message une fois terminé
    })
  })
  
  # Affichage des tableaux de résultats
  output$results_table1 <- renderTable({
    req(results1())
    results1() %>%
      mutate(
        mean_score = round(mean_score, 2),
        ci_lower = round(ci_lower, 2),
        ci_upper = round(ci_upper, 2)
      ) %>%
      rename("Gène" = gene,
             "Maladie" = disease,
             "Score Moyen" = mean_score,
             "IC Bas" = ci_lower,
             "IC Haut" = ci_upper)
  })
  
  output$results_table2 <- renderTable({
    req(results2())
    results2() %>%
      mutate(
        mean_score = round(mean_score, 2),
        ci_lower = round(ci_lower, 2),
        ci_upper = round(ci_upper, 2)
      ) %>%
      rename("Gène" = gene,
             "Maladie" = disease,
             "Score Moyen (0 to 10)" = mean_score,
             "IC Bas" = ci_lower,
             "IC Haut" = ci_upper)
  })
  
  # Affichage des raisons détaillées
  output$reasons1 <- renderText({
    req(reasons1())
    paste(reasons1(), collapse = "\n")
  })
  
  output$reasons2 <- renderText({
    req(reasons2())
    paste(reasons2(), collapse = "\n")
  })
  
  # Messages de progression
  output$progress1 <- renderText({
    progress_msg1()
  })
  
  output$progress2 <- renderText({
    progress_msg2()
  })
}

# Lancement de l'application
shinyApp(ui = ui, server = server)
