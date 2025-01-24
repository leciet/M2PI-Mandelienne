library(shiny)
library(ollamar)
library(tidyverse)
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

ollamar::pull("llama3")

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
      gène = gene,
      maladie = disease,
      score_moyen = mean_score,
      conf_inf = ci_lower,
      conf_sup = ci_upper,
      niv_conf = conf_level,
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

# Charger les matrices d'assignation depuis un fichier .RData
load("association_matrices_2.RData") 

ui <- fluidPage(
  titlePanel("Analyse de matrices d'association"),
  navbarPage('Menu',
             tabPanel("Matrices d'associations",
                      sidebarLayout(
                        sidebarPanel(
                          checkboxGroupInput("matrices", "Sélectionnez les matrices que vous souhaitez utiliser:", 
                                             choices = names(association_matrices)),
                          sliderInput("seuil", "Choisir le seuil (>=) :", min = 1, max = 1, value = 1,step = 1),
                          actionButton("apply_filter", "Création du dataframe"),
                          downloadButton("download_data", "Télécharger le dataframe")
                        ),
                        
                        mainPanel(
                          tableOutput("result_table")
                        )
                      )),
             tabPanel("Associations", 
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
                        "Note : Cette application nécessite d'avoir installé l'application Ollama en local pour fonctionner. \n Pour en savoir plus rendez-vous sur",tags$a(href="https://ollama.com/download","https://ollama.com/download",target = "_blank")  
                      ),
                      # Conditional panel to show message if data frame is missing
                      conditionalPanel(
                        condition = "!output.dataframe_created",  # Dynamic check for data frame
                        div(
                          style = "background-color: #fff3cd; color: #856404; padding: 10px; margin: 10px 0; border-radius: 5px;",
                          "Important : Vous devez d'abord créer le dataframe dans l'onglet 'Matrices d'associations' avant d'utiliser les analyses dans cet onglet."
                        )
                      ),
                      
                      tabsetPanel(
                        tabPanel("Recherche par Maladie",
                                 sidebarLayout(
                                   sidebarPanel(
                                     h4("Sélectionnez une maladie et au moins un gène puis cliquez sur 'Lancer l'analyse'"),
                                     selectInput("disease1", "Sélectionner une maladie:", choices = NULL),
                                     selectizeInput("genes1", "Sélectionner un ou plusieurs gènes:", choices = NULL, multiple = TRUE),
                                     numericInput("n_scores1", "Nombre de scores par association:", value = 5, min = 1, max = 10),
                                     sliderInput("conf_level1", "Niveau de confiance:", min = 0.8, max = 0.99, value = 0.95),
                                     actionButton("analyze1", "Lancer l'analyse", class = "btn-primary"),
                                     downloadButton("download_results1", "Télécharger les résultats")
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
                        
                        tabPanel("Recherche par Gène",
                                 sidebarLayout(
                                   sidebarPanel(
                                     h4("Sélectionnez un gène et au moins une maladie puis cliquez sur 'Lancer l'analyse'"),
                                     selectInput("gene2", "Sélectionner un gène:", choices = NULL),
                                     selectizeInput("diseases2", "Sélectionner une ou plusieurs maladies:", choices = NULL, multiple = TRUE),
                                     numericInput("n_scores2", "Nombre de scores par association:", value = 5, min = 1, max = 10),
                                     sliderInput("conf_level2", "Niveau de confiance:", min = 0.8, max = 0.99, value = 0.95),
                                     actionButton("analyze2", "Lancer l'analyse", class = "btn-primary"),
                                     downloadButton("download_results2", "Télécharger les résultats")
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
  )
)

server <- function(input, output, session) {
  
  # Store cumulative results
  cumulative_results1 <- reactiveVal(data.frame(
    gène = character(), maladie = character(),
    score_moyen = numeric(), conf_inf = numeric(),
    conf_sup = numeric(), niv_conf = numeric(),
    matrices_util = character(), seuil = numeric(),
    stringsAsFactors = FALSE
  ))
  
  cumulative_results2 <- reactiveVal(data.frame(
    gène = character(), maladie = character(),
    score_moyen = numeric(), conf_inf = numeric(),
    conf_sup = numeric(), niv_conf = numeric(),
    matrices_util = character(), seuil = numeric(),
    stringsAsFactors = FALSE
  ))
  
  reasons1 <- reactiveVal(character())
  reasons2 <- reactiveVal(character())
  
  # Helper to append matrix/threshold data
  process_results <- function(results, matrices, threshold) {
    results %>%
      mutate(
        matrices_util = paste(matrices, collapse = ", "),
        seuil = threshold
      )
  }
  # Matrices réactives
  matrices_reactive <- reactive({
    association_matrices
  })
  
  # Combiner les matrices sélectionnées
  combined_matrix <- reactive({
    req(input$matrices)
    selected_matrices <- matrices_reactive()[input$matrices]
    Reduce(`+`, lapply(selected_matrices, as.matrix))
  })
  
  # Dataframe filtré
  df_tidy <- eventReactive(input$apply_filter, {
    req(combined_matrix())
    combined_matrix() %>%
      as.data.frame() %>%
      mutate(maladie = rownames(.)) %>%
      pivot_longer(-maladie, names_to = "gène", values_to = "Association") %>%
      filter(Association >= input$seuil)
  })
  
  # Update "dataframe_created" status
  output$dataframe_created <- reactive({
    !is.null(df_tidy()) && nrow(df_tidy()) > 0
  })
  outputOptions(output, "dataframe_created", suspendWhenHidden = FALSE)
  
  # Display filtered data frame
  output$result_table <- renderTable({
    req(df_tidy())
    df_tidy()
  })
  
  # Download filtered data frame
  output$download_data <- downloadHandler(
    filename = function() {
      paste("filtered_dataframe_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(df_tidy(), file, row.names = FALSE)
    }
  )
  
  # Vérifier si Ollama est installé
  output$ollama_installed <- reactive({
    check_ollama()
  })
  outputOptions(output, "ollama_installed", suspendWhenHidden = FALSE)
  
  # Handle dynamic updates for inputs
  observe({
    req(df_tidy())
    unique_diseases <- unique(df_tidy()$maladie)
    updateSelectInput(session, "disease1", choices = unique_diseases)
  })
  
  observe({
    req(input$matrices)  # Ensure matrices are selected
    max_value <- length(input$matrices)
    updateSliderInput(session, "seuil", max = max_value)
  })
  
  observe({
    req(input$disease1, df_tidy())
    filtered_data <- df_tidy() %>% filter(maladie == input$disease1)
    unique_genes <- unique(filtered_data$gène)
    updateSelectizeInput(session, "genes1", choices = unique_genes)
  })
  
  observe({
    req(df_tidy())
    unique_genes <- unique(df_tidy()$gène)
    updateSelectInput(session, "gene2", choices = unique_genes)
  })
  
  observe({
    req(input$gene2, df_tidy())
    filtered_data <- df_tidy() %>% filter(gène == input$gene2)
    unique_diseases <- unique(filtered_data$maladie)
    updateSelectizeInput(session, "diseases2", choices = unique_diseases)
  })
  
  # Append results for Tab 1
  observeEvent(input$analyze1, {
    req(input$disease1, input$genes1)
    
    # Vérifier si Ollama est installé
    if (!check_ollama()) {
      showNotification("Ollama n'est pas installé. Veuillez l'installer avant de continuer.", type = "error")
      return()
    }
    
    # Afficher une progression avec message d'erreur si nécessaire
    tryCatch({
      withProgress(message = 'Analyse en cours...', value = 0, {
        all_results <- lapply(input$genes1, function(g) {
          incProgress(1 / length(input$genes1))
          analyze_gene_disease(
            gene = g,
            disease =  input$disease1,
            n_scores = input$n_scores1,
            conf_level = input$conf_level1
          )
        })
        
        # Sauvegarder les résultats
        new_results <- do.call(rbind, lapply(all_results, `[[`, "results"))
        new_results <- process_results(new_results, input$matrices, input$seuil)
        cumulative_results1(bind_rows(cumulative_results1(), new_results))
        
        # Ajouter des raisons
        reasons1(c(
          lapply(seq_along(all_results), function(i) {
            sprintf("Pour %s - %s :\n%s",
                    all_results[[i]]$results$gène,
                    all_results[[i]]$results$maladie,
                    paste(all_results[[i]]$reasons, collapse = "\n"))
          }),
          reasons1()
        ))
      })
    }, error = function(e) {
      # Si une erreur survient, afficher une notification
      showNotification(
        paste("Erreur pendant l'analyse :", conditionMessage(e)), 
        type = "error", 
        duration = 10
      )
    })
  })
  
  # Append results for Tab 2
  observeEvent(input$analyze2, {
    req(input$diseases2, input$gene2)
    
    # Vérifier si Ollama est installé
    if (!check_ollama()) {
      showNotification("Ollama n'est pas installé. Veuillez l'installer avant de continuer.", type = "error")
      return()
    }
    
    # Afficher une progression avec message d'erreur si nécessaire
    tryCatch({
      withProgress(message = 'Analyse en cours...', value = 0, {
        all_results <- lapply(input$gene2, function(g) {
          incProgress(1 / length(input$gene2))
          analyze_gene_disease(
            gene = g,
            disease = input$diseases2,
            n_scores = input$n_scores2,
            conf_level = input$conf_level2
          )
        })
        
        # Sauvegarder les résultats
        new_results <- do.call(rbind, lapply(all_results, `[[`, "results"))
        new_results <- process_results(new_results, input$matrices, input$seuil)
        cumulative_results2(bind_rows(cumulative_results2(), new_results))
        
        # Ajouter des raisons
        reasons2(c(
          lapply(seq_along(all_results), function(i) {
            sprintf("Pour %s - %s :\n%s",
                    all_results[[i]]$results$gène,
                    all_results[[i]]$results$maladie,
                    paste(all_results[[i]]$reasons, collapse = "\n"))
          }),
          reasons2()
        ))
      })
    }, error = function(e) {
      # Si une erreur survient, afficher une notification
      showNotification(
        paste("Erreur pendant l'analyse :", conditionMessage(e)), 
        type = "error", 
        duration = 10
      )
    })
  })
  
  # Display and download results for Tab 1
  output$results_table1 <- renderTable({
    req(cumulative_results1())
    validate(
      need(nrow(cumulative_results1()) > 0, "Aucun résultat disponible. Veuillez réessayer.")
    )
    cumulative_results1()
  })
  
  
  output$download_results1 <- downloadHandler(
    filename = function() {
      paste("results_table1_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(cumulative_results1(), file, row.names = FALSE)
    }
  )
  
  output$reasons1 <- renderText({
    req(reasons1())
    paste(reasons1(), collapse = "\n\n")  # Join entries with double newline
  })
  
  # Display and download results for Tab 2
  output$results_table2 <- renderTable({
    req(cumulative_results2())
    validate(
      need(nrow(cumulative_results2()) > 0, "Aucun résultat disponible. Veuillez réessayer.")
    )
    cumulative_results2()
  })
  
  output$download_results2 <- downloadHandler(
    filename = function() {
      paste("results_table2_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(cumulative_results2(), file, row.names = FALSE)
    }
  )
  output$reasons2 <- renderText({
    req(reasons2())
    paste(reasons2(), collapse = "\n\n")  # Double saut de ligne pour séparer chaque entrée
  })
}

shinyApp(ui = ui, server = server)