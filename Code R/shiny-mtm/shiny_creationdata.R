library(shiny)
library(dplyr)
library(tidyr)

# Charger les matrices d'assignation depuis un fichier .RData
load("association_matrices_2.RData") 

ui <- fluidPage(
  titlePanel("Analyse de matrices d'association"),
  
  sidebarLayout(
    sidebarPanel(
      checkboxGroupInput("matrices", "Sélectionnez les matrices que vous souhaitez utiliser:", 
                         choices = names(association_matrices)), # Ajout des noms directement
      numericInput("seuil", "Choisir le seuil (>=) :", min = 1, max = length(association_matrices), value = 1),
      actionButton("apply_filter", "Création du dataframe"),
      downloadButton("download_data", "Télécharger le dataframe")
    ),
    
    mainPanel(
      tableOutput("result_table")
    )
  )
)

server <- function(input, output, session) {
  
  # Regrouper les matrices dans une liste statique
  matrices <- reactive({
    association_matrices
  })
  
  # Combiner les matrices sélectionnées
  combined_matrix <- reactive({
    req(input$matrices) # Vérifie qu'au moins une matrice est sélectionnée
    selected_matrices <- matrices()[input$matrices]
    
    # Sommation élément par élément
    Reduce(`+`, lapply(selected_matrices, as.matrix))
  })
  
  # Filtrer selon le seuil et convertir au format long
  filtered_table <- eventReactive(input$apply_filter, {
    req(combined_matrix())
    
    # Conversion au format long
    combined_matrix() %>%
      as.data.frame() %>%
      mutate(disease = rownames(.)) %>%
      pivot_longer(-disease, names_to = "gene", values_to = "Association") %>%
      filter(Association >= input$seuil)
  })
  
  # Sauvegarder le tableau filtré comme un dataframe réactif 'df_tidy'
  df_tidy <- reactiveVal(NULL)  # Créer une variable réactive pour df_tidy
  
  observeEvent(input$apply_filter, {
    df_tidy(filtered_table())  # Mettre à jour 'df_tidy' avec les données filtrées
  })
  
  # Affichage dans l'interface utilisateur
  output$result_table <- renderTable({
    req(df_tidy())  # Attendre que df_tidy soit disponible
    df_tidy()  # Afficher df_tidy
  }, rownames = TRUE) # Inclure les noms des lignes dans le tableau
  
  # Fonction pour télécharger le dataframe filtré
  output$download_data <- downloadHandler(
    filename = function() {
      paste("filtered_dataframe_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(df_tidy(), file)  
    }
  )
}

shinyApp(ui = ui, server = server)
