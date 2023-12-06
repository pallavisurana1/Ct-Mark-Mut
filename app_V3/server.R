
# server.R

##----------------------------------------------------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------------------------------------------------

classify_mutation <- function(mutation) {
  if (grepl("^c\\.\\d+[ACGT]>[ACGT]", mutation)) {
    # Point mutations (substitutions)
    return("Point Mutation")
  } else if (grepl("del", mutation)) {
    # Deletions
    return("Deletion")
  } else if (grepl("ins", mutation)) {
    # Insertions
    return("Insertion")
  } else if (grepl("^c\\.\\d+\\+\\d+", mutation) || grepl("^c\\.\\d+\\-\\d+", mutation)) {
    # Splice site mutations
    return("Splice Site Mutation")
  } else if (grepl("dup", mutation)) {
    # Duplications
    return("Duplication")
  } else if (grepl("\\*", mutation)) {
    # Mutations at the stop codon
    return("Stop Codon Mutation")
  } else if (grepl("c\\.\\-\\d+", mutation)) {
    # Promoter region mutations
    return("Promoter Mutation")
  } else if (grepl("c\\.\\*\\d+", mutation)) {
    # Mutations in the untranslated region (UTR)
    return("UTR Mutation")
  } else if (grepl("c\\.[\\d_]+[ACGT]+>[ACGT]+", mutation)) {
    # Complex substitutions
    return("Complex Substitution")
  } else if (grepl("c\\.[\\d_]+del[ACGT]+", mutation)) {
    # Deletions of specific nucleotides
    return("Specific Deletion")
  } else if (grepl("c\\.[\\d_]+ins[ACGT]+", mutation)) {
    # Insertions of specific nucleotides
    return("Specific Insertion")
  } else {
    # Other or unspecified mutations
    return("Other")
  }
}

classify_amino_acid_mutation <- function(mutation) {
  if (grepl("fs", mutation)) {
    return("Frameshift")
  } else if (grepl("\\*", mutation)) {
    return("Nonsense")
  } else if (grepl("=", mutation)) {
    return("Silent")
  } else if (grepl("[A-Z]\\d+[A-Z]", mutation) && !grepl("del|ins", mutation)) {
    return("Missense")
  } else if (grepl("del", mutation) && !grepl("fs", mutation)) {
    return("In-Frame Deletion")
  } else if (grepl("ins", mutation) && !grepl("fs", mutation)) {
    return("In-Frame Insertion")
  } else {
    return("Other")
  }
}


##----------------------------------------------------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------------------------------------------------


# Define server 
server <- function(input, output, session) {

  log_messages <- reactiveVal(character(0))
  shinyjs::toggle("loadingPopup", anim = FALSE)
  Sys.sleep(3)
  
  ##-Function to add a log message
  add_log_message <- function(message) {
    current_messages <- log_messages()
    updated_messages <- c(current_messages, message)
    log_messages(updated_messages)
  }
  
  ##-Function to update loading popup message
  update_loading_message <- function(message) {
    shinyjs::html("loadingPopup", message)
  }
  
  ##-pick the organ and cell type of interest from user
  observeEvent(input$organ1, {
    choices <- cto_all$cell_type[cto_all$organ == input$organ1]
    updateSelectInput(session, "cell_type1", choices = choices)
  })
  
  
  ##-To get the corresponding cosmic_tissue
  reactive_cosmic_tissue <- reactive({
    req(input$organ1, input$cell_type1) 
    filtered_data <- cto_all[cto_all$organ == input$organ1 & cto_all$cell_type == input$cell_type1, ]
    return(filtered_data$cosmic_tissue)
  })
  
  observeEvent(input$run, {
    
    add_log_message("Retrieve Mutations button clicked")
    Sys.sleep(3)
    
    shinyjs::toggle("loadingPopup", anim = TRUE)
    update_loading_message("Loading data, please wait...")
    Sys.sleep(3)
    
    ##-get cell type marker information from panglaoDb 
    markers <- reactive({
      x = panglaoDB_csv  %>% subset(organ == input$organ1) %>% subset(cell_type == input$cell_type1)
      
      ##-take into account the species that the user picks for the cell types markers information
      if("Mm" %in% input$speciesInput && "Hs" %in% input$speciesInput){
        x <- x %>% subset(species %in% c("Mm Hs"))
      }else if("Mm" %in% input$speciesInput){
        x <- x %>% subset(species == "Mm")
      }else if("Hs" %in% input$speciesInput){
        x <- x %>% subset(species == "Hs")
      }
      
      x
    })
    
    ##-genes to query from COSMIC
    # markers()$official_gene_symbol

    ##- used the COSMIC v4 API here
    ##- query data for one gene from cosmic
    COSMIC_query <- list(
      format = 'json',
      maxList = 500,
      df = 'MutationID,GeneName,MutationCDS,MutationAA,AccessionNumber,PrimaryHistology,PrimarySite,PubmedPMID',
      page = 1,
      page_size = 2000
    )
    
    
    ##-fetch data using cosmic api
    fetch_data <- function(gene_nms, page_limit = 10) {
      all_data <- data.frame()  # initialize an empty data frame
      for (gene in gene_nms) {
        COSMIC_query$terms <- gene
        COSMIC_query$page_size <- 500  # set the page size to the maximum
        page <- 1
        while (TRUE) {
          if (page > page_limit) {
            break  # stop if the page limit has been reached
          }
          COSMIC_query$page <- page
          response <- GET('https://clinicaltables.nlm.nih.gov/api/cosmic/v4/search', query = COSMIC_query)
          data <- content(response, as = "text") %>% fromJSON()
          data.df <- data[[4]] %>% as.data.frame()
          if (nrow(data.df) == 0) {
            break  # stop if there are no more results
          }
          colnames(data.df) <- strsplit(COSMIC_query$df, ",")[[1]]
          all_data <- rbind(all_data, data.df)  # combine the data frames
          all_data = all_data %>% as.data.frame %>% 
                                  dplyr::filter(PrimarySite %in% reactive_cosmic_tissue()) %>% 
                                  as.data.frame
          page <- page + 1
        }
      }
      return(all_data)
    }
    
    
    ##-COSMIC mutation data for cell type markers
    cosmic_df <- reactive({
      req(input$organ1, input$cell_type1)
      
      fetched_df = fetch_data(markers()$official_gene_symbol, page_limit = 10)
      fetched_df = fetched_df[!grepl(";", fetched_df$PrimarySite), ]
      fetched_df = fetched_df[!grepl(";", fetched_df$PrimaryHistology), ]
      
      ##-compute mutation burden for primary site
      mutation_burden <- fetched_df %>%
        group_by(PrimarySite) %>%
        summarise(mut_per_mb_Site = n() / 1000000) %>%
        arrange(mut_per_mb_Site) %>%
        mutate(PrimarySite = forcats::fct_reorder(PrimarySite, mut_per_mb_Site))
      
      ##-merge mutation burden and cosmic fetched data
      merged_Df = merge(fetched_df, mutation_burden, by = "PrimarySite")
      merged_Df$PubmedPMID = NULL
      
      merged_Df$mut_per_mb_Site = -log10(merged_Df$mut_per_mb_Site)
      merged_Df
      
    })
    
    output$mutationTable <- renderDT({
      if (nrow(cosmic_df()) > 0) {
        cosmic_df() %>%
          lapply(function(x) ifelse(nchar(as.character(x)) > 20, paste0(substr(x, 1, 20), "..."), x)) %>%
          as.data.frame() %>%
          DT::datatable(., options = list(scrollX = TRUE, scrollY = "200px", pageLength = -1, dom = 't'))
      } else {
        stop("No mutation data found for marker genes for organism of interest")
        add_log_message("No mutations found")
        Sys.sleep(3)
      }
    })
    
    
    output$markerTable <- renderDT({
      if (nrow(markers()) > 0) {
        markers() %>%
          select_if(~is.character(.)) %>%
          lapply(function(x) ifelse(nchar(as.character(x)) > 20, paste0(substr(x, 1, 20), "..."), x)) %>%
          as.data.frame() %>%
          DT::datatable(., options = list(scrollX = TRUE, scrollY = "200px", pageLength = -1, dom = 't'))
      } else {
        stop("No markers found for cell type chosen in the organism of interest")
        add_log_message("No markers found")
        Sys.sleep(3)
      }
    })
    
    output$downloadMarkerTable <- downloadHandler(
      filename = function() {
        paste0("cell-type_marker_table_", Sys.Date(), ".csv")
      },
      content = function(file) {
        write.csv(markers(), file, row.names = FALSE)
      }
    )
    
    output$downloadMutationTable <- downloadHandler(
      filename = function() {
        paste0("mutations_table_", Sys.Date(), ".csv")
      },
      content = function(file) {
        write.csv(cosmic_df(), file, row.names = FALSE)
      }
    )
    
    ##-AA plot
    add_log_message("Data processing and plot generation started")
    Sys.sleep(3)
    
    output$mutationAA_heatmap <- renderPlotly({
      cosmic_data <- cosmic_df()
      
      aa_mutation_counts <- cosmic_data %>%
        as.data.frame() %>%
        group_by(GeneName, MutationAA) %>%
        summarise(n = n(), .groups = 'drop') %>%
        pivot_wider(names_from = MutationAA, values_from = n, values_fill = list(n = 0))
      
      # Convert to a matrix for the heatmap
      aa_mutation_matrix <- as.matrix(aa_mutation_counts[-1])
      rownames(aa_mutation_matrix) <- aa_mutation_counts$GeneName
      heatmap_data1 <- reshape2::melt(aa_mutation_matrix)
      heatmap_data1 = heatmap_data1 %>% mutate(Classification = sapply(Var2, classify_amino_acid_mutation))
      
      plot_ly(data = heatmap_data1, x = ~Classification, y = ~Var1, z = ~value, type = "heatmap",  colors = "Blues") %>%
        layout(
          title = "Amino Acid Mutation",
          xaxis = list(title = " ",
                       tickangle = 45),  # Custom X-axis label
          yaxis = list(title = " "))
      
    })
    
    ##-cds plot
    output$mutationCDS_heatmap <- renderPlotly({
      cosmic_data <- cosmic_df()
      
      # Count the frequency of each MutationCDS for each gene
      cds_mutation_counts <- cosmic_data %>%
        as.data.frame() %>%
        group_by(GeneName, MutationCDS) %>%
        summarise(n = n(), .groups = 'drop') %>%
        pivot_wider(names_from = MutationCDS, values_from = n, values_fill = list(n = 0))
      
      # Convert to a matrix for the heatmap
      cds_mutation_matrix <- as.matrix(cds_mutation_counts[-1])
      rownames(cds_mutation_matrix) <- cds_mutation_counts$GeneName
      heatmap_data <- reshape2::melt(cds_mutation_matrix)
      heatmap_data = heatmap_data %>% mutate(Classification = sapply(Var2, classify_mutation))
      
      # Create the heatmap
      plot_ly(data = heatmap_data, x = ~Classification, y = ~Var1, z = ~value, type = "heatmap",  colors = "Blues") %>%
        layout(
          title = "Coding Site Mutation",
          xaxis = list(title = " ",
                       tickangle = 45),  # Custom X-axis label
          yaxis = list(title = " "))
      
    })
    
    ##-primary histology
    output$histology_heatmap <- renderPlotly({
      cosmic_data <- cosmic_df()
      
      # Count the frequency of each MutationCDS for each gene
      mutation_counts <- cosmic_data %>%
        as.data.frame() %>%
        group_by(GeneName, PrimaryHistology) %>%
        summarise(n = n(), .groups = 'drop') %>%
        pivot_wider(names_from = PrimaryHistology, values_from = n, values_fill = list(n = 0))
      
      # Convert to a matrix for the heatmap
      mutation_counts_matrix <- as.matrix(mutation_counts[-1])
      rownames(mutation_counts_matrix) <- mutation_counts$GeneName
      heatmap_data2 <- reshape2::melt(mutation_counts_matrix)

      # Create the heatmap
      plot_ly(data = heatmap_data2, x = ~Var2, y = ~Var1, z = ~value, type = "heatmap",  colors = "Blues") %>%
        layout(
          title = "Primary Histology",
          xaxis = list(title = " ",
                       tickangle = 45),  # Custom X-axis label
          yaxis = list(title = " "))
       
    })
    
    add_log_message("Data processing and plot generation completed")
    Sys.sleep(3)
    
    observeEvent(input$viewMarkerData, {
      shinyjs::toggle("markerDataContent")
    })
    
    observeEvent(input$viewMutationData, {
      shinyjs::toggle("mutationDataContent")
    })
  
  })
}


##----------------------------------------------------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------------------------------------------------
