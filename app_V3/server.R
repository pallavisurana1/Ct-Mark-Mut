
# server.R

##----------------------------------------------------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------------------------------------------------


# Define server 
server <- function(input, output, session) {
  
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
    
    output$mutationAA_plot <- renderPlotly({
      cosmic_data <- cosmic_df()  
      
      aa_mutation_counts <- cosmic_data %>%
        group_by(MutationAA, GeneName) %>%
        summarise(Frequency = n()) %>%
        ungroup() %>%
        arrange(desc(Frequency))
      
      p_aa <- ggplot(aa_mutation_counts, aes(x = MutationAA, y = Frequency, fill = GeneName)) +
        geom_bar(stat = "identity") +
        labs(title = "Frequency of Amino Acid Mutations",
             x = "Amino Acid Mutation",
             y = "Frequency") +
        scale_fill_brewer(palette = "Set1") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
              legend.title = element_blank())
      
      ggplotly(p_aa)
    })
    
    
    output$mutationCDS_plot <- renderPlotly({
      cosmic_data <- cosmic_df()  
      
      cds_mutation_counts <- cosmic_data %>%
        group_by(MutationCDS, GeneName) %>%
        summarise(Frequency = n()) %>%
        ungroup() %>%
        arrange(desc(Frequency))
      
      p_cds <- ggplot(cds_mutation_counts, aes(x = MutationCDS, y = Frequency, fill = GeneName)) +
        geom_bar(stat = "identity") +
        labs(title = "Frequency of Coding DNA Sequence Mutations",
             x = "Coding DNA Sequence Mutation",
             y = "Frequency") +
        scale_fill_brewer(palette = "Set1") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
              legend.title = element_blank())
      
      ggplotly(p_cds)
    })
    
  
    
  })
}


##----------------------------------------------------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------------------------------------------------
