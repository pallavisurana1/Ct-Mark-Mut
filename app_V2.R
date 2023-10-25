
##----------------------------------------------------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------------------------------------------------

#-Idea of the shiny app
#-1. Give the user options to choose cell type of interest
#-2. Choose a species
#-3. User sees a list of marker genes of that cell type of interest in table from panglaoDB
#-4. fetch mutations landscape for marker genes and summary plot or table

##----------------------------------------------------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------------------------------------------------

pacman::p_load(httr, jsonlite, dplyr, readr, curl, data.table, shiny, DT, tidyr, janitor, scales, ggplot2, plotly)

##----------------------------------------------------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------------------------------------------------

url <- "https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz"

# Read data directly without downloading
panglaoDB_csv <- read_tsv(url) %>%
  clean_names() %>%
  unique() %>%
  drop_na(c(organ, cell_type, official_gene_symbol))

cto_all <- panglaoDB_csv %>%
  select(organ, cell_type) %>%
  unique()

##----------------------------------------------------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------------------------------------------------


# Define UI
ui <- fluidPage(
  
  titlePanel(
    h1("Ct-Mark-Mut"),
    tags$small("Mutation Data for Cell Type Marker Genes")
  ),
  
  sidebarLayout(
    
    sidebarPanel(
      
      wellPanel(
        tags$h4("Select Parameters"),
        
        # Organ Selection
        selectInput(
          inputId = "organ1",
          label = "Organ:",
          choices = unique(cto_all$organ),
          selected = NULL
        ),
        
        # Cell Type Selection
        selectInput(
          inputId = "cell_type1",
          label = "Cell Type:",
          choices = NULL,
          selected = NULL
        ),
        
        # Organism of Interest Selection
        checkboxGroupInput(
          "speciesInput", 
          label = "Species:", 
          choices = list("Hs" = "Hs", "Mm" = "Mm"),
          selected = c("Hs","Mm")
        ),
        
        tags$br(),
        
        # Run Button
        actionButton(
          "run", 
          "Retrieve Mutations"
        )
      )
      
    ),
    
    mainPanel(
      
      # Marker Table Display
      tags$h3("Marker Data (PanglaoDB)"),
      DTOutput("markerTable"),
      
      tags$br(),
      
      downloadButton(
        "downloadMarkerTable", 
        label = "Download Marker Data"
      ),
      
      tags$br(),
      tags$hr(),
      tags$br(),
      
      # Gene Selection
      tags$h3("Select Gene(s) from Marker Table"),
      selectInput(
        "selectedGenes", 
        choices = NULL, 
        multiple = TRUE
      ),
      
      tags$br(),
      tags$hr(),
      tags$br(),
      
      # Mutation Table Display
      tags$h3("Mutation Data (COSMIC)"),
      DTOutput("mutationTable"),
      
      tags$br(),
      
      downloadButton(
        "downloadMutationTable", 
        label = "Download Mutation Data"
      ),
      
      tags$br(),
      tags$hr(),
      tags$br(),
      
      # Plot Display
      tags$h3("Site Burden Plot"),
      plotlyOutput("site_burden_plot"),
      tags$br(),
      tags$hr()
    )
  )
)


##----------------------------------------------------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------------------------------------------------


# Define server 
server <- function(input, output, session) {
  
  ##-pick the organ and cell type of interest from user
  observeEvent(input$organ1, {
    choices <- cto_all$cell_type[cto_all$organ == input$organ1]
    updateSelectInput(session, "cell_type1", choices = choices)
  })
  
  
  observeEvent(input$run, {
    
    ##-get cell type marker information from panglaoDb tsv file
    markers <- reactive({
      x = panglaoDB_csv  %>% subset(organ == input$organ1) %>% subset(cell_type == input$cell_type1)
      
      ##-take inot account the species that the user picks for the cell types markers information
      if("Mm" %in% input$speciesInput && "Hs" %in% input$speciesInput){
        x <- x %>% subset(species %in% c("Mm Hs"))
      }else if("Mm" %in% input$speciesInput){
        x <- x %>% subset(species == "Mm")
      }else if("Hs" %in% input$speciesInput){
        x <- x %>% subset(species == "Hs")
      }
      
      x
    })
    
    # Observe changes in markerTable and update available genes for selection
    observe({
      available_genes <- markers()$official_gene_symbol
      updateSelectInput(session, "selectedGenes", choices = available_genes)
    })
    
    
    ##-genes to query from COSMIC
    gene_query1 <- reactive({
      input$selectedGenes
    })
    
    
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
          page <- page + 1
        }
      }
      return(all_data)
    }
    
    
    ##-COSMIC mutation data for cell type markers
    cosmic_df <- reactive({
      
      fetched_df = fetch_data(gene_query1(), page_limit = 10)
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
    
    
    
    ##-mutation burden plot
    output$site_burden_plot <- renderPlotly({
      
      ##-plot for primary site
      p <- ggplot(cosmic_df(), aes(x = reorder(PrimarySite, -mut_per_mb_Site), y = mut_per_mb_Site)) +
        geom_col(fill = "steelblue") +
        labs(title = "Mutation Burden by Primary Site",
             x = "Primary Site",
             y = "Mutations per Megabase (-log10)") +
        theme(plot.title = element_text(size = 18),
              axis.title = element_text(size = 14),
              panel.background = element_blank(),
              panel.grid.major = element_line(colour = "gray90"),
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"),
              axis.ticks = element_line(colour = "black")) 
      theme_classic()
      
      ##-convert to plotly
      ggplotly(p) %>% layout(xaxis = list(tickangle = 300))  
      
    })
    
  })
}


##----------------------------------------------------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------------------------------------------------


##-Run the app
shinyApp(ui, server)


##----------------------------------------------------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------------------------------------------------



