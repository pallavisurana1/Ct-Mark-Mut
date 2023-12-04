
# ui.R

##----------------------------------------------------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------------------------------------------------

#-Idea of the shiny app
#-1. Give the user options to choose cell type of interest
#-2. Choose a species
#-3. User sees a list of marker genes of that cell type of interest in table from panglaoDB
#-4. fetch mutations landscape for marker genes and summary plot or table

##----------------------------------------------------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------------------------------------------------

pacman::p_load(httr, jsonlite, dplyr, readr, curl, data.table, shiny, DT, tidyr, janitor, scales, ggplot2, pheatmap, plotly)

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

# map tissues across panglaoDB and COSMIC v4
tissue_mapping <- c(
  Pancreas = "pancreas",
  'Connective tissue' = "connective_tissue",
  Brain = "central_nervous_system",
  Lungs = "lung",
  'Smooth muscle' = "smooth_muscle",
  'Immune system' = "immune_system",
  Epithelium = "epithelium",
  Heart = "heart",
  Liver = "liver",
  'Adrenal glands' = "adrenal_gland",
  'GI tract' = "gastrointestinal_tract",
  Reproductive = "genital_tract",
  Kidney = "kidney",
  Zygote = "ns",  # DO NOT MAP
  Vasculature = "soft_tissue",  #
  Embryo = "ns",  # 
  Blood = "haematopoietic_and_lymphoid",
  Thyroid = "thyroid",
  Bone = "bone",
  Skin = "skin",
  'Mammary gland' = "breast",
  Eye = "eye",
  'Skeletal muscle' = "ns",  # 
  'Olfactory system' = "ns",  # 
  'Parathyroid glands' = "parathyroid",
  'Oral cavity' = "upper_aerodigestive_tract",
  Thymus = "thymus",
  Placenta = "placenta",
  'Urinary bladder' = "urinary_tract"
)

# Function to map tissue if a match is found
map_tissue_if_match <- function(organ_name) {
  if (organ_name %in% names(tissue_mapping)) {
    return(tissue_mapping[organ_name])
  } else {
    return(NA)
  }
}

cto_all <- cto_all %>%
  mutate(cosmic_tissue = sapply(organ, map_tissue_if_match))

##----------------------------------------------------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------------------------------------------------


# Define UI
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
      
      # Amino Acid Mutations Frequency Plot Display
      tags$h3("Amino Acid Mutations Frequency"),
      plotlyOutput("mutationAA_plot"),
      tags$br(),
      tags$hr(),
      tags$br(),
      
      # Coding DNA Sequence Mutations Frequency Plot Display
      tags$h3("Coding DNA Sequence Mutations Frequency"),
      plotlyOutput("mutationCDS_plot"),
      tags$br(),
      tags$hr()
    )
  )
)


##----------------------------------------------------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------------------------------------------------



