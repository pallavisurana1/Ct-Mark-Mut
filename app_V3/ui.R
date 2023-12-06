
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

pacman::p_load(httr, jsonlite, dplyr, readr, curl, data.table, shiny, shinyjs, DT, tidyr, janitor, scales, ggplot2, pheatmap, plotly)

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
ui <- fluidPage(
  useShinyjs(),  # Initialize shinyjs
  tags$head(tags$style(
    HTML(".well { background-color: #f7f7f7; border: none; padding: 20px; border-radius: 5px; }"),
    HTML("#loadingPopup { position: absolute; top: 10px; right: 10px; background-color: #ffffff; padding: 10px; border: 1px solid #ccc; border-radius: 5px; }"))
  ),
  
  titlePanel("Ct-Mark-Mut", windowTitle = "Ct-Mark-Mut"),
  
  # Horizontal layout for Selection Parameters
  wellPanel(
    fluidRow(
      column(4,
             selectInput("organ1", "Organ:", choices = unique(cto_all$organ), selected = NULL)
      ),
      column(4,
             selectInput("cell_type1", "Cell Type:", choices = NULL, selected = NULL)
      ),
      column(4,
             checkboxGroupInput("speciesInput", "Species:", choices = list("Hs" = "Hs", "Mm" = "Mm"), selected = c("Hs", "Mm"))
      )
    ),
    fluidRow(
      column(12,
             actionButton("run", "Retrieve Mutations", class = "btn-primary", style = "width: 100%;")
      )
    )
  ),
  
  # Heatmap and Data Display
  fluidRow(
    column(4, plotlyOutput("mutationAA_heatmap", height = "300px")),
    column(4, plotlyOutput("mutationCDS_heatmap", height = "300px")),
    column(4, plotlyOutput("histology_heatmap", height = "300px"))
  ),
  
  tags$hr(),
  
  fluidRow(
    column(6, 
           actionButton("viewMarkerData", "View Marker Data", class = "btn-primary", style = "width: 100%;"),
           hidden(
             div(
               id = "markerDataContent",
               DTOutput("markerTable"),
               div(
                 style = "display: flex; justify-content: center;",
                 downloadButton("downloadMarkerTable", "Download Marker Data")
               )
             )
           )
    ),
    column(6,
           actionButton("viewMutationData", "View Mutation Data", class = "btn-primary", style = "width: 100%;"),
           hidden(
             div(
               id = "mutationDataContent",
               DTOutput("mutationTable"),
               div(
                 style = "display: flex; justify-content: center;",
                 downloadButton("downloadMutationTable", "Download Mutation Data")
               )
             )
           )
    ),
    # Loading popup (initially hidden)
    div(
      id = "loadingPopup",
      class = "well",
      style = "display: none;",
      "Loading data, please wait..."
    )
  ),
  
  # Log messages display
  div(
    id = "logMessages",
    class = "well",
    style = "display: none;",
    tags$ul()
  )
)




##----------------------------------------------------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------------------------------------------------



