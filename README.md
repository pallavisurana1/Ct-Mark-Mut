## Ct-Mark-Mut

Ct-Mark-Mut is a Shiny app that allows users to explore mutations in marker genes for specific cell types in different species. The app provides access to a list of marker genes for a given cell type, fetched from the PanglaoDB database, and then allows the user to query the COSMIC mutation database for mutations in the selected marker genes.

### Usage

To use the app, follow these steps:

1.  Choose the organ of interest from the drop-down list.

2.  Choose the cell type corresponding to the selected organ.

3.  Choose one or more species from the check-box options to display mutations for the selected species.

4.  Click the "Retrieve mutations" button to fetch mutations in the selected marker genes.

6.  The app will display plots for Amino Acid and CDS region mutations in the data followed by optional viewing and downloading of tables.

7.  You can also download the tables by clicking the "Download Table" button once "View..." is clicked on.


### Requirements

The app requires the following R packages:

-   httr
-   jsonlite
-   dplyr
-   data.table
-   readr
-   shiny
-   shinyjs
-   DT
-   tidyr
-   janitor
-   scales
-   ggplot2
-   pheatmap
-   plotly

**To run the RShiny app from the Github repository, follow the steps below:**

1.  Clone the repository to your local machine using the following command in your terminal or command prompt:

    **`git clone https://github.com/PallaviSurana1/Ct-Mark-Mut.git`**

2.  Make sure you have R and RStudio installed on your machine.

3.  Open RStudio and set your working directory to the directory where you cloned the repository using the **`setwd()`** function. For example, if you cloned the repository to your desktop, you can set your working directory as follows:

    **`setwd("~/path/to/download/Ct-Mark-Mut")`**

4.  Install the required packages by running the following command in your R console:

    **`install.packages(c("shiny", "dplyr"))`**

5.  Load the RShiny library by running the following command in your R console:

    **`library(shiny)`**

6.  Run the app by running the following command in your R console:

    **`runApp("Ct-Mark-Mut")`**

    This will open the app in a new window in your default web browser.

Note: Before running the app, make sure to check the **`readme.md`** file in the repository to ensure you have the necessary data files in the correct directory.

### Functionality

There are some input widgets defined in the sidebarPanel() of the app.

1. `selectInput()` is used to select an organ of interest from a list of organs.
2. `selectInput()` is used to select a cell type of interest from a list of cell types that correspond to the selected organ.
3. `checkboxGroupInput()` is used to select the species for which the mutation data is needed.
4. `actionButton()` is used to initiate the process to fetch mutations in marker genes.

There are several output widgets defined in the `ui.R` code.

1. `plotlyOutput()` is used to display the Amino Acid Mutation heatmap.
2. `plotlyOutput()` is used to display the CDS Region Mutation heatmap.
3. `plotlyOutput()` is used to display the Primary Histology heatmap.
4. `actionButton()` is used to view the marker gene data.
5. `actionButton()` is used to view the mutation data.
6. Loading pop-up messages are displayed during data processing.
7. Log messages are displayed in a well panel.

### Acknowledgments

This app is based on data from the following sources:

The PanglaoDB database The COSMIC v4 database
