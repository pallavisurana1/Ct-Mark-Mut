## Ct.Mark.Mut

Ct-Mark-Mut is a Shiny app that allows users to explore mutations in marker genes for specific cell types in different species. The app provides access to a list of marker genes for a given cell type, fetched from the PanglaoDB database, and then allows the user to query the COSMIC mutation database for mutations in the selected marker genes.

### Usage

To use the app, follow these steps:

1.  Choose the organ of interest from the drop-down list.

2.  Choose the cell type corresponding to the selected organ.

3.  Choose one or more species from the check-box options to display mutations for the selected species.

4.  Set the number of pages you want to fetch from COSMIC. The number of mutations displayed will depend on the page limit.

5.  Click the "Show me the mutations in the cell type markers" button to fetch mutations in the selected marker genes.

6.  The app will display a table of the selected marker genes from PanglaoDB, followed by a table of mutations in the marker genes from COSMIC.

7.  You can also download the tables by clicking the "Download Table" button.

8.  Finally, a summary plot of the mutation burden is displayed at the bottom of the page.

### Requirements

The app requires the following R packages:

-   httr

-   jsonlite

-   dplyr

-   data.table

-   shiny

-   DT

-   tidyr

-   janitor

-   scales

-   ggplot2

-   plotly

### Acknowledgments

This app is based on data from the following sources:

The PanglaoDB database The COSMIC database

### Contributors

The app was developed by Pallavi Surana.
