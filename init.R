# init.R
#
# Example R code to install packages if not already installed
#

my_packages = c("shiny", "shinydashboard", "shinycssloaders", "shinyanimate", "shinydashboardPlus",
                "shinyEffects", "shinybusy", "shinyBS", "magrittr", "readr", 
                "tidyverse", "ggplot2", "dplyr", "magrittr", "formattable", "fontawesome", "shinyWidgets",
                "tibble", "utils", "leaflet", "tmap", "raster", "rasterVis", "sp", "terra", "sf", "ggplot2",
                "tidyverse", "tmaptools", "mapview", "mapdeck", "terra", "RColorBrewer", "leafsync",
                "plainview", "leafem", "rasterDT", "SGAT", "spdep", "ModelMetrics", "colourvalues", 
                "googlePolylines", "spatialwidget")

install_if_missing = function(p) {
  if (p %in% rownames(installed.packages()) == FALSE) {
    install.packages(p)
  }
}
invisible(sapply(my_packages, install_if_missing))
