#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(leaflet)
library(viridis)

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
 
 output$economic_es <- renderLeaflet({
     leaflet() %>%
         addTiles() %>%
         addRasterImage(economic_es_mean_raster, colors = pal, opacity = 0.5) %>%
         addLegend(values = values(economic_es_mean_raster), pal = pal, title = "Amount (Normalized)")
 })
 
 output$ecological_es <- renderLeaflet({
     leaflet() %>%
         addTiles() %>%
         addRasterImage(ecological_es_mean_raster, colors = pal) %>%
         addLegend(values = values(ecological_es_mean_raster), pal = pal, title = "Amount (Normalized)")
 })
 
 output$social_es <- renderLeaflet({
     leaflet() %>%
         addTiles() %>%
         addRasterImage(social_es_mean_raster, colors = pal) %>%
         addLegend(values = values(social_es_mean_raster), pal = pal, title = "Amount (Normalized)")
 })
 
 output$economic_priori <- renderTmap({
     hotspot_economic_priori_tm
 })
 
 output$ecological_priori <- renderMapview({
     mapview::mapview(ecological_prioritized_raster, query.type = 'mousemove', 
                      query.digits = 2, trim = TRUE, na.opacity = 1)
 }
    
 )

})
