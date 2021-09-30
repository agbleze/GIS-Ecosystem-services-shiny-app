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
        sync(m1, m2, m3, m4)
 output$economic_es <- renderMapview({
     # leaflet() %>%
     #     addTiles() %>%
     #     addRasterImage(economic_es_mean_raster, colors = pal, opacity = 0.5) %>%
     #     addLegend(values = values(economic_es_mean_raster), pal = pal, title = "Amount (Normalized)")
         
         mapview::mapview(norm_economic_es_mean_clip1_raster, layer.name = "Economic Ecosystem Services")
 })
 
 output$ecological_es <- renderMapview({
     # leaflet() %>%
     #     addTiles() %>%
     #     addRasterImage(ecological_es_mean_raster, colors = pal) %>%
     #     addLegend(values = values(ecological_es_mean_raster), pal = pal, title = "Amount (Normalized)")
     #     
         mapview::mapview(ecological_es_mean_raster, layer.name = "Ecological Ecosystem Services")
 })
 
 output$social_es <- renderMapview({
     # leaflet() %>%
     #     addTiles() %>%
     #     addRasterImage(social_es_mean_raster, colors = pal) %>%
     #     addLegend(values = values(social_es_mean_raster), pal = pal, title = "Amount (Normalized)")
     #     
         mapview::mapview(social_es_mean_raster, layer.name = "Social Ecosystem Services")
 })
 
 output$economic_priori <- renderMapview({
     mapview::mapview(norm_economic_prioritized_raster, layer.name = "Economic prioritization")
 })
 
 output$ecological_priori <- renderMapview({
     mapview::mapview(norm_ecological_prioritized_raster, layer.name = "Ecological prioritization", 
                      query.type = 'mousemove', 
                      query.digits = 2)
 })
 
 output$social_priori <- renderMapview({
     mapview::mapView(norm_social_prioritized_raster, layer.name = "Social prioritization") + 
                 mapview::mapview(fran_adm1_shp, zcol = "NAME_1", alpha.regions = 0, legend = FALSE,
                                  layer.name = "Regions of France")
 })
 
 
#  output$equal_priori_hotspot <- renderMapview({
#          #m1
#          mapview(equal_priori_hotspot_raster, query.type = 'mousemove', 
#                  query.digits = 2) + mapview(economic_priori_hotspot_raster, query.type = 'mousemove',
#                                              query.digits = 2) 
# 
#          #sync(m1, m2, m3, m4)
#  })
#  
#  output$economic_priori_hotspot <- renderMapview({
#          m2
# })
#  
#  output$ecological_priori_hotspot <- renderMapview(m3)
#  output$social_priori_hotspot <- renderMapview(m4)
 
 # output$try <- renderMapview(
 #         mapview(equal_priori_hotspot_raster, query.type = 'mousemove', 
 #                       query.digits = 2) + mapview(economic_priori_hotspot_raster, query.type = 'mousemove',
 #                       query.digits = 2) + mapview(ecological_priori_hotspot_raster, query.type = 'mousemove',
 #                       query.digits = 2) + mapview(social_priori_hotspot_raster, query.type = 'mousemove',
 #                       query.digits = 2) 
 # )
 
 output$hotspot <- renderUI({
         sync(list(m1, m2, m3, m4)) 
 })
 
 ###### Geocomputations
 output$correlation_es <- renderPlot({
         raster::pairs(all_es_mean_stack)
 })
 
 output$economic_spa_autocorr <- renderValueBox({
        # h4("Spatial autocorrelation")
         valueBox(value = h2(formattable::comma(economic_moran, digits = 3), style = "background: green"), 
                  subtitle = "Economic Ecosystem Services" 
                  # href = "https://agbleze.github.io/Profile/index.html", 
                  # HTML("style = color: blue; background: green")
         )
 })
 
 output$ecological_spa_autocorr <- renderValueBox({
         valueBox(value = h2(formattable::comma(ecological_moran, digits = 3), style = "background: green"), 
                  subtitle = "Ecological Ecosystem Services")
 })
 
 output$social_spa_autocorr <- renderValueBox({
         valueBox(value = h2(formattable::comma(social_moran, digits = 3), style = "background: green"),
                  subtitle = "Social Ecosystem Services"
                  )
 })

})

