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
   #     sync(m1, m2, m3, m4)
 output$economic_es <- renderMapview({
     # leaflet() %>%
     #     addTiles() %>%
     #     addRasterImage(economic_es_mean_raster, colors = pal, opacity = 0.5) %>%
     #     addLegend(values = values(economic_es_mean_raster), pal = pal, title = "Amount (Normalized)")
         
   economic_es_mean_clip1_raster <- raster(economic_es_mean_clip1_filepath)
         mapview::mapview(economic_es_mean_clip1_raster, layer.name = "Economic Ecosystem Services")
 })
 
 output$ecological_es <- renderMapview({
     # leaflet() %>%
     #     addTiles() %>%
     #     addRasterImage(ecological_es_mean_raster, colors = pal) %>%
     #     addLegend(values = values(ecological_es_mean_raster), pal = pal, title = "Amount (Normalized)")
     #     
    ecological_es_mean_raster <- raster(ecological_es_mean_file_path)     
    mapview::mapview(ecological_es_mean_raster, layer.name = "Ecological Ecosystem Services")
 })
 
 output$social_es <- renderMapview({
     # leaflet() %>%
     #     addTiles() %>%
     #     addRasterImage(social_es_mean_raster, colors = pal) %>%
     #     addLegend(values = values(social_es_mean_raster), pal = pal, title = "Amount (Normalized)")
     #     
          social_es_mean_raster <- raster(social_es_mean_filepath)
         mapview::mapview(social_es_mean_raster, layer.name = "Social Ecosystem Services")
 })
 
 output$economic_priori <- renderMapview({
   economic_prioritized_raster <- raster::raster(economic_prioritized_filepath)
   norm_economic_prioritized_raster <- norm_function(economic_prioritized_raster)
  
   mapview::mapview(norm_economic_prioritized_raster, layer.name = "Economic prioritization")
 })
 
 output$ecological_priori <- renderMapview({
   ecological_prioritized_raster <- raster::raster(ecological_prioritized_filepath)
   norm_ecological_prioritized_raster <- norm_function(ecological_prioritized_raster)
   
     mapview::mapview(norm_ecological_prioritized_raster, layer.name = "Ecological prioritization", 
                      query.type = 'mousemove', 
                      query.digits = 2)
 })
 
 output$social_priori <- renderMapview({
   social_prioritized_raster <- raster::raster(social_prioritized_filepath)
   norm_social_prioritized_raster <- norm_function(social_prioritized_raster)
   
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
   
   
   raster(ecological_priori_hotspot_filepath) -> ecological_priori_hotspot_raster
   raster(economic_priori_hotspot_filepath) -> economic_priori_hotspot_raster
   raster(social_priori_hotspot_filepath) -> social_priori_hotspot_raster
   raster(equal_priori_hotspot_filepath) -> equal_priori_hotspot_raster
   
   m1 <- mapview(equal_priori_hotspot_raster, query.type = 'mousemove', 
                 query.digits = 2, legend.opacity = 0, layer.name = "Equal prioritization") 
   
   m2 <- mapview(economic_priori_hotspot_raster, query.type = 'mousemove',
                 query.digits = 2, layer.name = "Economic prioritization")
   
   m3 <- mapview(ecological_priori_hotspot_raster, query.type = 'mousemove',
                 query.digits = 2, layer.name = "Ecological prioritization")
   
   m4 <- mapview(social_priori_hotspot_raster, query.type = 'mousemove',
                 query.digits = 2, layer.name = "Social prioritization", na.alpha = 0.1) # + paste0("Social prioritization hotspot")
   
   sync(list(m1, m2, m3, m4)) 
   
       #  sync(list(m1, m2, m3, m4)) 
 })
 
 ###### Geocomputations
 output$correlation_es <- renderPlot({
   
   ## read raster files
   #ecological_es_sum_raster <- raster(ecological_es_sum_file_path)
    ecological_es_mean_raster <- raster(ecological_es_mean_file_path)
   #economic_es_mean_raster <- raster(economic_es_mean_filepath)
   #economic_es_mean_resam_raster <- raster(economic_es_mean_resam_filepath)
    social_es_mean_raster <- raster(social_es_mean_filepath)
    economic_es_mean_clip1_raster <- raster(economic_es_mean_clip1_filepath)
   
   # resam_equal_weight_raster <- resample(equal_weight_raster, norm_ecological_es_mean_raster)
   # norm_resam_equal_weight_raster <- norm_function(resam_equal_weight_raster)
   
   resam_social_es_mean_raster <- resample(social_es_mean_raster, ecological_es_mean_raster)
   resam_economic_es_mean_raster <- resample(economic_es_mean_clip1_raster, ecological_es_mean_raster)
   
   norm_resam_social_es_mean_raster <- norm_function(resam_social_es_mean_raster)
   norm_resam_economic_es_mean_raster <- norm_function(resam_economic_es_mean_raster)
    norm_ecological_es_mean_raster <- norm_function(ecological_es_mean_raster)
   
    all_es_mean_stack <- raster::stack(norm_resam_economic_es_mean_raster, norm_ecological_es_mean_raster, norm_resam_social_es_mean_raster)
    
   raster::pairs(all_es_mean_stack)
 })
 
 output$economic_spa_autocorr <- renderValueBox({
        # h4("Spatial autocorrelation")
   economic_es_mean_clip1_raster <- raster(economic_es_mean_clip1_filepath)
#   norm_economic_es_mean_clip1_raster <- norm_function(economic_es_mean_clip1_raster)
   
   economic_moran <- Moran(economic_es_mean_clip1_raster)
   
         valueBox(value = h2(formattable::comma(economic_moran, digits = 3), style = "background: green"), 
                  subtitle = "Economic Ecosystem Services" 
                  # href = "https://agbleze.github.io/Profile/index.html", 
                  # HTML("style = color: blue; background: green")
         )
 })
 
 output$ecological_spa_autocorr <- renderValueBox({
   ecological_es_mean_raster <- raster(ecological_es_mean_file_path)
   ecological_moran <- Moran(ecological_es_mean_raster)
   
         valueBox(value = h2(formattable::comma(ecological_moran, digits = 3), style = "background: green"), 
                  subtitle = "Ecological Ecosystem Services")
 })
 
 output$social_spa_autocorr <- renderValueBox({
   
   social_es_mean_raster <- raster(social_es_mean_filepath)
   
   social_moran <- Moran(social_es_mean_raster)
   
         valueBox(value = h2(formattable::comma(social_moran, digits = 3), style = "background: green"),
                  subtitle = "Social Ecosystem Services"
                  )
 })

})

