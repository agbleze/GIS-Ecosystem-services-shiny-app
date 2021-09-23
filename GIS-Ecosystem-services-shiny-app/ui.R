#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(leaflet)

# Define UI for application that draws a histogram
shinyUI(fluidPage(

    # Application title
    titlePanel("Ecosystem service app"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            sliderInput("bins",
                        "Number of bins:",
                        min = 1,
                        max = 50,
                        value = 30)
        ),

        # Show a plot of the generated distribution
        mainPanel(
            tabsetPanel(
                tabPanel(title = h3(strong(em("Ecosystem Services"))),
                         tabsetPanel(
                             tabPanel(title = "Economic Ecosystem Services",
                                      fluidPage(
                                          fluidRow(
                                          #    leafletOutput("leaflet_test"),
                                              mapviewOutput("economic_es")
                                          )
                                      )),
                             tabPanel(title = "Ecological Ecosystem Services",
                                      fluidPage(
                                          fluidRow(
                                            #  leafletOutput("ecological_es")
                                              mapviewOutput("ecological_es")
                                          )
                                      )),
                             tabPanel(title = "Social Ecosystem Services",
                                      fluidPage(
                                          fluidRow(
                                            #  leafletOutput("social_es")
                                              mapviewOutput("social_es")
                                          )
                                      )
                                      )
                             )),
                tabPanel(title = h3(strong(em("Prioritization"))),
                         tabsetPanel(
                             tabPanel(title = "Economic Prioritization",
                                 fluidPage(
                                     fluidRow(
                                         mapviewOutput("economic_priori")
                                     )
                             )
                             
                         ),
                         tabPanel(title = "Ecological Prioritization",
                             fluidPage(
                                 fluidRow(
                                     mapviewOutput("ecological_priori")
                                 ) 
                             )
                         ),
                         tabPanel(title = "Social Prioritization",
                                  fluidPage(
                                      fluidRow(
                                          mapviewOutput("social_priori")
                                      )
                                  )
                                  )
                         
                         )),
                tabPanel(title = h3(strong(em("Geocomputations"))),
                         tabsetPanel(
                             tabPanel(title = "Hotspot Analysis",
                                      fluidPage(
                                          uiOutput("hotspot")
                                          )
                                      #     fluidRow(
                                      #        column(width = 6, mapviewOutput("equal_priori_hotspot")),
                                      #        column(width = 6, mapviewOutput("economic_priori_hotspot"))
                                      #     )
                                      # ),
                                      # fluidPage(
                                      #     fluidRow(
                                      #         column(width = 6, mapviewOutput("ecological_priori_hotspot")),
                                      #         column(width = 6, mapviewOutput("social_priori_hotspot"))
                                      #     )
                                      # )
                                      ),
                             tabPanel(title = "Spatial Autocorrelation",
                                      fluidPage(
                                          uiOutput("try")
                                      )),
                             tabPanel(title = "Clustering"),
                             tabPanel(title = "Modelling")
                         ))
            )
        )
    )
))
