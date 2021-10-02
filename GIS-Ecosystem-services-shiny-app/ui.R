# UI of app

library(shiny)
library(leaflet)

# Define UI for application that draws a histogram
shinyUI(fluidPage(

    # Application title
    titlePanel("Ecosystem service app"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            # sliderInput("bins",
            #             "Number of bins:",
            #             min = 1,
            #             max = 50,
            #             value = 30)
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
                                              mapviewOutput("economic_es") %>%
                                                  withSpinner()
                                          )
                                      )),
                             tabPanel(title = "Ecological Ecosystem Services",
                                      fluidPage(
                                          fluidRow(
                                            #  leafletOutput("ecological_es")
                                              mapviewOutput("ecological_es") %>%
                                                  withSpinner()
                                          )
                                      )),
                             tabPanel(title = "Social Ecosystem Services",
                                      fluidPage(
                                          fluidRow(
                                            #  leafletOutput("social_es")
                                              mapviewOutput("social_es") %>%
                                                  withSpinner()
                                          )
                                      )
                                      )
                             )),
                tabPanel(title = h3(strong(em("Prioritization"))),
                         tabsetPanel(
                             tabPanel(title = "Economic Prioritization",
                                 fluidPage(
                                     fluidRow(
                                         mapviewOutput("economic_priori") %>%
                                             withSpinner()
                                     )
                             )
                             
                         ),
                         tabPanel(title = "Ecological Prioritization",
                             fluidPage(
                                 fluidRow(
                                     mapviewOutput("ecological_priori") %>%
                                         withSpinner()
                                 ) 
                             )
                         ),
                         tabPanel(title = "Social Prioritization",
                                  fluidPage(
                                      fluidRow(
                                          mapviewOutput("social_priori") %>%
                                              withSpinner()
                                      )
                                  )
                                  )
                         
                         )),
                tabPanel(title = h3(strong(em("Geocomputations"))),
                         tabsetPanel(
                             tabPanel(title = "Exploratory Data Analysis",
                                      fluidPage(
                                          fluidRow(
                                              plotOutput("correlation_es") %>%
                                                  withSpinner()
                                          )
                                      )),
                             tabPanel(title = "Hotspot Analysis",
                                      fluidPage(
                                          uiOutput("hotspot") %>%
                                              withSpinner()
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
                                          fluidRow(
                                              valueBoxOutput("economic_spa_autocorr") %>% withSpinner(),
                                              valueBoxOutput("ecological_spa_autocorr") %>% withSpinner(),
                                              valueBoxOutput("social_spa_autocorr") %>%
                                                  withSpinner()
                                              
                                          )
                                      )) #,
                             # tabPanel(title = "Clustering",
                             #          fluidPage(
                             #              fluidRow(
                             #                  
                             #              )
                             #              
                             #          )),
                             # tabPanel(title = "Modelling")
                         ))
            )
        )
    )
))
