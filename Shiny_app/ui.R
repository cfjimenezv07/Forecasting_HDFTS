library(shiny)
library(tidyverse)
library(leaflet)

Ui <- fluidPage(
  fluidRow(
    column(4,
           div(
             style = "display: flex; justify-content: center; align-items: center",
             h3(
               "Forecasting high-dimensional functional time series: Application to subnational age-specific mortality"
             )),
           wellPanel(
             selectInput(
               "countrySelector",
               "Country",
               c("USA", "France", "Japan")
             ),
             conditionalPanel(
               "input.countrySelector == 'USA'",
               selectInput("stateSelector",
                           "State",
                           sapply(readRDS("names/names_states.rds"), 
                                  function(x) x)
               )
             ),
             conditionalPanel(
               "input.countrySelector == 'France'",
               selectInput("departmentSelector",
                           "Department",
                           sapply(readRDS("names/names_departments.rds"), 
                                  function(x) x)
               )
             ),
             conditionalPanel(
               "input.countrySelector == 'Japan'",
               selectInput("prefectureSelector",
                           "Prefecture",
                           sapply(readRDS("names/names_prefectures.rds"), 
                                  function(x) x)
               )
             ),
             conditionalPanel(
               "input.countrySelector == 'Australia'",
               selectInput("regionSelector",
                           "Region",
                           sapply(readRDS("names/names_regions.rds"), 
                                  function(x) x)
               )
             ),
             selectInput("genderSelector",
                         "Gender",
                         c("Male", "Female")
             )
           ),
           plotOutput("rainbowPlot")
    ),
    column(4,
           leafletOutput("map")
    ),
    column(4,
           div(
             style = "display: flex; justify-content: center; align-items: center",
             h4("Point forecast errors")
           ),
           DT::dataTableOutput("metricsTable"),
           DT::dataTableOutput("summaryMetricsTable"),
           div(
             style = "display: flex; justify-content: center; align-items: center",
             h4("Interval forecast measurements")
           ),
           DT::dataTableOutput("measurementsTable"),
           div(
             style = "display: flex; justify-content: center; align-items: center; padding-top:5%",
             h4(
               "Authors",
               br(),
               br(),
               "Cristian F. Jiménez-Varón and Ying Sun",
               br(),
               "CEMSE Division",
               br(),
               "King Abdullah University of Science and Technology",
               br(),
               br(),
               "Han Lin Shang",
               br(),
               "Department of Actuarial Studies and Business Analytics",
               br(),
               "Macquarie University")
           )
    )
  ),
  tags$style(type = "text/css",
             ".dataTables_length, .dataTables_filter {display:none;}
             h4 {text-align:center;}
             #map {height: calc(100vh - 30px) !important;}")
)