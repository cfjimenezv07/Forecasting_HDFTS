source("mapRender.R")
source("renderForecastTable.R")
source("rainbowRender.R")
source("intervalMeasurements.R")

get.state <- function(input) {
  if(input$countrySelector == "USA") {
    state <- input$stateSelector
  } else if(input$countrySelector == "France") {
    state <- input$departmentSelector
  } else if(input$countrySelector == "Australia") {
    state <- input$regionSelector
  } else {
    state <- input$prefectureSelector
  }
  state
}

SERVER <- function(input, output) {
  output$map <- renderLeaflet({
    render.map(input$countrySelector, get.state(input))
  })
  output$metricsTable <- DT::renderDataTable({
    ltable <- point.forecast.table(input$countrySelector, 
                                   input$genderSelector, "RMSPE", get.state(input))
    rtable <- point.forecast.table(input$countrySelector, 
                                   input$genderSelector, "MAPE", get.state(input), F)
    out <- cbind(ltable, rtable)
    sketch.upper <- htmltools::withTags(table(
      class = 'display',
      thead(
        tr(
          th(colspan = 1, ''),
          th(colspan = 2, 'RMSPE'),
          th(colspan = 2, 'MAPE')
        ),
        tr(
          lapply(c('h', rep(c('FMP', 'FANOVA'), 2)), th)
        )
      )
    ))
    DT::datatable(round(out, 4),
                  options = list(pageLength = 5), rownames = F, 
                  container = sketch.upper)
  })
  output$rainbowPlot <- renderPlot({
    rainbow.generator(input$countrySelector, get.state(input), input$genderSelector)
  })
  output$summaryMetricsTable <- DT::renderDataTable({
    sketch.lower <- htmltools::withTags(table(
      class = 'display',
      thead(
        tr(
          th(colspan = 1, ''),
          th(colspan = 2, 'RMSPE'),
          th(colspan = 2, 'MAPE')
        ),
        tr(
          lapply(c('', rep(c('FMP', 'FANOVA'), 2)), th)
        )
      )
    ))
    ltable <- point.forecast.table(input$countrySelector, 
                                   input$genderSelector, "RMSPE", get.state(input))
    rtable <- point.forecast.table(input$countrySelector, 
                                   input$genderSelector, "MAPE", get.state(input), F)
    out <- cbind(ltable, rtable)[, -1]
    outmeans <- colMeans(out)
    outsd <- apply(out, 2, sd)
    outsummary <- rbind(outmeans, outsd)
    rownames(outsummary) <- c("Mean", "Sd")
    DT::datatable(round(outsummary, 4),
                  options = list(pageLength = 5, bInfo = F, paging = F),  
                  container = sketch.lower)
    
  })
  output$measurementsTable <- DT::renderDataTable({
    sketch.measurements <- htmltools::withTags(table(
      class = 'display',
      thead(
        tr(
          th(colspan = 1, ''),
          th(colspan = 2, 'FMP'),
          th(colspan = 2, 'FANOVA')
        ),
        tr(
          lapply(c('', rep(c('80%', '95%'), 2)), th)
        )
      )
    ))
    DT::datatable(measurementsTableGenerator(get.state(input), input$genderSelector),
                  options = list(bInfo = F, paging = F), 
                  container = sketch.measurements)
  })
}
