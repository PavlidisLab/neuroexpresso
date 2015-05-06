library(shiny)


shinyUI(fluidPage(
    titlePanel("cell type expression"),
    sidebarLayout(
        sidebarPanel(
            textInput(inputId = 'geneSelect',value = '',
                      label = 'Select Gene')
        ),
        mainPanel(
            plotOutput('expressionPlot')
        )
    )
))