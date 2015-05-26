library(shiny)
library(RCurl)

# runbefore --------------------
eval( expr = parse( text = getURL(
    "https://raw.githubusercontent.com/oganm/toSource/master/ogbox.R",
    ssl.verifypeer=FALSE) ))
sourceGithub(oganm,masterOfCellTypes,runVars)
groupNames = 'PyramidalDeep'
mouseDes = read.design('Data/meltedDesign.tsv')

regions =
    trimNAs(
        trimElement(
            unique(
                unlist(
                    strsplit(as.character(mouseDes[,regionNames]),',')))
            ,c('ALL','All','all','Cerebrum'))) #S pecial names

# UI ------------------------

shinyUI(fluidPage(
    titlePanel("cell type expression"),
    sidebarLayout(
        sidebarPanel(
            p('Wait until the plot renders then enter a gene symbol.'),
            textInput(inputId = 'geneSearch',value = 'Ogn',
                      label = 'Select Gene'),
            
            selectInput(inputId = "regionChoice",
                        label= 'Select region',
                        choices = c(regions,'All','.messy details')),
            
            checkboxInput(inputId = 'jitterBox',
                          label = 'Jitter?', 
                          value = F)
        ),
        mainPanel(
            plotOutput('expressionPlot')
        )
    )
))