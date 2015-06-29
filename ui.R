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
            textOutput(outputId = 'didYouMean'),
            
            selectInput(inputId = "regionChoice",
                        label= 'Select region',
                        choices = c(regions,'All','.messy details')),
            
            checkboxInput(inputId = 'jitterBox',
                          label = 'Jitter?', 
                          value = F),
            sliderInput(inputId = 'pointSize',
                        label = 'Point size',
                        min = 4,
                        max = 20,
                        value = 7),
            sliderInput(inputId = 'ySize',
                        label = 'y text size',
                        min = 0,
                        max = 36,
                        value = 14),
            sliderInput(inputId = 'yTitleSize',
                        label = 'y title size',
                        min = 0,
                        max = 36,
                        value = 20),
            sliderInput(inputId = 'xSize',
                        label = 'x text size',
                        min = 0,
                        max = 36,
                        value = 20),
            checkboxInput(inputId = 'color',
                          label = 'Color?', 
                          value = T),
            textInput(inputId = 'additionalGG',value = '',
                      label = 'More layers')
            
        ),
        mainPanel(
            plotOutput('expressionPlot')
        )
    )
))
