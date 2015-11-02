library(shiny)
library(RCurl)
library(ggvis)
# runbefore --------------------
library(ogbox)
# sourceGithub(oganm,masterOfCellTypes,runVars)
groupNames = 'PyramidalDeep'
regionNames = 'Region'
mouseDes = read.design('Data/meltedDesign.tsv')

regions =
    trimNAs(
        trimElement(
            unique(
                unlist(
                    strsplit(as.character(mouseDes[,regionNames]),',')))
            ,c('ALL','All','all','Cerebrum'))) #S pecial names

# user ID system
inputUserid <- function(inputId, value='') {
    #   print(paste(inputId, "=", value))
    tagList(
        singleton(tags$head(tags$script(src = "js/md5.js", type='text/javascript'))),
        singleton(tags$head(tags$script(src = "js/shinyBindings.js", type='text/javascript'))),
        tags$body(onload="setvalues()"),
        tags$input(id = inputId, class = "userid", value=as.character(value), type="text", style="display:none;")
    )
    
}

inputIp <- function(inputId, value=''){
    
    tagList(
        singleton(tags$head(tags$script(src = "js/md5.js", type='text/javascript'))),
        singleton(tags$head(tags$script(src = "js/shinyBindings.js", type='text/javascript'))),
        tags$body(onload="setvalues()"),
        tags$input(id = inputId, class = "ipaddr", value=as.character(value), type="text", style="display:none;")
    )
    
}

# UI ------------------------

shinyUI(fluidPage(
    titlePanel("Expression in brain cell types"),
    sidebarLayout(
        sidebarPanel(
            tags$head(tags$script('$(function () { $("#expressionPlot").click(function(e){ $("#ggvis-tooltip").hide(); }); })')),
            inputIp("ipid"),
            inputUserid("fingerprint"),
            h1('About'),
            p('This application aims to make it easier to visualize gene expression in mouse brain cell types. The data here is compiled for a project aiming to select cell type specific genes in brain.'),
            p('It is compiled by combining data from', 
                    a(href="http://www.chibi.ubc.ca/Gemma/arrays/showArrayDesign.html?id=7", 'GPL339'),'and',
                    a(href="http://www.chibi.ubc.ca/Gemma/arrays/showArrayDesign.html?id=3", 'GPL1261')),
            p('To see genes that are only available for GPL1261, choose that platform below. This will remove some of the samples'),
            p('Click on data points to see their sources'),
            a(href="https://github.com/oganm/cellTypeExpression", 'Source Code'),
            br(),
            a(href="https://github.com/oganm/brainCellTypeSpecificGenes", 'Project Page'),
            br(),
            p('Wait until the plot renders then enter a gene symbol.'),
            textInput(inputId = 'geneSearch',value = 'Ogn',
                      label = 'Select Gene'),
            textOutput(outputId = 'didYouMean'),
            
            selectInput(inputId = "regionChoice",
                        label= 'Select region',
                        # choices = c(regions,'All','.messy details')),
                        choices = c(regions,'All')),
            selectInput(inputId = 'platform',
                        label = 'Select Platform',
                        choices = c('GPL339','GPL1261')),
            
           # checkboxInput(inputId = 'privacyBox',
        #                  label = "We keep an anonymized log of valid searches. Uncheck this box before you leave if you are feeling paranoid so we won't.", 
            #              value = T),
        
            uiOutput('expressionUI'),
            checkboxInput(inputId = 'color',
                          label = 'Color?', 
                          value = T)
            #textInput(inputId = 'additionalGG',value = '',
            #          label = 'More layers')
            
        ),
        mainPanel(
            ggvisOutput('expressionPlot')
            #plotOutput('expressionPlot')
        )
    )
))
