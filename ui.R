library(shiny)
library(RCurl)
library(ggvis)
library(ogbox)
library(shinyTree)
library(shinyjs)
library(V8)

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

javaScript = "shinyjs.coinPlot = function(){
$.each($('#g64').children(),function(){$(this).attr('d','M0,4.675581A6.6755811781245455,1.6755811781245455 0 1,1 0,-3.675581A6.6755811781245455,1.6755811781245455 0 1,1 0,4.6755811781245455Z')})
}
shinyjs.changeTree = function(params){
eval(\"$('#tree').jstree(true).settings.core.data=\"+params);
$('#tree').jstree(true).refresh();
}

shinyjs.open = function(){
$('#tree').jstree(true).open_all();
}

shinyjs.deselect = function(){
$('#tree').jstree(true).deselect_all();
}

"

# UI ------------------------

shinyUI(fluidPage(
    useShinyjs(),
    includeCSS('www/style.css'),
    extendShinyjs(text = javaScript),
    titlePanel("Neuroexpresso"),
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
            fluidRow(
                column(4,htmlOutput(outputId = 'geneSearchHtml')),
                column(4, htmlOutput(outputId ='regionSelectHtml')),
                column(4,selectInput(inputId = 'platform',
                                     label = 'Select Platform',
                                     choices = c('GPL339','GPL1261')))
            ),
            
            
            textOutput(outputId = 'didYouMean'),
            textOutput(outputId = 'synonyms'),
            
            uiOutput('expressionUI'),
            fluidRow(
                column(2,
                       checkboxInput(inputId = 'color',
                                     label = 'Color?', 
                                     value = T)),
                column(3,
                       radioButtons(inputId = 'ordering',
                                    label = 'Order by',
                                    choices = c('Cell type','A-Z'), selected = NULL, inline = FALSE, width = NULL)),
                column(7, 
                       htmlOutput(outputId = 'selectTree'),
                       #htmlOutput('tree'))
                       shinyTree("tree",search = TRUE))
            )
        ),
        mainPanel(
            htmlOutput('warning'),
            ggvisOutput('expressionPlot')
        )
    )
))
