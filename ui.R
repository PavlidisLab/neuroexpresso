library(shiny)
library(RCurl)
library(ggvis)
library(ogbox)
library(shinyTree)
library(shinyjs)

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

javaScript = "shinyjs.changeTree = function(params){
    eval(\"$('#tree').jstree(true).settings.core.data=\"+params);
    $('#tree').jstree(true).refresh();
    //$('#tree').jstree(true).open_all();
    //$('#tree').jstree(true).deselect_all();
}

shinyjs.changeTreeTry = function(){
 $('#tree').jstree(true).settings.core.data=[{'text' : 'root1', 'state' : {'selected' : true }, 'icon' : 'signal', 'children' : [{'text' : 'hede'},{'text' : 'hebe'}]},{'text' : 'root2', 'children' : [{'text' : 'SubListA', 'children' : [{'text' : 'leaf1'},{'text' : 'leaf2'},{'text' : 'leaf3'}]},{'text' : 'SubListB', 'state' : {'disabled' : true,'selected' : true }, 'children' : [{'text' : 'leafA'},{'text' : 'leafB'}]}]}];
 $('#tree').jstree(true).refresh();
 $('#tree').jstree(true).open_all();
 $('#tree').jstree(true).deselect_all();
}

shinyjs.changeTreeTrySimple = function(){
 $('#tree').jstree(true).settings.core.data=[{'text' : 'goygoy', 'state' : {'selected' : true }, 'icon' : 'signal', 'children' : [{'text' : 'hede'},{'text' : 'hebe'}]},{'text' : 'root2', 'children' : [{'text' : 'SubListA', 'children' : [{'text' : 'leaf1'},{'text' : 'leaf2'},{'text' : 'leaf3'}]},{'text' : 'SubListB', 'state' : {'disabled' : true,'selected' : true }, 'children' : [{'text' : 'leafA'},{'text' : 'leafB'}]}]}];
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
