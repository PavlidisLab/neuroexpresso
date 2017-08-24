#library(V8)
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

shinyjs.setDefaultTree = function(){
$.jstree.defaults.checkbox.three_state= false
}

"

# UI ------------------------

shinyUI(fluidPage(theme = shinytheme('lumen'),
                  #tags$head(includeScript("www/js/google-analytics.js")),
                  useShinyjs(),
                  includeCSS('www/style.css'),
                  extendShinyjs(text = javaScript,functions = c('changeTree',
                                                                'open',
                                                                'deselect',
                                                                'setDefaultTree')),
                  titlePanel("NeuroExpresso"),
                  fluidRow(column(4,wellPanel(
                      tags$head(tags$script('$(function () { $("#expressionPlot").click(function(e){ $("#ggvis-tooltip").hide(); }); })')),
                      #inputIp("ipid"),
                      #inputUserid("fingerprint"),
                      h1('About'),
                      p('This application aims to make it easier to visualize gene expression in mouse brain cell types. The data here is compiled for a project aiming to select cell type specific genes in brain.'),
                      p('It is compiled by combining data from', 
                        a(href="http://www.chibi.ubc.ca/Gemma/arrays/showArrayDesign.html?id=7",target= '_blank', 'GPL339'),',',
                        a(href="http://www.chibi.ubc.ca/Gemma/arrays/showArrayDesign.html?id=3",target= '_blank' ,'GPL1261'), 'and RNA-seq data by',
                        a(href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71585",target= '_blank' ,'Tasic et al.')),
                      p('RNA-seq data is normalized to match distribution of microarray data for plotting purposes.'),
                      p('To see genes that are only available for GPL1261 or RNA-seq data, choose that platform below. This will remove some of the samples'),
                      p('Click on data points to see their sources'),
                      a(href="https://github.com/oganm/neuroexpresso",target= '_blank', 'Source Code'),
                      br(),
                      a(href="https://github.com/oganm/brainGenesManuscript",target= '_blank', 'Project Page'),
                      br(),
                      p('Wait until the plot renders then enter a gene symbol.'),
                      tabsetPanel(id = 'tabs', 
                                  tabPanel('Gene Search', value = 'genes',
                                           fluidRow(
                                               column(4,
                                                      textInput(inputId = 'searchGenes',
                                                                     label = 'Select Gene',
                                                                     value = 'Dok5')
                                                                ),
                                               column(4,  selectInput(inputId = "regionChoice",
                                                                      label= 'Select region',
                                                                      selected = 'Cortex',
                                                                      choices = names(regionGroups$GPL339))),
                                               column(4,selectInput(inputId = 'platform',
                                                                    label = 'Select Platform',
                                                                    choices =  names(exprs)))
                                           ),
                                           
                                           
                                           htmlOutput(outputId = 'didYouMean'),
                                           textOutput(outputId = 'synonyms'),
                                           
                                           uiOutput('expressionUI')
                                  ),
                                  tabPanel('Diff. Expr', value = 'difExp',
                                           h4('Select groups on the graph to perform differential expression'),
                                           h4(id = "select_text", "Please select first group"),
                                           actionButton(inputId = "group1Selected", label = "Save group 1"),
                                           actionButton(inputId = "group2Selected", label = "Save group 2"),
                                           actionButton(inputId = 'newSelection', label = 'New Selection'),
                                           downloadButton(outputId = 'downloadDifGenes', label = 'Download')),
                                  tabPanel('Marker Genes', value = 'marker')),
                      conditionalPanel(condition = "input.tabs!='marker'",
                                       #br(),
                                       fluidRow(
                                           column(7, 
                                                  htmlOutput(outputId = 'selectTree'),
                                                  #htmlOutput('tree'))
                                                  shinyTree("tree",search = TRUE, checkbox = TRUE)),
                                           column(5,
                                                  fluidRow(column(5,
                                                                  checkboxGroupInput(inputId = 'graphSettings',
                                                                                     label = 'Options',
                                                                                     selected = c('Color'),
                                                                                     choices = c('Fixed Y axis',
                                                                                                 'Color')),
                                                                  checkboxGroupInput(inputId = 'display',
                                                                                     label = 'Display',
                                                                                     selected =  c('Microarray','RNAseq'),
                                                                                     choices =  c('Microarray','RNAseq'))),
                                                           column(7,
                                                                  radioButtons(inputId = 'ordering',
                                                                               label = 'Order by',
                                                                               choices = c('Cell type','A-Z'), selected = NULL, inline = FALSE, width = NULL))
                                                  )# ,
                                                  # sliderInput('','width',
                                                  #             min = 20, max  = 1500, sep = '', post = ' px',
                                                  #             value = 750),
                                                  # sliderInput('plotHeight','height',
                                                  #             min = 20, max  = 1500, sep = '', post = ' px',
                                                  #             value = 700)
                                                  )
                                       ))
                      # tabPanel('Bulk Search',value = 'bulk', p('yo')))
                  ),
                  wellPanel(
                      p('Developed and maintained by',
                        a(href="https://github.com/oganm/", 'Ogan Mancarci'), 'at',
                        a(href="http://www.chibi.ubc.ca/faculty/paul-pavlidis/pavlidis-lab/", 'Pavlidis Lab')
                      ),
                      br(),
                      p('Supported by'),
                      a(href = 'http://neurodevnet.ca/',
                        img(src = 'NeuroDevNet.png',height = '100')),
                      a(href= 'http://www.nih.gov/',
                        img(src='NIH.png', height = '100')),
                      a(href = 'http://www.cihr-irsc.gc.ca/e/193.html',
                        img(src = 'cihr.png', height = '100'))
                  )),
                  column(7,
                         conditionalPanel(condition = "input.tabs!='marker'",
                                          htmlOutput('warning'),
                                          ggvisOutput('difPlot'),
                                          ggvisOutput('expressionPlot'),
                                          wellPanel(id = 'difGenePanel',type='hidden',dataTableOutput('difGeneTable'))),
                         conditionalPanel(condition = "input.tabs == 'marker'",
                                          wellPanel(h3('Marker Genes'),
                                                    p('The list of marker genes identified in the study can be accessed ',
                                                      a(href = 'http://www.chibi.ubc.ca/supplement-to-mancarci-et-al-neuroexpresso/', 
                                                        target="_blank",'here.')))),
                         br(),
                         wellPanel(h3('How to cite'),
                                   p('If using NeuroExpresso or the data provided, please cite our pre-print:'),
                                   p('Mancarci BO, Toker L, Tripathy SJ, Li B, Rocco B, Sibille E, et al. Cross-laboratory analysis of brain cell type transcriptomes with applications to interpretation of bulk tissue data. bioRxiv. 2017 May 18;89219.'),
                                   p('Or wait for us to publish the paper. (This box will be updated)')
                                   )
                  )
)))
