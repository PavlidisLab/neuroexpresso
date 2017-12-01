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

shinyjs.openStaticTree = function(){
    $('#staticRegionTree').jstree(true).open_all();
}

shyinjs.setStaticTree = function(params){
    eval(\"$('#staticRegionTree').jstree(true).settings.core.data=\"+params);
$('#staticRegionTree').jstree(true).refresh();
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
                                                                'setDefaultTree',
                                                                'openStaticTree',
                                                                'setStaticTree')),
                  titlePanel("NeuroExpresso"),
                  fluidRow(column(4,wellPanel(
                      tags$head(tags$script('$(function () { $("#expressionPlot").click(function(e){ $("#ggvis-tooltip").hide(); }); })')),
                      #inputIp("ipid"),
                      #inputUserid("fingerprint"),
                      tabsetPanel(id = 'tabs', 
                                  tabPanel('Gene Search', value = 'genes',
                                           fluidRow(
                                               column(4,
                                                      textInput(inputId = 'searchGenes',
                                                                     label = 'Select Gene',
                                                                     value = 'Dok5')
                                                                ),
                                               column(4,  selectInput(inputId = "regionChoice",
                                                                      label= 'Select Region',
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
                                  tabPanel('Help', value = 'help')),
                      conditionalPanel(condition = "input.tabs!='help'",
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
                         conditionalPanel(condition = "input.tabs!='help'",
                                          htmlOutput('warning'),
                                          ggvisOutput('difPlot'),
                                          ggvisOutput('expressionPlot'),
                                          wellPanel(id = 'difGenePanel',type='hidden',dataTableOutput('difGeneTable'))),
                         conditionalPanel(condition = "input.tabs == 'help'",
                                         helpPage()
                                          ),
                         br(),
                         wellPanel(
                             h3('NeuroExpresso data and marker genes'),
                             p('The data in NeuroExpresso compiled as part of a project aiming to select marker genes for brain cell types
                               and calculate marker gene profiles.'),
                             p('The data in neuroexpresso and marker genes identified in the study can be accessed ',
                               a(href = 'http://pavlab.msl.ubc.ca/supplement-to-mancarci-et-al-neuroexpresso/', 
                                 target="_blank",'here.')),
                             p('R package for calculation of marker gene profiles can be found ',
                               a(href='https://github.com/oganm/markerGeneProfile',
                                 target = '_blank','here.')),
                             h3('How to cite'),
                             p('If using NeuroExpresso or the data provided, please cite:'),
                             a(href = 'http://www.eneuro.org/content/early/2017/11/20/ENEURO.0212-17.2017',
                               target="_blank", 
                               'B. Ogan Mancarci et al., “Cross-Laboratory Analysis of Brain Cell Type Transcriptomes with Applications to Interpretation of Bulk Tissue Data,” ENeuro, November 20, 2017, ENEURO.0212-17.2017, https://doi.org/10.1523/ENEURO.0212-17.2017.'),
                             h3('Contact'),
                             p('If you have questions or problems, mail', 
                               a(href="mailto:ogan.mancarci@msl.ubc.ca",
                                 target= '_blank', 'Ogan Mancarci'),
                               ". Please mention NeuroExpresso by name to ensure avoiding spam detectors"),
                             p('To report bugs, open an issue on the',a(href="https://github.com/oganm/neuroexpresso/issues",target= '_blank', 'github repo'))
                         )
                  )
                  )))
