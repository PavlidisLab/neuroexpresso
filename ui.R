# UI ------------------------

shinyUI(fluidPage(theme = shinytheme('lumen'),
                  #tags$head(includeScript("www/js/google-analytics.js")),
                  useShinyjs(),
                  includeCSS('www/style.css'),
                  extendShinyjs(script = 'www/js/shinyjsFunctions.js',functions = c('changeTree',
                                                                'open',
                                                                'deselect',
                                                                'setDefaultTree',
                                                                'openStaticTree',
                                                                'setStaticTree',
                                                                'hidePlotTooltip')),
                  titlePanel("NeuroExpresso"),
                  fluidRow(column(4,wellPanel(
                      tags$head(tags$script('$(function () { $("#expressionPlot").click(function(e){ $("#ggvis-tooltip").hide(); }); })'),
                                tags$script('$(function () { $("#difPlot").click(function(e){ $("#ggvis-tooltip").hide(); }); })'),
                                tags$script('$(function () { $("#tabs").click(function(e){ $("#ggvis-tooltip").hide(); }); })')),
                      tabsetPanel(id = 'tabs', 
                                  tabPanel('Gene Search', value = 'genes'),
                                  tabPanel('Diff. Expr', value = 'difExp'),
                                  tabPanel('Help', value = 'help')
                      ),
                      
                      conditionalPanel(condition = "input.tabs=='difExp'",
                                       h4('Select groups on the graph to perform differential expression'),
                                       h4(id = "select_text", "Please select first group"),
                                       actionButton(inputId = "group1Selected", label = "Save group 1"),
                                       actionButton(inputId = "group2Selected", label = "Save group 2"),
                                       actionButton(inputId = 'newSelection', label = 'New Selection'),
                                       downloadButton(outputId = 'downloadDifGenes', label = 'Download')),
                      
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
                      conditionalPanel(condition = 'input.graphSettings.indexOf("Is marker?")>-1',htmlOutput(outputId = 'tisAMarker')),
                      htmlOutput(outputId ='allenBrain'),
                      br(),
                      #br(),
                      fluidRow(
                          column(6, 
                                 htmlOutput(outputId = 'selectTree'),
                                 #htmlOutput('tree'))
                                 shinyTree("tree",search = TRUE, checkbox = TRUE)),
                          column(6,
                                 fluidRow(column(5,
                                                 checkboxGroupInput(inputId = 'graphSettings',
                                                                    label = 'Options',
                                                                    selected = c('Color'),
                                                                    choices = c('Fixed Y axis',
                                                                                'Color',
                                                                                'Is marker?')),
                                                 checkboxGroupInput(inputId = 'display',
                                                                    label = 'Display',
                                                                    selected =  c('Microarray','RNAseq'),
                                                                    choices =  c('Microarray','RNAseq'))),
                                          column(7,
                                                 radioButtons(inputId = 'ordering',
                                                              label = 'Order by',
                                                              choices = c('Cell type','A-Z'), selected = NULL, inline = FALSE, width = NULL))
                                 )
                          )
                      )
                  ),
                  acknowledge()
                  ),
                  column(8,
                         conditionalPanel(condition = "input.tabs!='help'",
                                          htmlOutput('warning'),
                                          ggvisOutput('difPlot'),
                                          ggvisOutput('expressionPlot'),
                                          wellPanel(id = 'difGenePanel',type='hidden',dataTableOutput('difGeneTable'))),
                         conditionalPanel(condition = "input.tabs == 'help'",
                                          helpPage()
                         ),
                         br(),
                         bottomInfo()
                  )
                  )))
