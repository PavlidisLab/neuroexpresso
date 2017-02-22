
print('starting stuff')

# beginning of server -----------
shinyServer(function(input, output, session) {
    
    lb = linked_brush2(keys = NULL, "red")
    
    vals =  reactiveValues(fingerprint = '', # will become user's fingerprint hash
                           ipid = '',  # will become user's ip address
                           searchedGenes = 'Ogn Cortex GPL339', # a history of searched genes to be saved
                           gene = 'Ogn',
                           region = 'Cortex',
                           platform = 'GPL339',
                           treeChoice = NULL,
                           treeSelected = list(),
                           new = TRUE, # did the app just start?
                           querry = 'NULL', # will become the querry string if something is being querried
                           hierarchies=NULL, # which hierarchy is selected
                           hierarchInit = NULL, # will be set to the initial hierarchy
                           difGroup1= NULL,
                           difGroup2 = NULL,
                           differentiallyExpressed = NULL)
    hide(id="group2Selected")
    hide(id = 'newSelection')
    hide(id = 'downloadDifGenes')
    hide(id = 'difGenePanel')
    
    observe({
        if(input$tabs=='genes'){
            hide(id = 'difPlot-container')
            show(id = 'expressionPlot-container')
            enable(selector = "input[value = 'Fixed Y axis']")
        } else if(input$tabs == 'difExp'){
            hide(id = 'expressionPlot-container')
            show(id = 'difPlot-container')
            disable(selector = "input[value = 'Fixed Y axis']")
        }
    })
    
    # querry parsing ------
    observe({
        vals$querry = parseQueryString(session$clientData$url_search)
        print(vals$querry)
    })
    
    # observe({
    #     if (!is.null(vals$querry$gene)){
    #         updateTextInput(session,
    #                         inputId = 'geneSearch',value = vals$querry$gene,
    #                         label = 'Select Gene')
    #     }
    # })
    
    # this shouldn't be necesarry but first plot does not show the axis labels for no real reason.
    observe({
        if(vals$new & is.null(vals$querry$gene)){
            # this replaces the initial configuration of the textInput so that the plot will be drawn twice
            # this shouldn't be necesarry, alas it seems to be...
            updateTextInput(session,
                            inputId = 'geneSearch',value = 'Ogn',
                            label = 'Select Gene')
            vals$new= TRUE
        } else if (vals$new & !is.null(vals$querry$gene)){
            updateTextInput(session,
                            inputId = 'geneSearch',value = vals$querry$gene,
                            label = 'Select Gene')
            vals$new= TRUE
        }

    })
    
    
    observe({
        if (!is.null(vals$querry$region)){
            # browser()
            regQuerry  = names(regionGroups)[tolower(names(regionGroups)) %in% tolower(vals$querry$region)]
            
            updateSelectInput(session,
                              inputId = "regionChoice",
                              label= 'Select region',
                              selected = regQuerry,
                              choices = names(regionGroups))
            # #planned feature
            # selected = nametreeVector(regionHierarchy)[
            #     str_extract(regionHierarchy %>% nametreeVector,'([A-Z]|[a-z])*?(?=\n)') %in% regQuerry
            #     ],
            # choices = regionHierarchy %>% nametreeVector)
        }
    })
    
    
    # check validity of input -----------
    observe({
        if (input$platform == 'GPL339'){
            selected = mouseGene$Gene.Symbol[tolower(mouseGene$Gene.Symbol) %in% tolower(input$geneSearch)]
        } else if (input$platform == 'GPL1261'){
            selected = mouseGene2$Gene.Symbol[tolower(mouseGene2$Gene.Symbol) %in% tolower(input$geneSearch)]
        }
        if (len(selected)>0){
            vals$searchedGenes <- c(isolate(vals$searchedGenes), paste(selected,input$regionChoice,input$platform))
            vals$gene = strsplit(vals$searchedGenes[len(vals$searchedGenes)],' ')[[1]][1]
            vals$region = strsplit(vals$searchedGenes[len(vals$searchedGenes)],' ')[[1]][2]
            vals$platform = strsplit(vals$searchedGenes[len(vals$searchedGenes)],' ')[[1]][3]
            print(vals$searchedGenes)
        }
    })
    
    # validity of trees
    observe({
        treeSelected = get_selected(input$tree)
        treeChoice =  input$treeChoice

        validTree = hierarchies %>% sapply(function(x){
            hiearNames = x %>% unlist %>% names %>% strsplit('\\.(?![ ,])',perl=TRUE) %>% unlist %>% unique
            all(treeSelected %in% hiearNames) | length(treeSelected)==0
        })
        
        if(!is.null(treeChoice)){
            if(validTree[treeChoice]){
                vals$treeSelected = treeSelected
                vals$treeChoice =  treeChoice
            }
        }
        # vals$treeSelected = treeSelected
        # vals$treeChoice =  treeChoice
    })
    
    # create frame as a reactive object to pass to ggvis
    frame = reactive({
        selected = vals$gene
        region =  vals$region
        platform = vals$platform
        
        if (platform == 'GPL339'){
            geneList = mouseGene
            expression = mouseExpr
            design = mouseDes
            regionSelect = regionGroups[[region]]
        } else if (platform == 'GPL1261'){
            geneList = mouseGene2
            expression = mouseExpr2
            design = mouseDes2
            regionSelect = regionGroups2[[region]]
        }
        frame = createFrame(selected,
                            geneList,
                            expression,
                            design,
                            prop,
                            'Reference',
                            "PMID",
                            coloring,
                            'Gene.Symbol', 
                            regionSelect,
                            'Color' %in% input$graphSettings,
                            input$ordering,
                            vals$treeChoice, vals$treeSelected)
        
        # oldFrame <<- frame
        lb$set_keys(1:nrow(frame))
        return(frame)
    })
    
    # differential expression -----------------
    observe({
        print(input$group1Selected)
        isolate({
            if(input$group1Selected>=1){
                if(!any(lb$selected())){
                    selected = TRUE
                } else {
                    selected = lb$selected()
                }
                print('save 1')
                hide(id = 'group1Selected')
                show(id = 'group2Selected')
                vals$difGroup1 = frame()[selected,]$GSM
                print(vals$difGroup1)
            }
        })
    })
    
    observe({
        print(input$group2Selected)
        isolate({
            if(input$group2Selected>=1){
                if(!any(lb$selected())){
                    selected = TRUE
                } else {
                    selected = lb$selected()
                }
                print('save 2')
                hide(id = 'group2Selected')
                show(id = 'newSelection')
                show(id = 'downloadDifGenes')
                vals$difGroup2 = frame()[selected,]$GSM
                print(vals$difGroup2)
                if(all(c(vals$difGroup1,vals$difGroup2) %in% colnames(mouseExpr2))){
                    toDif = mouseExpr2[c(vals$difGroup1,vals$difGroup2)]
                } else if (all(c(vals$difGroup1,vals$difGroup2) %in% colnames(mouseExpr))){
                    toDif = mouseExpr[c(vals$difGroup1,vals$difGroup2)]
                } else{
                    stop('WTF')
                }
                mm = model.matrix(~ groups,
                                  data.frame(groups = c(rep('a',length(vals$difGroup1)),
                                                        rep('b',length(vals$difGroup2)))))
                
                fit <- limma::lmFit(toDif, mm)
                fit <- limma::eBayes(fit)
                dif = limma::topTable(fit, coef=colnames(fit$design)[2],
                                      number = Inf)
                vals$differentiallyExpressed = data.frame(Symbol = rownames(dif),dif)
                
                show('difGenePanel')
            }
        })
    })
    output$difGeneTable = renderDataTable({
        vals$differentiallyExpressed %<>% 
            mutate(logFC =  round(logFC, digits=3),
                   AveExpr =  round(AveExpr, digits=3),
                   t =  round(t, digits=3),
                   P.Value =  round(P.Value, digits=3),
                   adj.P.Val =  round(adj.P.Val, digits=3),
                   B =  round(B, digits=3))
        datatable(vals$differentiallyExpressed,selection = 'single')
    })
    
    
    observe({
        input$difGeneTable_rows_selected
        isolate({
            if(!is.null(input$difGeneTable_rows_selected)){
                tableGene = vals$differentiallyExpressed[input$difGeneTable_rows_selected,'Symbol']
                updateTextInput(session,'geneSearch','Select Gene', value = tableGene)
                }
        })
    })
    
    output$downloadDifGenes = downloadHandler(
        filename = 'difGenes.tsv',
        content = function(file) {
            write_tsv(vals$differentiallyExpressed, file)
        })
    
    observe({
        print(input$newSelection)
        isolate({
            if(input$newSelection>=1){
                print('new selection')
                hide(id = 'newSelection')
                hide(id = 'downloadDifGenes')
                hide(id = 'difGenePanel')
                show(id = 'group1Selected')
                print(lb$selected())
            }
        })
    })
    
    
    frame %>%  #hede %>%
        ggvis(~prop,~gene,fill := ~color,
              key := ~id,size :=140 ,
              stroke := 'black',
              opacity := 0.7,
              size.brush := 300) %>%
        layer_points() %>% lb$input() %>%
        set_options(height = 700, width = 750) %>%
        add_axis('x',
                 #shame shame shame! but its mostly ggvis' fault
                 title='                                                                                                                                                                                                                            ',
                 properties = axis_props(labels = list(angle=45,
                                                       align='left',
                                                       fontSize = 20)),
                 title_offset = 150) %>%
        add_axis('y',
                 title = 'log2 expression',
                 properties = axis_props(title = list(fontSize = 17),
                                         labels = list(fontSize =14)),
                 title_offset = 50) %>%
        bind_shiny('difPlot')
    
    
    # the plot has to be reactively created for its title to be prone to changes. anything that is not part of the frame
    # has to be passed inside this reactive block and will cause plot to be refreshed.
    reactive({
        gene = vals$gene
        p = frame %>%
            ggvis(~prop,~gene,fill := ~color,key := ~id,size :=140 , stroke := 'black', opacity := 0.7) %>%
            layer_points() %>%
            add_tooltip(function(x){
                # get links to GSM
                if (!grepl('GSM',frame()$GSM[x$id])){
                    src = paste0('<p>Contact authors</p>',frame()$GSM[x$id])
                } else {
                    src = paste0("<p><a target='_blank' href=",
                                 "'http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",
                                 frame()$GSM[x$id],
                                 "'>",
                                 frame()$GSM[x$id],
                                 '</a></p>')
                }
                print(src)
                
                # get link to pubmed
                if(is.na(frame()$PMID[x$id])){
                    paper = paste0('<p>',as.character(frame()[x$id,]$reference),'</p>')
                } else {
                    paper = paste0("<p><a target='_blank' href=",
                                   "'http://www.ncbi.nlm.nih.gov/pubmed/",
                                   frame()$PMID[x$id],"'>",
                                   as.character(frame()[x$id,]$reference),
                                   '</a></p>')
                }
                
                print(paper)
                
                return(paste0(paper,src))
            }, on = 'click') %>% 
            add_axis('y',
                     title = paste(gene,"log2 expression"),
                     properties = axis_props(title = list(fontSize = 17),
                                             labels = list(fontSize =14)),
                     title_offset = 50) %>%
            add_axis('x',title='', properties = axis_props(labels = list(angle=45,
                                                                         align='left',
                                                                         fontSize = 20))) %>%
            set_options(height = 700, width = 750)
        if('Fixed Y axis' %in% input$graphSettings){
            p = p %>% scale_numeric('y',domain = c(minValue,maxValue))
        }
        
        return(p)
    }) %>%  bind_shiny('expressionPlot', 'expressionUI')
    
    output$didYouMean = renderText({
        if(any(tolower(mouseGene$Gene.Symbol) %in% tolower(input$geneSearch)) & input$platform =='GPL339'){
            return('')
        } else if (any(tolower(mouseGene2$Gene.Symbol) %in% tolower(input$geneSearch)) & input$platform =='GPL1261'){
            return('')
        }else {
            if (input$platform =='GPL339'){
                symbolList = mouseGene$Gene.Symbol
            } else if (input$platform =='GPL1261'){
                symbolList = mouseGene2$Gene.Symbol
            }
            return(paste('Did you mean:\n',paste(symbolList[order(adist(tolower(input$geneSearch), 
                                                                                   tolower(symbolList)))[1:5]],
                                                 collapse=', ')))
        }
    })
    
    # find if entered gene is a synonym of something else
    output$synonyms = renderText({
        synos = mouseSyno(input$geneSearch)[[1]]
        synos = sapply(synos, function(x){
            if(!x[1] == input$geneSearch){
                return(x[1])
            } else{
                return(NULL)
            }
        }
        ) %>% unlist
        
        if(len(synos) > 0 ){
            return(paste('Synonym of:', paste(synos,collapse=',')))
        } else{
            return('')
        }
    })
    
    
    output$warning = renderText({
        if(any(tolower(mouseGene$Gene.Symbol) %in% tolower(input$geneSearch)) & input$platform =='GPL339'){
            return('')
        } else if (any(tolower(mouseGene2$Gene.Symbol) %in% tolower(input$geneSearch)) & input$platform =='GPL1261'){
            return('')
        }else {
            stop('Gene symbol is not found!')
        }
    })
    
    
    
    observe({
        # print(get_selected(input$tree))
        # browser()
        region = 'Cortex'
        if (!is.null(vals$querry$region)){
            region = names(regionGroups)[tolower(names(regionGroups)) %in% tolower(vals$querry$region)]
        }
        levels = hierarchyNames[[1]]
        vals$hierarchInit = hierarchize(levels, mouseDes[!is.na(mouseDes[,levels[len(levels)]]) & !is.na(regionGroups[[region]]),])
    })
    
    # generate new hierarchies every time region changes
    observe({
        # browser()
        vals$hierarchies = switch(vals$platform,
               GPL339 =  lapply(hierarchyNames, function(levels){
                   hierarchize(levels,mouseDes[!is.na(mouseDes[,levels[len(levels)]]) & !is.na(regionGroups[[vals$region]]),])
               }),
               GPL1261 = lapply(hierarchyNames, function(levels){
                   hierarchize(levels,mouseDes2[!is.na(mouseDes2[,levels[len(levels)]]) & !is.na(regionGroups2[[vals$region]]),])
               }))
        # browser()
        
        # vals$hierarchies = lapply(hierarchyNames, function(levels){
        #     hierarchize(levels,mouseDes[!is.na(mouseDes[,levels[len(levels)]]) & !is.na(regionGroups[[region()]]),])
        # })
    })
    
    #     observe({
    #         frame()
    #         print('coinize')
    #         delay(1000,
    #               js$coinPlot())
    #     })
    
    observe({
        if (!is.null(input$treeChoice)){
            # browser()
            jsInput = toTreeJSON(vals$hierarchies[[input$treeChoice]])
            js$changeTree(jsInput) 
            delay(500,{
                js$open()
                js$deselect()
            })
        }
    })
    
    
    # choice of tree
    output$selectTree = renderUI({
        #browser()
        selectInput(inputId = "treeChoice",
                    label= 'Select hierarchy',
                    selected = names(hierarchyNames)[1],
                    choices = names(hierarchyNames))
    })
    
    
    output$tree = renderTree({
        # browser()
        js$setDefaultTree()
        #vals$hierarchies[[input$treeChoice]]
        vals$hierarchInit
    })
    
})
