source('helpers.R')
library(shinyjs)
#library(xkcd)


print('starting stuff')

# beginning of server -----------
shinyServer(function(input, output, session) {
    
    vals =  reactiveValues(fingerprint = '', # will become user's fingerprint hash
                           ipid = '',  # will become user's ip address
                           searchedGenes = 'Ogn Cortex', # a history of searched genes to be saved
                           querry = 'NULL', # will become the querry string if something is being querried
                           hierarchies=NULL, # which hierarchy is selected
                           hierarchInit = NULL # will be set to the initial hierarchy
                           )
    
    observe({
        vals$querry = parseQueryString(session$clientData$url_search)
        print(vals$querry)
    })
    
    output$geneSearchHtml = renderUI({
        if (!is.null(vals$querry$gene)){
            textInput(inputId = 'geneSearch',value = vals$querry$gene,
                      label = 'Select Gene')
        } else {
            textInput(inputId = 'geneSearch',value = 'Ogn',
                      label = 'Select Gene')
        }
    })
    
    output$regionSelectHtml = renderUI({
        if (!is.null(vals$querry$region)){
            # browser()
            regQuerry  = names(regionGroups)[tolower(names(regionGroups)) %in% tolower(vals$querry$region)]

            selectInput(inputId = "regionChoice",
                        label= 'Select region',
                        selected = regQuerry,
                        # choices = c(regions,'All','.messy details')),
                        choices = names(regionGroups))
        } else{
            selectInput(inputId = "regionChoice",
                        label= 'Select region',
                        selected = 'Cortex',
                        # choices = c(regions,'All','.messy details')),
                        choices = names(regionGroups))
        }
    })
    
    
    # doesn't seem to work right. will have to remove it to add to the package anyway
    # session$onSessionEnded(function(){
    #     privacy = TRUE # this was there when logging was opt-outable. no more...
    #     
    #     if (privacy){
    #         isolate({
    #             fingerprint = vals$fingerprint
    #             ipid = vals$ipid
    #             searchedGenes = vals$searchedGenes
    #         })
    #         
    #         print(fingerprint)
    #         print(ipid)
    #         
    #         print(searchedGenes)
    #         
    #         files = drop_dir(outputDir, n = 0)
    #         if ((ncol(files)>0)&&(paste0(ipid,'.',fingerprint) %in% unlist(apply(files[,1],2,basename)))){
    #             toWrite = drop_read_csv(unlist(files[unlist(apply(files[,1],2,basename)) %in% paste0(ipid,'.',fingerprint),1]), header=F)
    #         } else{
    #             toWrite = data.frame(V1= character(0))
    #         }
    #         toWrite = rbind(toWrite, data.frame(V1=searchedGenes))
    #         write.table(toWrite,paste0(ipid,'.',fingerprint) ,quote=F, row.names=F,col.names=F)
    #         drop_upload(file=paste0(ipid,'.',fingerprint), dest = outputDir)
    #     }
    # }) 
    # # for user fingerprinting
    # observe({
    #     print('Fingerprinting done')
    #     vals$fingerprint = input$fingerprint
    #     vals$ipid = input$ipid
    # })
    
    
    
    observe({
        vals$searchedGenes <- c(isolate(vals$searchedGenes), paste(as.character(gene()),region()))
        print(vals$searchedGenes)
    })
    
    
    
    # create frame as a reactive object to pass to ggvis
    frame = reactive({
        selected = gene()
        
        if (input$platform == 'GPL339'){
            geneList = mouseGene
            expression = mouseExpr
            design = mouseDes
            regionSelect = regionGroups[[region()]]
        } else if (input$platform == 'GPL1261'){
            geneList = mouseGene2
            expression = mouseExpr2
            design = mouseDes2
            regionSelect = regionGroups2[[region()]]
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
                            input$color,
                            input$ordering,
                            input$treeChoice, get_selected(input$tree))
        
        # oldFrame <<- frame
        return(frame)
    })
    
    
    
    gene = reactive({
        if (input$platform == 'GPL339'){
            selected = mouseGene$Gene.Symbol[tolower(mouseGene$Gene.Symbol) %in% tolower(input$geneSearch)]
        } else if (input$platform == 'GPL1261'){
            selected = mouseGene2$Gene.Symbol[tolower(mouseGene2$Gene.Symbol) %in% tolower(input$geneSearch)]
        }
        if (len(selected)==0){
            isolate({
                return(strsplit(vals$searchedGenes[len(vals$searchedGenes)],' ')[[1]][1])
            })
        }
        return(selected)
    })
    
    region = reactive({
        if (len(input$regionChoice)==0){
            print('ello?')
            return('Cortex')
        } else{
            input$regionChoice
        }
    })
    
    # the plot has to be reactively created for its title to be prone to changes. anything that is not part of the frame
    # has to be passed inside this reactive block and will cause plot to be refreshed.
    reactive({
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
                     title = paste(gene(),"log2 expression"),
                     properties = axis_props(title = list(fontSize = 17),
                                             labels = list(fontSize =14)),
                     title_offset = 50) %>%
            add_axis('x',title='', properties = axis_props(labels = list(angle=45,
                                                                         align='left',
                                                                         fontSize = 20))) %>%
            set_options(height = input$plotHeight, width = input$plotWidth)
        return(p)
    }) %>%  bind_shiny('expressionPlot', 'expressionUI')
    
    output$didYouMean = renderText({
        if(any(tolower(mouseGene$Gene.Symbol) %in% tolower(input$geneSearch))){
            return('')
        } else {
            return(paste('Did you mean:\n',paste(mouseGene$Gene.Symbol[order(adist(tolower(input$geneSearch), 
                                                                                   tolower(mouseGene$Gene.Symbol)))[1:5]],
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
        if(any(tolower(mouseGene$Gene.Symbol) %in% tolower(input$geneSearch))){
            return('')
        } else {
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
        vals$hierarchies = switch(input$platform,
               GPL339 =  lapply(hierarchyNames, function(levels){
                   hierarchize(levels,mouseDes[!is.na(mouseDes[,levels[len(levels)]]) & !is.na(regionGroups[[region()]]),])
               }),
               GPL1261 = lapply(hierarchyNames, function(levels){
                   hierarchize(levels,mouseDes2[!is.na(mouseDes2[,levels[len(levels)]]) & !is.na(regionGroups2[[region()]]),])
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
