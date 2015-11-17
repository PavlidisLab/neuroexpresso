source('helpers.R')
#library(xkcd)

# functions ---------

createFrame = function(gene,
                       geneList,
                       expression,
                       design,
                       prop,
                       reference = 'Reference',
                       pmid = 'PMID',
                       coloring,
                       field = 'Gene.Symbol',
                       regionSelect,
                       color = T){
    mouseExpr = expression[,!is.na(regionSelect),with = F]
    mouseDes = design[!is.na(regionSelect),]
    mouseGene = geneList
    isNeuron = unique(cbind(mouseDes$MajorType,mouseDes[,prop]))
    frame = data.frame(t(mouseExpr[mouseGene[,field] %in% gene,]),mouseDes[,prop],mouseDes$MajorType,mouseDes[,reference],mouseDes[,pmid])
    
    names(frame) = c('gene','prop','Type','reference','PMID')
    frame$prop = as.character(frame$prop)
    frame$prop = factor(frame$prop,levels = isNeuron[order(isNeuron[,1],isNeuron[,2]),2])
    
    colors = toColor(mouseDes[,prop],coloring)
    frame$color = apply(col2rgb(colors$cols),2,function(x){
        x = x/255
        rgb(x[1],x[2],x[3])
    })
    if (!color){
        frame$color = "#000000"
    }
    # this has to be unique because of gviss' key requirement
    frame$id = 1:nrow(frame)
    return(frame)
}


print('starting stuff')

# beginning of server -----------
shinyServer(function(input, output, session) {

    vals =  reactiveValues(fingerprint = '', ipid = '', searchedGenes = 'Ogn Cortex', querry = 'NULL')
    
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
            regQuerry  = names(regionGroups)[tolower(names(regionGroups)) %in% tolower(vals$querry$region)]
            
            selectInput(inputId = "regionChoice",
                        label= 'Select region',
                        selected = regQuerry,
                        # choices = c(regions,'All','.messy details')),
                        choices = c(names(regionGroups),'All'))
        } else{
            selectInput(inputId = "regionChoice",
                        label= 'Select region',
                        # choices = c(regions,'All','.messy details')),
                        choices = c(names(regionGroups),'All'))
        }
    })
    


    
    
    session$onSessionEnded(function(){
        privacy = TRUE
        
        if (privacy){
            isolate({
                fingerprint = vals$fingerprint
                ipid = vals$ipid
                searchedGenes = vals$searchedGenes
            })
            
            print(fingerprint)
            print(ipid)
        
            print(searchedGenes)
            
            files = drop_dir(outputDir)
            if ((ncol(files)>0)&&(paste0(ipid,'.',fingerprint) %in% unlist(apply(files[,1],2,basename)))){
                toWrite = drop_read_csv(unlist(files[unlist(apply(files[,1],2,basename)) %in% paste0(ipid,'.',fingerprint),1]), header=F)
            } else{
                toWrite = data.frame(V1= character(0))
            }
            toWrite = rbind(toWrite, data.frame(V1=searchedGenes))
            write.table(toWrite,paste0(ipid,'.',fingerprint) ,quote=F, row.names=F,col.names=F)
            drop_upload(file=paste0(ipid,'.',fingerprint), dest = outputDir)
        }
    }) 
    # for user fingerprinting
     observe({
         print('Fingerprinting done')
         vals$fingerprint = input$fingerprint
         vals$ipid = input$ipid
     })
    
    
    
    observe({
        vals$searchedGenes <- c(isolate(vals$searchedGenes), paste(as.character(gene()),region()))
        print(vals$searchedGenes)
    })
    
    
    
    # create frame as a reactive object to pass to ggvis
    frame = reactive({
#         if (input$platform == 'GPL339'){
#             selected = mouseGene$Gene.Symbol[tolower(mouseGene$Gene.Symbol) %in% tolower(input$geneSearch)]
#         } else if (input$platform == 'GPL1261'){
#             selected = mouseGene2$Gene.Symbol[tolower(mouseGene2$Gene.Symbol) %in% tolower(input$geneSearch)]
#         }
#         if (len(selected)==0){
#             return(oldFrame)
#             # stop('Gene symbol not in the list')
#         }
#         searchedGenes <<- c(searchedGenes, paste(as.character(selected),input$regionChoice))
        
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
        if (region() =='All'){
            frame = createFrame(selected,geneList,expression,design,prop,'Reference','PMID',coloring,'Gene.Symbol', design[,prop],input$color)
        
        } else {
            frame = createFrame(selected,geneList,expression,design,prop,'Reference',"PMID",coloring,'Gene.Symbol', regionSelect,input$color)
        }
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
        print('ello?')
        if (len(input$regionChoice)==0){
            return('Cortex')
        } else{
            input$regionChoice
        }
    })
    
    
    reactive({
    p = frame %>%
        ggvis(~prop,~gene,fill := ~color,key := ~id,size :=140 , stroke := 'black') %>%
        layer_points() %>%
        add_tooltip(function(x){
            # get links to GSM
            if (!grepl('GSM',rownames(frame())[x$id])){
                src = paste0('<p>Contact authors</p>',rownames(frame())[x$id])
            } else {
                src = paste0("<p><a target='_blank' href=",
                             "'http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",
                             rownames(frame())[x$id],
                             "'>",
                             rownames(frame())[x$id],
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
                                         labels = list(fontSize =14))) %>%
        add_axis('x',title='', properties = axis_props(labels = list(angle=45,
                                                                     align='left',
                                                                     fontSize = 20))) %>%
        set_options(height = 700, width = 750)
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
    output$warning = renderText({
        if(any(tolower(mouseGene$Gene.Symbol) %in% tolower(input$geneSearch))){
            return('')
        } else {
            stop('Gene symbol is not found!')
        }
    })
})
