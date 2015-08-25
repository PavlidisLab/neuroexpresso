source('helpers.R')
#library(xkcd)

# functions ---------

plotPretty = function(gene,geneList,expression,design, prop, coloring, field = 'Gene.Symbol', regionSelect, jitter, pointSize, color,xSize ,ySize, yTitleSize, additionalGG){
    if (jitter){
        jitter = 'jitter'
    } else {
        jitter= 'identity'
    }
    mouseExpr = expression[,!is.na(regionSelect),with = F]
    mouseDes = design[!is.na(regionSelect),]
    
    mouseGene = geneList
    isNeuron = unique(cbind(mouseDes$MajorType,mouseDes[,prop]))
    frame = data.frame(t(mouseExpr[mouseGene[,field] %in% gene,]),mouseDes[,prop],mouseDes$MajorType)
    
    names(frame) = c('gene','prop','Type')
    frame$prop = as.character(frame$prop)
    frame$prop = factor(frame$prop,levels = isNeuron[order(isNeuron[,1],isNeuron[,2]),2])
    
    colors = toColor(mouseDes[,prop],coloring)
    manualColor = scale_fill_manual(name='prop', values = colors$palette)
    pal = colors$palette[order(names(colors$palette))]
    # yrange = c(min(frame$gene)-2,max(frame$gene)+2)
    # xrange = c(0,len(unique(frame$prop))+3)
    if (color){
        p = ggplot(frame, aes(x = prop,y=gene, fill = prop, group=prop))
        p = p + geom_point(color='black',pch=21,size = pointSize, position = jitter)
    } else {
        p =  ggplot(frame, aes(x = prop,y=gene))
        p = p + geom_point(fill='black',color ='white', pch=21,size = pointSize, position = jitter)
    }
    # p = p + geom_point(size=5)
    
    p = p + manualColor +
        theme(panel.background = element_rect(fill = "gray80"),legend.key = element_rect(fill = "gray80"))+
        theme(panel.grid.major=element_blank())+
        xlab('')+
        ylab(paste(gene,"log2 expression"))+
        theme(legend.title=element_blank(),
              legend.text = element_text(size=14)) + 
        theme(axis.ticks = element_blank(), axis.text.x = element_blank()) + 
        theme(axis.title.x = element_text(size=20),
              axis.text.x  = element_text(angle=90, vjust=0.5, size=xSize)) + 
        theme(axis.title.y = element_text(size=yTitleSize),
              axis.text.y = element_text(size=ySize))
    p = p +  theme(legend.box = "horizontal")  
    p = p + teval(additionalGG)
    return(p)
}


print('starting stuff')

# beginning of server -----------
searchedGenes = NULL
shinyServer(function(input, output, session) {
    session$onSessionEnded(function(){
        if (privacy){
            
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
    
    observe({
        fingerprint <<- input$fingerprint
        ipid <<- input$ipid
      #  privacy<<-input$privacyBox
      privacy <<- T
    })

#     reactive({
#         write.table( ,file = 'searchLog', append=T)
#     })
    output$expressionPlot = renderPlot({
        if (input$platform == 'GPL339'){
            selected = mouseGene$Gene.Symbol[tolower(mouseGene$Gene.Symbol) %in% tolower(input$geneSearch)]
        } else if (input$platform == 'GPL1261'){
            selected = mouseGene2$Gene.Symbol[tolower(mouseGene2$Gene.Symbol) %in% tolower(input$geneSearch)]
        }
        if (len(selected)==0){
            stop('Gene symbol not in the list')
        }
        
        print(as.character(selected))
        searchedGenes <<- c(searchedGenes, paste(as.character(selected),input$regionChoice))
        

        if (input$regionChoice =='All'){
            if (input$platform == 'GPL339'){
                plotPretty(selected,mouseExpr,mouseDes, prop, coloring, field = 'Gene.Symbol', mouseDes[,prop], input$jitterBox, input$pointSize, input$color,input$xSize, input$ySize, input$yTitleSize, '')
            } else if (input$platform == 'GPL1261'){
                plotPretty(selected,mouseExpr2,mouseDes2, prop, coloring, field = 'Gene.Symbol', mouseDes2[,prop], input$jitterBox, input$pointSize, input$color,input$xSize, input$ySize, input$yTitleSize, '')
            }
        }else {
            if (input$platform == 'GPL339'){
                plotPretty(selected,mouseGene,mouseExpr,mouseDes, prop, coloring, field = 'Gene.Symbol', regionGroups[[input$regionChoice]], input$jitterBox, input$pointSize, input$color,input$xSize, input$ySize, input$yTitleSize, '')
            } else if (input$platform == 'GPL1261'){
                plotPretty(selected,mouseGene2,mouseExpr2,mouseDes2, prop, coloring, field = 'Gene.Symbol', regionGroups2[[input$regionChoice]], input$jitterBox, input$pointSize, input$color,input$xSize, input$ySize, input$yTitleSize, '')   
            }
            
        }
    }
    , height = 700, width = 900)
    
    output$didYouMean = renderText({
        if(any(tolower(mouseGene$Gene.Symbol) %in% tolower(input$geneSearch))){
            return('')
        } else {
            return(paste('Did you mean:\n',paste(mouseGene$Gene.Symbol[order(adist(tolower(input$geneSearch), 
                                                           tolower(mouseGene$Gene.Symbol)))[1:5]],
                         collapse=', ')))
        }
    })
})
