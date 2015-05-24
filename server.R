# initialization -----------
library(shiny)
library(ggplot2)
library(RCurl)
library(GGally)


eval( expr = parse( text = getURL(
    "https://raw.githubusercontent.com/oganm/toSource/master/ogbox.R",
    ssl.verifypeer=FALSE) ))
sourceGithub(oganm,masterOfCellTypes,runVars)

# loading data ----------
# allDataPre = read.csv(paste0(outFolder,'/','finalExp.csv'), header = T)  
# allDataPre = allDataPre[!grepl('[|]',allDataPre$Gene.Symbol),]
# list[mouseGene, mouseExpr]=sepExpr(allDataPre)
# rm(allDataPre)

mouseExpr = read.csv('Data/mouseExpr')
mouseGene = read.csv('Data/mouseGene')
mouseDes = read.design('Data/meltedDesign.tsv')
mouseDes = mouseDes[match(colnames(mouseExpr),make.names(mouseDes$sampleName),),]
rownames(mouseExpr) = mouseGene$Gene.Symbol

coloring = heatColors$CellType[4:len(heatColors$CellType)]
coloringDeep = c(heatColors$CellType[4:len(heatColors$CellType)],
                 GabaVIPReln = 'firebrick4',
                 GabaRelnCalb = 'firebrick3',
                 GabaSSTReln = 'firebrick1',
                 GabaReln = 'firebrick'
                 #GabaOxtr = 'orange',
                 #GabaHtr3a= 'indianred'
)

# some settings required for the plotting function -----
sourceGithub(oganm,masterOfCellTypes,runVars)
coloring = heatColors$CellType[4:len(heatColors$CellType)]

coloringDeep = c(heatColors$CellType[4:len(heatColors$CellType)],
                 GabaVIPReln = 'firebrick4',
                 GabaRelnCalb = 'firebrick3',
                 GabaSSTReln = 'firebrick1',
                 GabaReln = 'firebrick'
                 #GabaOxtr = 'orange',
                 #GabaHtr3a= 'indianred'
)

coloringDeep = coloringDeep[names(coloringDeep) != 'Gaba']


coloringPyramidal = coloringDeep[names(coloringDeep) != 'Pyramidal']

coloringPyramidal=c(coloringPyramidal,           Pyramidal_Thy1 = 'turquoise',
                    PyramidalCorticoThalam = 'blue',
                    Pyramidal_Glt_25d2 = 'blue4',
                    Pyramidal_S100a10 ='deepskyblue3')

coloring = coloringPyramidal
region="Region"
prop ='PyramidalDeep'

plotSingle = function(gene, prop, coloring, region, field = 'Gene.Symbol'){
    mouseExpr = mouseExpr[,!is.na(mouseDes[,prop])]
    mouseDes = mouseDes[!is.na(mouseDes[,prop]),]
    
    mouseGene = mouseGene
    isNeuron = unique(cbind(mouseDes$MajorType,mouseDes[,prop]))
    frame = data.frame(t(mouseExpr[mouseGene[,field] %in% gene,]),mouseDes[,prop],mouseDes[,region],mouseDes$MajorType)
    
    names(frame) = c('gene','prop','region','Type')
    colors = toColor(mouseDes[,prop],coloring)
    manualColor = scale_colour_manual(name='prop', values = colors$palette)
    pal = colors$palette[order(names(colors$palette))]
    p =  ggplot(frame, aes(x = prop,y=gene))
    p = p + geom_point(aes(color = prop,shape=region),size =4)
    
    p = p + manualColor +
        theme(panel.background = element_rect(fill = "gray60"),legend.key = element_rect(fill = "gray60"))+
        theme(panel.grid.major=element_blank())+
        xlab('')+
        ylab(paste(gene,"log2 expression"))+
        theme(legend.title=element_blank()) + 
        theme(axis.ticks = element_blank(), axis.text.x = element_blank()) + 
        theme(axis.title.x = element_text(size=20),
              axis.text.x  = element_text(angle=90, vjust=0.5, size=16)) + 
        theme(axis.title.y = element_text(size=20)) + 
        xlim(isNeuron[order(isNeuron[,1],isNeuron[,2]),2])
    p = p +  scale_shape_manual(values=c(0:10,35,12:18)) # someone will call me antisemite for this...
    p = p +  theme(legend.box = "horizontal")
    return(p)
}



    
shinyServer(function(input, output) {
    
    output$expressionPlot = renderPlot({
        plotSingle(input$geneSelect, prop, coloring, region, field = 'Gene.Symbol')
    }, height = 800, width = 800)
})