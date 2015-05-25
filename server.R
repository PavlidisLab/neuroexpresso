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

# load the region data -------
groupNames = 'PyramidalDeep'
regions =
    trimNAs(
        trimElement(
            unique(
                unlist(
                    strsplit(as.character(mouseDes[,regionNames]),',')))
            ,c('ALL','All','all','Cerebrum'))) #S pecial names
regionBased = expand.grid(groupNames, regions)
regionGroups = vector(mode = 'list', length = nrow(regionBased))
names(regionGroups) = paste0(regionBased$Var2,'_',regionBased$Var1)

for (i in 1:nrow(regionBased)){
    regionGroups[[i]] = mouseDes[,as.character(regionBased$Var1[i])]
    
    # remove everything except the region and ALL labeled ones. for anything but cerebellum, add Cerebrum labelled ones as well
    if (regionBased$Var2[i] == 'Cerebellum'){
        regionGroups[[i]][!grepl(paste0('(^|,)((',regionBased$Var2[i],')|((A|a)(L|l)(l|l)))($|,)'),mouseDes[,regionNames])] = NA
    } else {
        # look for cerebrums
        cerebrums = unique(regionGroups[[i]][grepl('(Cerebrum)',mouseDes[,regionNames])])
        
        # find which cerebrums are not represented in the region
        cerebString = paste(cerebrums[!cerebrums %in% regionGroups[[i]][grepl(paste0('(^|,)((',regionBased$Var2[i],')|((A|a)(L|l)(l|l)))($|,)'),mouseDes[,regionNames])]],
                            collapse = ')|(')
        
        # add them as well (or not remove them as well) with all the rest of the region samples
        regionGroups[[i]][(!grepl(paste0('(^|,)((',regionBased$Var2[i],')|((A|a)(L|l)(l|l)))($|,)'),mouseDes[,regionNames])
                           & !(grepl(paste0('(',cerebString,')'),mouseDes[,as.character(regionBased$Var1[i])]) & grepl('Cerebrum',mouseDes[,regionNames])))] =  NA   
    }
}

names(regionGroups) = regions

# some settings required for the plotting function -----
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

coloringPyramidal=c(coloringPyramidal,          
                    Pyramidal_Thy1 = 'turquoise',
                    PyramidalCorticoThalam = 'blue',
                    Pyramidal_Glt_25d2 = 'blue4',
                    Pyramidal_S100a10 ='deepskyblue3')

coloring = coloringPyramidal
prop ='PyramidalDeep'

plotSingle = function(gene, prop, coloring, field = 'Gene.Symbol'){
    mouseExpr = mouseExpr[,!is.na(mouseDes[,prop])]
    mouseDes = mouseDes[!is.na(mouseDes[,prop]),]
    
    mouseGene = mouseGene
    isNeuron = unique(cbind(mouseDes$MajorType,mouseDes[,prop]))
    frame = data.frame(t(mouseExpr[mouseGene[,field] %in% gene,]),mouseDes[,prop],mouseDes$Region,mouseDes$MajorType)
    
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

plotPretty = function(gene, prop, coloring, field = 'Gene.Symbol', regionSelect){
    mouseExpr = mouseExpr[,!is.na(regionSelect)]
    mouseDes = mouseDes[!is.na(regionSelect),]
    
    mouseGene = mouseGene
    isNeuron = unique(cbind(mouseDes$MajorType,mouseDes[,prop]))
    frame = data.frame(t(mouseExpr[mouseGene[,field] %in% gene,]),mouseDes[,prop],mouseDes$MajorType)
    
    names(frame) = c('gene','prop','Type')
    colors = toColor(mouseDes[,prop],coloring)
    manualColor = scale_colour_manual(name='prop', values = colors$palette)
    pal = colors$palette[order(names(colors$palette))]
    p =  ggplot(frame, aes(x = prop,y=gene))
    p = p + geom_point(aes(color = prop),size =4)
    
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
    p = p +  theme(legend.box = "horizontal")
    return(p)
}


print('starting stuff')
    
shinyServer(function(input, output) {
    
    
    output$expressionPlot = renderPlot({
        if (input$regionChoice =='.messy details'){
            plotSingle(input$geneSearch, prop, coloring, field = 'Gene.Symbol')
        } else if (input$regionChoice =='All'){
            plotPretty(input$geneSearch, prop, coloring, field = 'Gene.Symbol', mouseDes[,prop])
        }else {
            plotPretty(input$geneSearch, prop, coloring, field = 'Gene.Symbol', regionGroups[[input$regionChoice]])
        }
    }
    , height = 800, width = 800)
})