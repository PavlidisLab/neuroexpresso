# initialization -----------
library(shiny)
library(ggvis)
library(ggplot2)
library(RCurl)
library(GGally)
library(data.table)
library(rdrop2)
library(httr)
library(ogbox)
library(geneSynonym)
library(dplyr)
library(magrittr)
library(shinyTree)

print('now it begins')


sourceGithub(OganM,brainCellTypeSpecificGenes,'R/regionize.R')

print('now we start')
token <- readRDS("droptoken.rds")

if ((Sys.info()["nodename"])=='kent.pavlab.chibi.ubc.ca'){
    set_config(config(cainfo = '/home/omancarci/R/x86_64-unknown-linux-gnu-library/3.1/httr/cacert.pem')) 
}

print('one hand')


drop_acc(dtoken = token)

outputDir = "Gene Searches"
prop ='ShinyNames'

regionNames = 'Region'

hierarchyNames = list(NeuronTypes = c('MajorType','Neurotransmitter1','ShinyNames'))

# loading data ----------
# allDataPre = read.csv(paste0(outFolder,'/','finalExp.csv'), header = T)  
# allDataPre = allDataPre[!grepl('[|]',allDataPre$Gene.Symbol),]
# list[mouseGene, mouseExpr]=sepExpr(allDataPre)
# rm(allDataPre)

# mouseExpr = read.csv('Data/mouseExpr')
mouseExpr =  fread('Data/mouseExpr',data.table = FALSE)
mouseGene = read.csv('Data/mouseGene')
mouseDes = read.design('Data/meltedDesign.tsv')
mouseDes = mouseDes[match(colnames(mouseExpr),make.names(mouseDes$sampleName)),]
rownames(mouseExpr) = mouseGene$Gene.Symbol
print('data loaded 1')
mouseExpr2 =  fread('Data/mouseExpr2',data.table = FALSE)
mouseGene2 = read.csv('Data/mouseGene2')
mouseDes2 = read.design('Data/meltedDesign2.tsv')
mouseDes2 = mouseDes2[match(colnames(mouseExpr2),make.names(mouseDes2$sampleName)),]
rownames(mouseExpr2) = mouseGene2$Gene.Symbol
# mouseExpr2 = mouseExpr2[,!mouseDes2$PyramidalDeep %in% 'Layer5Pyra',with=F]
# mouseDes2 = mouseDes2[!mouseDes2$PyramidalDeep %in% 'Layer5Pyra',]
print('data loaded 2')
print('one heart')

# load the region data -------
regionGroups = regionize(mouseDes,regionNames,prop)
names(regionGroups) = sapply(names(regionGroups), function(x){
    strsplit(x,split = '_')[[1]][1]
})
regionGroups$All =  mouseDes[,prop]
# second regions
regionGroups2 = regionize(mouseDes2,regionNames,prop)
names(regionGroups2) = sapply(names(regionGroups2), function(x){
    strsplit(x,split = '_')[[1]][1]
})
regionGroups2$All =  mouseDes2[,prop]


hierarchize = function(levels,design){
    out = vector(mode = 'list', length = len(unique(design[levels[1]]) %>% trimNAs))
    
    out = lapply(out,function(x){structure('',stselected = TRUE)})
    names(out) = unique(design[levels[1]]) %>% trimNAs %>% sort

    if (len(levels)>1){
        out = lapply(names(out),function(x){
            hierarchize(levels[-1] ,design[design[,levels[1]] %in% x,])
        })
        names(out) = unique(design[levels[1]]) %>% trimNAs %>% sort
        for(i in 1:len(out)){
            if (len(out[[i]])==1 && names(out[[i]]) == names(out[i])){
                out[[i]] = structure('',stselected = TRUE)}
        }
    }
    return(out)
}


# deal with hierarchies 
hierarchies = lapply(hierarchyNames, function(levels){
    hierarchize(levels,mouseDes[!is.na(mouseDes[,levels[len(levels)]]),])
})



# some settings required for the plotting function -----

coloring = c(Oligo = 'darkgreen',
             Oligodendrocyte = 'darkgreen',
             Bergmann = 'palegreen',
             MotorCholin = 'darkorange4',
             Cholin = 'darkorange',
             Cholinergic = 'darkorange',
             Spiny = 'blanchedalmond',
             Gluta = 'slategray',
             Basket = 'mediumpurple4',
             Golgi = 'orchid',
             Pyramidal = 'turquoise',
             Purkinje = 'purple',
             Inter = 'pink',
             CerebGranule = 'thistle',
             DentateGranule = 'thistle3',
             Microglia = 'white',
             # Gaba = 'firebrick4',
             Astrocyte = 'yellow',
             GabaPV = 'firebrick2',
             Stem = 'blue' ,
             Ependymal = 'orange',
             Serotonergic = 'darkolivegreen',
             Hypocretinergic = 'cadetblue',
             Dopaminergic = 'gray0',
             Th_positive_LC = 'blueviolet',
             GabaVIPReln = 'firebrick4',
             GabaRelnCalb = 'firebrick3',
             GabaSSTReln = 'firebrick1',
             GabaReln = 'firebrick',
             GabaVIPReln = 'firebrick4',
             GabaReln = 'firebrick',
             GabaOxtr = 'firebrick2',
             GabaHtr3a = 'darkred',
             Pyramidal_Thy1 = 'turquoise',
             PyramidalCorticoThalam = 'blue',
             Pyramidal_Glt_25d2 = 'blue4',
             Pyramidal_S100a10 ='deepskyblue3',
             Layer5Pyra = 'blue3'
)



print("Even death won't part us now.")


# frame output function --------

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
                       color = T,
                       order = 'Cell type'){
    # browser()
    mouseExpr = expression[,!is.na(regionSelect),]
    mouseDes = design[!is.na(regionSelect),]
    mouseGene = geneList
    frame = data.frame(mouseExpr %>% colnames ,t(mouseExpr[mouseGene[,field] %in% gene,]),mouseDes[,prop],mouseDes$MajorType,mouseDes$Neurotransmitter1,mouseDes[,reference],mouseDes[,pmid])
    # amygdala fix. if a region doesnt exist returns an empty matrix
    if (nrow(frame)==0){
        frame[,2] = integer(0)
    }
    
    names(frame) = c('GSM','gene','prop','Type', 'shinyNames2','reference','PMID')
    # browser()
    if (order=='Cell type'){
        frame %<>% arrange(Type, shinyNames2, prop)
    } else if (order =='A-Z'){
        frame %<>% arrange(prop)
    }
    frame$prop %<>% as.char %>% factor(levels = unique(frame$prop))
    
    
    colors = toColor(frame$prop,coloring)
    frame$color = apply(col2rgb(colors$cols),2,function(x){
        x = x/255
        rgb(x[1],x[2],x[3])
    })
    # if color is false, set all to black.
    if (!color){
        frame$color = "#000000"
    }
    # amygdala fix again
    if(nrow(frame)==0){
        frame = cbind(frame,data.frame(id=character(0)))
        return(frame)
    }
    # this has to be unique because of gviss' key requirement
    
    frame$id = 1:nrow(frame)
    return(frame)
}

