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
library(reshape2)
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

# creates the tree to input for treejs given levels and a design file subset
hierarchize = function(levels,design){
    out = vector(mode = 'list', length = len(unique(design[levels[1]]) %>% trimNAs))
    
    out = lapply(out,function(x){structure('',stselected = TRUE)})
    names(out) = unique(design[levels[1]]) %>% trimNAs %>% sort
    
    if ((len(levels)>1) & (nrow(design)>0)){
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
             'FS Basket (G42)' = 'firebrick2',
             Stem = 'blue' ,
             Ependymal = 'orange',
             Serotonergic = 'darkolivegreen',
             Hypocretinergic = 'cadetblue',
             Dopaminergic = 'gray0',
             Th_positive_LC = 'blueviolet',
             GabaVIPReln = 'firebrick4',
             'VIPReln (G30)' = 'firebrick4',
             GabaRelnCalb = 'firebrick3',
             'Martinotti (GIN)' = 'firebrick3',
             GabaSSTReln = 'firebrick1',
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
                       order = 'Cell type',
                       treeChoice,
                       treeSelected){
    # browser()
    treeSelectedNames  = sapply(treeSelected,function(x){x[1]})
    names(treeSelected) = treeSelectedNames
    mouseExpr = expression[,!is.na(regionSelect),]
    mouseDes = design[!is.na(regionSelect),]
    mouseGene = geneList
    
    # if selection boxes are not yet loaded send an empty data frame
    if (len(treeChoice)==0){
        return(data.frame(GSM = character(0),
                          gene = double(0),
                          prop = character(0),
                          color = character(0),
                          reference = character(0),
                          PMID = character(0),
                          id = integer(0)))
    }
    
    if (len(treeSelected)==0){
        # treeSelected = design[hierarchyNames[[treeChoice]][len(hierarchyNames[[treeChoice]])]] %>% unique %>% trimNAs
        treeSelected = hierarchies[[treeChoice]] %>% unlist %>% names %>% gsub("^.*[.]",'',.)
    }
    
    # if (order == 'A-Z'){
    #     treeSelectedNames = sort(treeSelectedNames)
    # }
    
    tree = hierarchyNames[[treeChoice]]
    
    # to create groups to display have the fields relevant to the selected tree and find indexes of the choices in it
    selectFrom = mouseDes %>% select_(.dots=tree)
    # groups = lapply(treeSelectedNames, function(x){
    #     selectFrom %>% apply(1,function(y){x %in% y}) %>% which
    # })
    
    groups = lapply(treeSelected, function(x){
        if (is.null(attr(x,'ancestry'))){
            selectFrom[,len(tree)] %in% x %>% which
        } else{
            selectFrom[,len(attr(x,'ancestry'))+1] %in% x[1] %>% which
        }
        # selectFrom %>% apply(1,function(y){x %in% y}) %>% which
    })
    names(groups) = treeSelected # in case 
    
    
    while(groups %>% names %>% duplicated %>% any){
        names(groups)[groups %>% names %>% duplicated] %<>% paste0(' ')
    }
    
    if (order == 'A-Z'){
        groups = groups[groups %>% names %>% order]
    }
    expression = t(mouseExpr[mouseGene[,field] %in% gene,])
    
    frame = groups %>% melt
    colors = toColor(mouseDes$ShinyNames[frame$value],coloring)$col
    frame %<>% mutate(GSM = mouseDes$sampleName[value], 
                      gene = expression[value,], 
                      prop= L1,  
                      color = colors ,
                      reference = mouseDes$Reference[value], 
                      PMID = mouseDes$PMID[value]) %>% 
        select(GSM,gene,prop,color, reference,PMID)
    
    # amygdala fix. if a region doesnt exist returns an empty matrix
    if (nrow(frame)==0){
        frame[,2] = integer(0)
    }
    
    frame$color = apply(col2rgb(frame$color),2,function(x){
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

# turns a list acceptable by shinyTree into JSON format accetpable by jsTree.
toTreeJSON = function(list){
    if(length(list)==0){
        return("[{'text' : 'No cells in group'}]")
    }
    outString = '['
    for (i in 1:length(list)){
        outString %<>% paste0("{'text' : '",  names(list)[i], "'")
        attribs = attributes(list[[i]])
        
        stateAttribs = attribs[grepl('opened|disabled|selected',names(attribs))]
        children = attribs[grepl('names',names(attribs))]
        others = attribs[!grepl('opened|disabled|selected|names',names(attribs))]
        
        if (length(stateAttribs) >0){
            outString %<>% paste0(", 'state' : {")
            for (j in 1:length(stateAttribs)){
                outString %<>% paste0("'",gsub('st','',names(stateAttribs)[j]),"' : ", tolower(stateAttribs[j]))
                if (j < length(stateAttribs)){
                    outString %<>% paste0(",")
                }
            }
            outString %<>% paste('}')
        }
        
        if(length(others)>0){
            for (j in 1:length(others)){
                outString %<>% paste0( ", '",gsub('st','',names(others)[j]),"' : '", others[j],"'")
            }
        }
        
        if (class(list[[i]]) == 'list'){
            outString %<>% paste0(", 'children' : ",toTreeJSON(list[[i]]))
        }
        
        outString %<>% paste0("}")
        if (i < length(list)){
            outString %<>% paste0(",")
        }
        
    }
    outString %<>% paste0(']')
    return(outString)
}


