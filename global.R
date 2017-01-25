# initialization -----------
library(shiny)
library(ggvis)
library(RCurl)
library(GGally)
library(data.table)
library(httr)
library(ogbox)
library(geneSynonym)
library(dplyr)
library(magrittr)
library(shinyTree)
library(reshape2)
library(viridis)
library(lazyeval)
library(memoise)
library(data.tree)

print('now it begins')

#sourceGithub(OganM,toSource,'regionize.R')

#sourceGithub(OganM,brainGenesManuscript,'R/regionize.R')
load('memoReg.rda')

print('now we start')
#token <- readRDS("droptoken.rds")

if ((Sys.info()["nodename"])=='kent.pavlab.chibi.ubc.ca'){
    set_config(config(cainfo = '/home/omancarci/R/x86_64-unknown-linux-gnu-library/3.1/httr/cacert.pem')) 
}

print('one hand')


#drop_acc(dtoken = token)

outputDir = "Gene Searches"
prop ='ShinyNames'

regionNames = 'Region'

hierarchyNames = list(NeuronTypes = c('MajorType','Neurotransmitter','ShinyNames'),
                      Methodology = c('Method', 'Reference'))

# loading data ----------
# allDataPre = read.csv(paste0(outFolder,'/','finalExp.csv'), header = T)  
# allDataPre = allDataPre[!grepl('[|]',allDataPre$Gene.Symbol),]
# list[mouseGene, mouseExpr]=sepExpr(allDataPre)
# rm(allDataPre)

# mouseExpr = read.csv('Data/mouseExpr')
#mouseExpr =  fread('Data/mouseExpr',data.table = FALSE)
#mouseGene = read.csv('Data/mouseGene')
#mouseDes = read.design('Data/meltedDesign.tsv')

mouseExpr = readRDS('Data/mouseExpr.rds')
mouseGene = readRDS('Data/gene.rds')
mouseDes = readRDS('Data/mouseDes.rds')
mouseDes = mouseDes[match(colnames(mouseExpr),make.names(mouseDes$sampleName)),]
rownames(mouseExpr) = mouseGene$Gene.Symbol
print('data loaded 1')

#mouseExpr2 =  fread('Data/mouseExpr2',data.table = FALSE)
#mouseGene2 = read.csv('Data/mouseGene2')
#mouseDes2 = read.design('Data/meltedDesign2.tsv')

mouseExpr2 = readRDS('Data/mouseExpr2.rds')
mouseGene2 = readRDS('Data/gene2.rds')
mouseDes2 = readRDS('Data/mouseDes2.rds')
mouseDes2 = mouseDes2[match(colnames(mouseExpr2),make.names(mouseDes2$sampleName)),]


rownames(mouseExpr2) = mouseGene2$Gene.Symbol
# mouseExpr2 = mouseExpr2[,!mouseDes2$PyramidalDeep %in% 'Layer5Pyra',with=F]
# mouseDes2 = mouseDes2[!mouseDes2$PyramidalDeep %in% 'Layer5Pyra',]
print('data loaded 2')
print('one heart')

# load the region data -------
regionGroups = memoReg(mouseDes,regionNames,prop,
                         regionHierarchy = regionHierarchy
                         )
names(regionGroups) = sapply(names(regionGroups), function(x){
    strsplit(x,split = '_')[[1]][1]
})
# second regions
regionGroups2 = memoReg(mouseDes2,regionNames,prop, 
                          regionHierarchy= regionHierarchy
                          )
names(regionGroups2) = sapply(names(regionGroups2), function(x){
    strsplit(x,split = '_')[[1]][1]
})

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


# add region as an hierarchy currently working wierdly. check later
# hierarchDummy = regionHierarchy
# 
# members = hierarchDummy %>% unlist %>% names
# depth = members %>% strsplit(split = '\\.') %>% sapply(len) %>% max
# members = members %>% strsplit(split = '\\.')
# 
# for (i in 1:depth){
#     mouseDes[paste0('region',i)] = NA
#     layerMembers = members%>% sapply(function(x){x[i]})
#     for (j in unique(layerMembers)){
#         mouseDes[!is.na(regionGroups[[j]]),paste0('region',i)] = j
#     }
#     # clear na's for shorter branches
#     mouseDes[,paste0('region',i)][is.na(mouseDes[,paste0('region',i)])] =
#         mouseDes[is.na(mouseDes[,paste0('region',i)]),paste0('region',i-1)]
# }
# 
# for (i in 1:depth){
#     mouseDes2[paste0('region',i)] = NA
#     layerMembers = members%>% sapply(function(x){x[i]})
#     for (j in unique(layerMembers)){
#         mouseDes2[!is.na(regionGroups2[[j]]),paste0('region',i)] = j
#     }
#     # clear na's for shorter branches
#     mouseDes2[,paste0('region',i)][is.na(mouseDes2[,paste0('region',i)])] =
#         mouseDes2[is.na(mouseDes2[,paste0('region',i)]),paste0('region',i-1)]
# }
#
#hierarchyNames$BrainRegions = c(paste0('region',1:depth),'ShinyNames')

# deal with hierarchies 
hierarchies = lapply(hierarchyNames, function(levels){
    hierarchize(levels,mouseDes[!is.na(mouseDes[,levels[len(levels)]]),])
})





# some settings required for the plotting function -----
sourceGithub('oganm/brainGenesManuscript/R/cellColors.R')

coloring = cellColors()
coloring = c(coloring,
             ShreejoyGabaergic = 'pink',
             ShreejoyPurkinje = 'pink',
             "*Purkinje" = 'pink',
             ShreejoyPyramidal = 'pink',
             ShreejoyOligo = 'pink',
             ShreejoyAstrocyte= 'pink',
             ShreejoyThPosLC = 'pink')


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
        treeSelected = hierarchies[[treeChoice]] %>% unlist %>% names %>% gsub("^.*[.](?![ ,])",'',.,perl = T)
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
            out = selectFrom[,len(attr(x,'ancestry'))+1] %in% x[1] %>% which
            # limit selection to its ancestor in case there are leaves with the same name
            if (len(attr(x,'ancestry'))>0){
                limitToParent = lapply(1:len(attr(x,'ancestry')),function(i){
                    selectFrom[,tree[i]] %in% attr(x,'ancestry')[i] %>% which
                })
                out = intersectList(c(list(out),limitToParent))
            }
            return(out)
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




linked_brush2 <- function(keys, fill = "red") {
    stopifnot(is.character(fill), length(fill) == 1)
    
    rv <- shiny::reactiveValues(under_brush = character(), keys = character())
    rv$keys <- isolate(keys)
    
    input <- function(vis) {
        handle_brush(vis, fill = fill, on_move = function(items, ...) {
            rv$under_brush <- items$key__
        })
    }
    
    set_keys <- function(keys) {
        rv$keys <- keys
    }
    
    set_brush <- function(ids) {
        rv$under_brush <- ids
    }
    
    selected_r <- reactive(rv$keys %in% rv$under_brush)
    fill_r <- reactive(c("black", fill)[selected_r() + 1])
    
    list(
        input = input,
        selected = create_broker(selected_r),
        fill = create_broker(fill_r),
        set_keys = set_keys,
        set_brush = set_brush
    )
}