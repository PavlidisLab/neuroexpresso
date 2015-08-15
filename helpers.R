# initialization -----------
library(shiny)
library(ggplot2)
library(RCurl)
library(GGally)
library(data.table)
library(rdrop2)
library(httr)
token <- readRDS("droptoken.rds")
if ((Sys.info()["nodename"])=='kent.pavlab.chibi.ubc.ca'){
    set_config(config(cainfo = '/home/omancarci/R/x86_64-unknown-linux-gnu-library/3.1/httr/cacert.pem')) 
}

drop_acc(dtoken = token)

outputDir = "Gene Searches"

library(ogbox)
sourceGithub(oganm,masterOfCellTypes,runVars)

# loading data ----------
# allDataPre = read.csv(paste0(outFolder,'/','finalExp.csv'), header = T)  
# allDataPre = allDataPre[!grepl('[|]',allDataPre$Gene.Symbol),]
# list[mouseGene, mouseExpr]=sepExpr(allDataPre)
# rm(allDataPre)

# mouseExpr = read.csv('Data/mouseExpr')
mouseExpr =  fread('Data/mouseExpr')
mouseGene = read.csv('Data/mouseGene')
mouseDes = read.design('Data/meltedDesign.tsv')
mouseDes = mouseDes[match(colnames(mouseExpr),make.names(mouseDes$sampleName)),]
rownames(mouseExpr) = mouseGene$Gene.Symbol
print('data loaded 1')
mouseExpr2 =  fread('Data/mouseExpr2')
mouseGene2 = read.csv('Data/mouseGene2')
mouseDes2 = read.design('Data/meltedDesign2.tsv')
mouseDes2 = mouseDes2[match(colnames(mouseExpr2),make.names(mouseDes2$sampleName)),]
rownames(mouseExpr2) = mouseGene2$Gene.Symbol
mouseExpr2 = mouseExpr2[,!mouseDes2$PyramidalDeep %in% 'Layer5Pyra',with=F]
mouseDes2 = mouseDes2[!mouseDes2$PyramidalDeep %in% 'Layer5Pyra',]
print('data loaded 2')

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
##### second regions
groupNames = 'PyramidalDeep'
regions =
    trimNAs(
        trimElement(
            unique(
                unlist(
                    strsplit(as.character(mouseDes2[,regionNames]),',')))
            ,c('ALL','All','all','Cerebrum'))) #S pecial names
regionBased = expand.grid(groupNames, regions)
regionGroups2 = vector(mode = 'list', length = nrow(regionBased))
names(regionGroups2) = paste0(regionBased$Var2,'_',regionBased$Var1)

for (i in 1:nrow(regionBased)){
    regionGroups2[[i]] = mouseDes2[,as.character(regionBased$Var1[i])]
    
    # remove everything except the region and ALL labeled ones. for anything but cerebellum, add Cerebrum labelled ones as well
    if (regionBased$Var2[i] == 'Cerebellum'){
        regionGroups2[[i]][!grepl(paste0('(^|,)((',regionBased$Var2[i],')|((A|a)(L|l)(l|l)))($|,)'),mouseDes2[,regionNames])] = NA
    } else {
        # look for cerebrums
        cerebrums = unique(regionGroups2[[i]][grepl('(Cerebrum)',mouseDes2[,regionNames])])
        
        # find which cerebrums are not represented in the region
        cerebString = paste(cerebrums[!cerebrums %in% regionGroups2[[i]][grepl(paste0('(^|,)((',regionBased$Var2[i],')|((A|a)(L|l)(l|l)))($|,)'),mouseDes2[,regionNames])]],
                            collapse = ')|(')
        
        # add them as well (or not remove them as well) with all the rest of the region samples
        regionGroups2[[i]][(!grepl(paste0('(^|,)((',regionBased$Var2[i],')|((A|a)(L|l)(l|l)))($|,)'),mouseDes2[,regionNames])
                           & !(grepl(paste0('(',cerebString,')'),mouseDes2[,as.character(regionBased$Var1[i])]) & grepl('Cerebrum',mouseDes2[,regionNames])))] =  NA   
    }
}

names(regionGroups2) = regions



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
