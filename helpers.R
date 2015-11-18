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

hierarchyNames = 'ShinyNames'

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
# second regions
regionGroups2 = regionize(mouseDes2,regionNames,prop)
names(regionGroups2) = sapply(names(regionGroups2), function(x){
    strsplit(x,split = '_')[[1]][1]
})


# deal with hierarchies 
# hierarchies = lapply(hierarchyNames, function(x){
#     des
# })

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
