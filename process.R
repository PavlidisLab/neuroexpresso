library(dplyr)
library(magrittr)
library(ogbox)
library(memoise)
library(allenBrain)
prop = 'ShinyNames'
# mouseDes2 =
#     read.design('/home/omancarci/wholeOtto/omancarci/brainGenesManuscript/data-raw/Mouse_Cell_Type_Data/n_expressoSamples2.tsv')
# expression2 =
#     read.exp('/home/omancarci/wholeOtto/omancarci/brainGenesManuscript/data-raw/Mouse_Cell_Type_Data/n_expressoExpr2.csv')
# 
# mouseDes =
#     read.design('/home/omancarci/wholeOtto/omancarci/brainGenesManuscript/data-raw/Mouse_Cell_Type_Data/n_expressoSamples.tsv')
# expression =
#     read.exp('/home/omancarci/wholeOtto/omancarci/brainGenesManuscript/data-raw/Mouse_Cell_Type_Data/n_expressoExpr.csv')

ogbox::getGithubFile(githubPath = 'PavlidisLab/neuroExpressoAnalysis/data-raw/Mouse_Cell_Type_Data/n_expressoSamplesWithRNAseq.tsv',downloadPath = 'Data/n_expressoSamplesWithRNAseq.tsv')
ogbox::getGithubFile(githubPath = 'PavlidisLab/neuroExpressoAnalysis/data-raw/Mouse_Cell_Type_Data/n_expressoSamplesWithRNAseq2.tsv',downloadPath = 'Data/n_expressoSamplesWithRNAseq2.tsv')
ogbox::getGithubFile(githubPath = 'PavlidisLab/neuroExpressoAnalysis/data-raw/Mouse_Cell_Type_Data/meltedSingleCells.tsv',downloadPath = 'Data/meltedSingleCells.tsv')

ogbox::getGithubFile(githubPath = 'PavlidisLab/neuroExpressoAnalysis/data-raw/Mouse_Cell_Type_Data/n_expressoExprWithRNAseq.csv',downloadPath = 'Data/n_expressoExprWithRNAseq.csv')
ogbox::getGithubFile(githubPath = 'PavlidisLab/neuroExpressoAnalysis/data-raw/Mouse_Cell_Type_Data/n_expressoExprWithRNAseq2.csv',downloadPath = 'Data/n_expressoExprWithRNAseq2.csv')
ogbox::getGithubFile(githubPath = 'PavlidisLab/neuroExpressoAnalysis/data-raw/Mouse_Cell_Type_Data/TasicPrimaryMeanComparable.csv',downloadPath = 'Data/TasicPrimaryMeanComparable.csv')


designs = list(GPL339 = read.design('Data/n_expressoSamplesWithRNAseq.tsv'),
               GPL1261 = read.design('Data/n_expressoSamplesWithRNAseq2.tsv'),
               RNAseq = read.design('Data/meltedSingleCells.tsv'))

exprs  = list(GPL339 = read.exp('Data/n_expressoExprWithRNAseq.csv',check.names=FALSE),
              GPL1261 =  read.exp('Data/n_expressoExprWithRNAseq2.csv',check.names=FALSE),
              RNAseq =  read.exp('Data/TasicPrimaryMeanComparable.csv',check.names=FALSE))

file.remove('Data/n_expressoSamplesWithRNAseq.tsv')
file.remove('Data/n_expressoSamplesWithRNAseq2.tsv')
file.remove('Data/meltedSingleCells.tsv')
file.remove('Data/n_expressoExprWithRNAseq.csv')
file.remove('Data/n_expressoExprWithRNAseq2.csv')
file.remove('Data/TasicPrimaryMeanComparable.csv')


exprs %<>% lapply(function(x){
    x %<>% filter(!grepl('[|]',Gene.Symbol))
    x %>% sepExpr
}) 
genes = exprs %>% lapply(function(x){
    out = x[[1]] 
    if('Probe' %in% names(out)){
       out  %<>% select(Probe, Gene.Symbol,GemmaIDs,GeneNames,NCBIids)
    }
    return(out)
})
exprs %<>% lapply(function(x){
    x[[2]]
})


sourceGithub('oganm/markerGeneProfile/R/regionize.R')
loadGithub('oganm/markerGeneProfile/data/mouseRegionHierarchy.rda')
loadGithub('oganm/markerGeneProfile/data/mouseMarkerGenesCombined.rda')
loadGithub('oganm/neuroExpressoAnalysis/data/publishableNameDictionary.rda')

mouseMarkerGenesCombined %<>% lapply(function(x){
    names(x) = replaceElement(names(x),dictionary = publishableNameDictionary$ShinyNames,labels = publishableNameDictionary$PyramidalDeep) %$% newVector
    return(x)
}) 

saveRDS(mouseMarkerGenesCombined,'Data/mouseMarkerGenesCombined')

regionHierarchy = mouseRegionHierarchy
prop ='ShinyNames'
regionNames = 'Region'

for(i in 1:length(genes)){
    rownames(exprs[[i]]) = genes[[i]]$Gene.Symbol
    designs[[i]] = designs[[i]][match(make.names(colnames(exprs[[i]])),make.names(designs[[i]]$sampleName)),]
    
    regionGroups = memoReg(designs[[i]],regionNames,prop,
            regionHierarchy = regionHierarchy)
}

save(memoReg,regionHierarchy,file = 'memoReg.rda')



# stuff to map tasic data to original source ------

ogbox::getGithubFile(githubPath = 'PavlidisLab/neuroExpressoAnalysis/data-raw/Mouse_Cell_Type_Data/singleCellMatchings.tsv',downloadPath = 'Data/singleCellMatchings.tsv')

rnaSeqMap = read.design('Data/singleCellMatchings.tsv')

saveRDS(rnaSeqMap, file = 'Data/rnaSeqMap.rds')



# saving of data ------------------
saveRDS(exprs, file = 'Data/exprs.rds')
saveRDS(genes, file = 'Data/genes.rds')
saveRDS(designs, file = 'Data/designs.rds')

minValue = exprs %>% sapply(function(x){x %>% unlist %>% trimNAs %>% min}) %>% min
maxValue = exprs %>% sapply(function(x){x %>% unlist %>% trimNAs %>% max}) %>% max

saveRDS(minValue, file = 'Data/minValue.rds')
saveRDS(maxValue, file = 'Data/maxValue.rds')

# get allen brain institute datasets for the genes ---------------
allGenes = genes %>% purrr::map('Gene.Symbol') %>% unlist %>% unique
allenIDs = seq_along(allGenes) %>% lapply(getGeneDatasets,function(x)planeOfSection = 'both')
pb = txtProgressBar(min = 0, max = length(allGenes), initial = 0,style = 3) 
allenIDs = seq_along(allGenes) %>% lapply(function(i){
    setTxtProgressBar(pb,i)
    getGeneDatasets(allGenes[[i]],planeOfSection = 'both')
})
allenIDs = allenIDs[sapply(allenIDs,length)>0]

saveRDS(allenIDs,file = 'Data/allenIDs.rds')

# re-creation of token
#token <- drop_auth()
#saveRDS(token, "droptoken.rds")

# to speed up regionization on startup
