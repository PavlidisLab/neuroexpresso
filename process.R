library(dplyr)
library(magrittr)
library(ogbox)
library(memoise)
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

designs = list(GPL339 = read.design('/home/omancarci/wholeOtto/omancarci/brainGenesManuscript/data-raw/Mouse_Cell_Type_Data/n_expressoSamplesWithRNAseq.tsv'),
               GPL1261 = read.design('/home/omancarci/wholeOtto/omancarci/brainGenesManuscript/data-raw/Mouse_Cell_Type_Data/n_expressoSamplesWithRNAseq2.tsv'),
               RNAseq = read.design('/home/omancarci/wholeOtto/omancarci/brainGenesManuscript/data-raw/Mouse_Cell_Type_Data/meltedSingleCells.tsv'))

exprs  = list(GPL339 = read.exp('/home/omancarci/wholeOtto/omancarci/brainGenesManuscript/data-raw/Mouse_Cell_Type_Data/n_expressoExprWithRNAseq.csv',check.names=FALSE),
              GPL1261 =  read.exp('/home/omancarci/wholeOtto/omancarci/brainGenesManuscript/data-raw/Mouse_Cell_Type_Data/n_expressoExprWithRNAseq2.csv',check.names=FALSE),
              RNAseq =  read.exp('/home/omancarci/wholeOtto/omancarci/brainGenesManuscript/data-raw/Mouse_Cell_Type_Data/TasicPrimaryMeanComparable.csv',check.names=FALSE))


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


sourceGithub('oganm/brainGenesManuscript/R/regionize.R')
prop ='ShinyNames'
regionNames = 'Region'

for(i in 1:length(genes)){
    rownames(exprs[[i]]) = genes[[i]]$Gene.Symbol
    designs[[i]] = designs[[i]][match(make.names(colnames(exprs[[i]])),make.names(designs[[i]]$sampleName)),]
    
    regionGroups = memoReg(designs[[i]],regionNames,prop,
            regionHierarchy = regionHierarchy)
}

save(memoReg,regionHierarchy,file = 'memoReg.rda')


saveRDS(exprs, file = 'Data/exprs.rds')
saveRDS(genes, file = 'Data/genes.rds')
saveRDS(designs, file = 'Data/designs.rds')



# re-creation of token
#token <- drop_auth()
#saveRDS(token, "droptoken.rds")

# to speed up regionization on startup
