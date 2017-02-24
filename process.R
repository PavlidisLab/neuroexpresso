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

mouseDes2 =
    read.design('/home/omancarci/wholeOtto/omancarci/brainGenesManuscript/data-raw/Mouse_Cell_Type_Data/n_expressoSamplesWithRNAseq2.tsv')
expression2 =
    read.exp('/home/omancarci/wholeOtto/omancarci/brainGenesManuscript/data-raw/Mouse_Cell_Type_Data/n_expressoExprWithRNAseq2.csv',check.names=FALSE)

mouseDes =
    read.design('/home/omancarci/wholeOtto/omancarci/brainGenesManuscript/data-raw/Mouse_Cell_Type_Data/n_expressoSamplesWithRNAseq.tsv')
expression =
    read.exp('/home/omancarci/wholeOtto/omancarci/brainGenesManuscript/data-raw/Mouse_Cell_Type_Data/n_expressoExprWithRNAseq.csv',check.names=FALSE)

mouseDes3 = 
    read.design('/home/omancarci/wholeOtto/omancarci/brainGenesManuscript/data-raw/Mouse_Cell_Type_Data/meltedSingleCells.tsv')
expression3 =
    read.exp('/home/omancarci/wholeOtto/omancarci/brainGenesManuscript/data-raw/Mouse_Cell_Type_Data/TasicPrimaryMeanComparable.csv',check.names=FALSE)

list[gene,exp] = sepExpr(expression)
list[gene2,exp2] = sepExpr(expression2)
list[gene3,exp3] = sepExpr(expression3)



exp = exp[!grepl('[|]',gene$Gene.Symbol),]
gene = gene[!grepl('[|]',gene$Gene.Symbol),]

exp2 = exp2[!grepl('[|]',gene2$Gene.Symbol),]
gene2 = gene2[!grepl('[|]',gene2$Gene.Symbol),]


gene %<>% select(Probe, Gene.Symbol,GemmaIDs,GeneNames,NCBIids)
gene2 %<>% select(Probe, Gene.Symbol, GemmaIDs,GeneNames,NCBIids)

exp = exp[,!is.na(mouseDes[prop])]
mouseDes = mouseDes[!is.na(mouseDes[prop]),]

exp2 = exp2[,!is.na(mouseDes2[prop])]
mouseDes2 = mouseDes2[!is.na(mouseDes2[prop]),]

exp3 =exp3[,!is.na(mouseDes3[prop])] 
mouseDes3 = mouseDes3[!is.na(mouseDes3[prop]),]


saveRDS(exp, file = 'Data/mouseExpr.rds')
saveRDS(mouseDes, file = 'Data/mouseDes.rds')
saveRDS(gene, file = 'Data/gene.rds')

saveRDS(exp2, file = 'Data/mouseExpr2.rds')
saveRDS(mouseDes2, file = 'Data/mouseDes2.rds')
saveRDS(gene2 ,file = 'Data/gene2.rds')


saveRDS(exp3, file = 'Data/mouseExpr3.rds')
saveRDS(mouseDes3, file = 'Data/mouseDes3.rds')
saveRDS(gene3 ,file = 'Data/gene3.rds')
#write.table(gene,'Data/mouseGene',sep=',',row.names=F)
#write.table(format(exp,digits=3),'Data/mouseExpr',sep=',',row.names=F,quote=F)

#write.design(mouseDes,'Data/meltedDesign.tsv')
#write.design(mouseDes2,'Data/meltedDesign2.tsv')


#write.table(gene2,'Data/mouseGene2',sep=',',row.names=F)
#write.table(format(exp2,digits=3),'Data/mouseExpr2',sep=',',row.names=F,quote=F)

# file.remove('Data/finalExp.csv')
# file.remove('Data/finalExp2.csv')

# re-creation of token
#token <- drop_auth()
#saveRDS(token, "droptoken.rds")

# to speed up regionization on startup
sourceGithub('oganm/brainGenesManuscript/R/regionize.R')
prop ='ShinyNames'
regionNames = 'Region'

#mouseDes = read.design('Data/meltedDesign.tsv')
mouseDes = mouseDes[match(colnames(exp),mouseDes$sampleName),]
regionGroups = memoReg(mouseDes,regionNames,prop,
                       regionHierarchy = regionHierarchy
)

#mouseDes2 = read.design('Data/meltedDesign2.tsv')
mouseDes2 = mouseDes2[match(colnames(exp2),mouseDes2$sampleName),]

regionGroups2 = memoReg(mouseDes2,regionNames,prop, 
                        regionHierarchy= regionHierarchy
)

mouseDes3 = mouseDes3[match(colnames(exp3),mouseDes3$sampleName),]
regionGroups3 = memoReg(mouseDes3,regionNames,prop,
                        regionHierarchy= regionHierarchy
)

save(memoReg,regionHierarchy,file = 'memoReg.rda')
purge()
