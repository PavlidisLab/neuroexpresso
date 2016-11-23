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


list[gene,exp] = sepExpr(expression)
list[gene2,exp2] = sepExpr(expression2)

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

saveRDS(exp, file = 'Data/mouseExpr.rds')
saveRDS(mouseDes, file = 'Data/mouseDes.rds')
saveRDS(gene, file = 'Data/gene.rds')

saveRDS(exp2, file = 'Data/mouseExpr2.rds')
saveRDS(mouseDes2, file = 'Data/mouseDes2.rds')
saveRDS(gene2 ,file = 'Data/gene2.rds')
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
sourceGithub(OganM,brainGenesManuscript,'R/regionize.R')
prop ='ShinyNames'
regionNames = 'Region'

#mouseDes = read.design('Data/meltedDesign.tsv')
mouseDes = mouseDes[match(colnames(exp),make.names(mouseDes$sampleName)),]
regionGroups = memoReg(mouseDes,regionNames,prop,
                       regionHierarchy = regionHierarchy
)

#mouseDes2 = read.design('Data/meltedDesign2.tsv')
mouseDes2 = mouseDes2[match(colnames(exp2),make.names(mouseDes2$sampleName)),]

regionGroups2 = memoReg(mouseDes2,regionNames,prop, 
                        regionHierarchy= regionHierarchy
)

save(memoReg,regionHierarchy,file = 'memoReg.rda')
purge()
