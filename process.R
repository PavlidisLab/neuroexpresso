library(dplyr)
library(magrittr)
#library(ogbox)
#mouseDes2 = read.design('/home/omancarci/wholeOtto/omancarci/brainCellTypeSpecificGenes/data/meltedDesign2.tsv')
#expression2 = read.exp('/home/omancarci/wholeOtto/omancarci/brainCellTypeSpecificGenes/data/finalExp2.csv')

#mouseDes = read.design('/home/omancarci/wholeOtto/omancarci/brainCellTypeSpecificGenes/data/meltedDesign.tsv')
#expression = read.exp('/home/omancarci/wholeOtto/omancarci/brainCellTypeSpecificGenes/data/finalExp.csv')


list[gene,exp] = sepExpr(expression)
list[gene2,exp2] = sepExpr(expression2)

exp = exp[!grepl('[|]',gene$Gene.Symbol),]
gene = gene[!grepl('[|]',gene$Gene.Symbol),]

exp2 = exp2[!grepl('[|]',gene2$Gene.Symbol),]
gene2 = gene2[!grepl('[|]',gene2$Gene.Symbol),]


gene %<>% select(Probe, Gene.Symbol,GemmaIDs,GeneNames,NCBIids)
gene2 %<>% select(Probe, Gene.Symbol, GemmaIDs,GeneNames,NCBIids)

write.table(gene,'Data/mouseGene',sep=',',row.names=F)
write.table(format(exp,digits=3),'Data/mouseExpr',sep=',',row.names=F,quote=F)

write.design(mouseDes,'Data/meltedDesign.tsv')
write.design(mouseDes2,'Data/meltedDesign2.tsv')


write.table(gene2,'Data/mouseGene2',sep=',',row.names=F)
write.table(format(exp2,digits=3),'Data/mouseExpr2',sep=',',row.names=F,quote=F)

file.remove('Data/finalExp.csv')
file.remove('Data/finalExp2.csv')

purge()
