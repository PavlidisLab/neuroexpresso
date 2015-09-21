library(gplots)
mouseDes2 = read.design('Data/meltedDesign2.tsv')
expression2 = read.exp('Data/finalExp2.csv')

mouseDes = read.design('Data/meltedDesign.tsv')
expression = read.exp('Data/finalExp.csv')

list[gene,exp] = sepExpr(expression)
list[gene2,exp2] = sepExpr(expression2)

exp = exp[!grepl('[|]',gene$Gene.Symbol),]
gene = gene[!grepl('[|]',gene$Gene.Symbol),]

exp2 = exp2[!grepl('[|]',gene2$Gene.Symbol),]
gene2 = gene2[!grepl('[|]',gene2$Gene.Symbol),]

write.table(gene,'Data/mouseGeneTemp',sep=',',row.names=F)
write.table(format(exp,digits=3),'Data/mouseExprTemp',sep=',',row.names=F,quote=F)

write.table(gene2,'Data/mouseGene2',sep=',',row.names=F)
write.table(format(exp2,digits=3),'Data/mouseExpr2',sep=',',row.names=F,quote=F)

purge()
