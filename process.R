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
# mouseExpr =  fread('Data/mouseExpr')
# mouseGene = read.csv('Data/mouseGene')
# mouseDes = read.design('Data/meltedDesign.tsv')
# 
# mouseExpr2 =  fread('Data/mouseExpr2')
# mouseGene2 = read.csv('Data/mouseGene2')
# mouseDes2 = read.design('Data/meltedDesign2.tsv')
# 
# 
# tempGene2 = mouseGene2[match(mouseGene$Gene.Symbol,mouseGene2$Gene.Symbol),]
# tempExp2 = mouseExpr2[match(mouseGene$Gene.Symbol,mouseGene2$Gene.Symbol),]
#     
# tempGene1 = mouseGene[!is.na(tempGene2$Gene.Symbol),]
# tempExp1 = mouseExpr[!is.na(tempGene2$Gene.Symbol),]
# 
# tempExp2 = tempExp2[!is.na(tempGene2$Gene.Symbol),]
# tempGene2 = tempGene2[!is.na(tempGene2$Gene.Symbol),]
# 
# sameProbes = tempGene1$Probe[as.char(tempGene2$Probe) == as.char(tempGene1$Probe)]
# 
# remove = (mouseGene2$Probe %in% sameProbes)
# 
# mouseGene2 = mouseGene2[!remove,]
# mouseExpr2 = mouseExpr2[!remove,]
# write.table(mouseGene2,'Data/mouseGene22',sep=',',row.names=F)
# write.table(mouseExpr2,'Data/mouseExpr22',sep=',',row.names=F,quote=F)
# 

