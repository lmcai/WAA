
#awk -v OFS='\t' '{print $1,$2}' WAA1ReadsPerGene.out.tab >WAA1.edgeR.gene.count
#awk -v OFS='\t' '{print $1,$2}' WAA2ReadsPerGene.out.tab >WAA2.edgeR.gene.count
#awk -v OFS='\t' '{print $1,$2}' WAA3ReadsPerGene.out.tab >WAA3.edgeR.gene.count
#awk -v OFS='\t' '{print $1,$2}' WAA_SG1ReadsPerGene.out.tab >WAA_SG1.edgeR.gene.count
#awk -v OFS='\t' '{print $1,$2}' WAA_SG2ReadsPerGene.out.tab >WAA_SG2.edgeR.gene.count
#awk -v OFS='\t' '{print $1,$2}' WAA_SG3ReadsPerGene.out.tab >WAA_SG3.edgeR.gene.count

library(edgeR)

#######################
#salivary glands vs whoel body tissue

files <- dir(pattern="*\\.count$")
RG <- readDGE(files)
group=factor(c(1,1,1,2,2,2))
y <- DGEList(counts=RG,group=group)

#filter by counts
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)

#visualization
plotMDS(y)
plotBCV(y)

#glm test
design <- model.matrix(~0+group, data=y$samples)
y <- estimateDisp(y,design)
fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit,contrast=c(-1,1))
topTags(qlf)
summary(decideTests(qlf))
       group2
Down     2116
NotSig   8648
Up       2870
#plotMD(qlf)

#design <- model.matrix(~group, data=y$samples)
#fit <- glmQLFit(y,design)
#qlf <- glmQLFTest(fit,coef=2)
#topTags(qlf)

lrt <- glmLRT(fit,contrast=c(-1,1))
topTags(lrt)
summary(decideTests(lrt))
       -1*group1 1*group2
Down                 2019
NotSig               8498
Up                   3117

###
#output
write.table(qlf$table,'WAA_SG_qlf.tsv')
write.table(lrt$table,'WAA_SG_lrt.tsv')
write.table(decideTests(qlf),'WAA_SG_qlf_updownregulation.tsv')

#######################
#host plant resistance 
files <- dir(pattern="EL*")
group=factor(c(1,1,1,2,2,2,3,3,3))
RG <- readDGE(files)
y <- DGEList(counts=RG,group=group)

#filter by counts
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
design <- model.matrix(~0+group, data=y$samples)
y <- estimateDisp(y,design)

#visualization
plotMDS(y)
plotBCV(y)

#glm test

y <- estimateDisp(y,design)
fit <- glmQLFit(y,design)
qlf_16_87 <- glmQLFTest(fit,contrast=c(-1,1,0))
summary(decideTests(qlf_16_87))
none of the pairwise comparisons result in significant result




