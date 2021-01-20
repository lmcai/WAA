
#awk -v OFS='\t' '{print $1,$2}' WAA1ReadsPerGene.out.tab >WAA1.edgeR.gene.count
#awk -v OFS='\t' '{print $1,$2}' WAA2ReadsPerGene.out.tab >WAA2.edgeR.gene.count
#awk -v OFS='\t' '{print $1,$2}' WAA3ReadsPerGene.out.tab >WAA3.edgeR.gene.count
#awk -v OFS='\t' '{print $1,$2}' WAA_SG1ReadsPerGene.out.tab >WAA_SG1.edgeR.gene.count
#awk -v OFS='\t' '{print $1,$2}' WAA_SG2ReadsPerGene.out.tab >WAA_SG2.edgeR.gene.count
#awk -v OFS='\t' '{print $1,$2}' WAA_SG3ReadsPerGene.out.tab >WAA_SG3.edgeR.gene.count

library(edgeR)

#######################
#salivary glands vs whole body tissue

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

#####################################
#Malus analysis

#prefiltered genes (significantly expressed in the experiemnt set but not the control set)
setwd('malus_sig_exp_genes/')
files =c("1p.edgeR.gene.count","25p.edgeR.gene.count","37p.edgeR.gene.count","40p.edgeR.gene.count","46p.edgeR.gene.count","48p.edgeR.gene.count","4p.edgeR.gene.count","15p.edgeR.gene.count","16p.edgeR.gene.count","20p.edgeR.gene.count","27p.edgeR.gene.count","2p.edgeR.gene.count","33p.edgeR.gene.count","36p.edgeR.gene.count","39p.edgeR.gene.count","42p.edgeR.gene.count","6p.edgeR.gene.count","13p.edgeR.gene.count","19p.edgeR.gene.count","23p.edgeR.gene.count","32p.edgeR.gene.count")
group=factor(c(rep(1,7),rep(2,10),rep(3,4)))
RG <- readDGE(files)
y <- DGEList(counts=RG,group=group)
#plotMDS(y,cex=0.3,col=c(rep("black",7), rep("red",10),rep("green",4)) )
#filter by counts
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)

plotMDS(y,cex=0.3,col=c(rep("black",7), rep("red",10),rep("green",4)) )


#add control
setwd('..')
files =c("1p.edgeR.gene.count","25p.edgeR.gene.count","37p.edgeR.gene.count","40p.edgeR.gene.count","46p.edgeR.gene.count","48p.edgeR.gene.count","4p.edgeR.gene.count","15p.edgeR.gene.count","16p.edgeR.gene.count","20p.edgeR.gene.count","27p.edgeR.gene.count","2p.edgeR.gene.count","33p.edgeR.gene.count","36p.edgeR.gene.count","39p.edgeR.gene.count","42p.edgeR.gene.count","6p.edgeR.gene.count","13p.edgeR.gene.count","19p.edgeR.gene.count","23p.edgeR.gene.count","32p.edgeR.gene.count","18c.edgeR.gene.count","29c.edgeR.gene.count","43c.edgeR.gene.count","50c.edgeR.gene.count","51c.edgeR.gene.count","10c.edgeR.gene.count","22c.edgeR.gene.count","24c.edgeR.gene.count","41c.edgeR.gene.count","45c.edgeR.gene.count","26c.edgeR.gene.count","31c.edgeR.gene.count","35c.edgeR.gene.count","49c.edgeR.gene.count","5c.edgeR.gene.count")
group=factor(c(rep(1,7),rep(2,10),rep(3,4),rep(4,15)))
RG <- readDGE(files)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
plotMDS(y,cex=0.3,col=c(rep("black",7), rep("red",10),rep("green",4),rep("blue",15)))

###########
#exluding negative treatment and identify candidates
files =c("18c.edgeR.gene.count","29c.edgeR.gene.count","43c.edgeR.gene.count","50c.edgeR.gene.count","51c.edgeR.gene.count","10c.edgeR.gene.count","22c.edgeR.gene.count","24c.edgeR.gene.count","41c.edgeR.gene.count","45c.edgeR.gene.count","26c.edgeR.gene.count","31c.edgeR.gene.count","35c.edgeR.gene.count","49c.edgeR.gene.count","5c.edgeR.gene.count","1p.edgeR.gene.count","37p.edgeR.gene.count","40p.edgeR.gene.count","46p.edgeR.gene.count","4p.edgeR.gene.count","15p.edgeR.gene.count","20p.edgeR.gene.count","33p.edgeR.gene.count","36p.edgeR.gene.count","39p.edgeR.gene.count","42p.edgeR.gene.count","6p.edgeR.gene.count","19p.edgeR.gene.count","23p.edgeR.gene.count","32p.edgeR.gene.count")
group=factor(c(rep(1,15),rep(2,15)))
RG <- readDGE(files)
y <- DGEList(counts=RG,group=group)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
plotMDS(y,cex=0.7,col=c(rep("black",15), rep("red",15)))

design <- model.matrix(~0+group, data=y$samples)
y <- estimateDisp(y,design)
fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit,contrast=c(-1,1))
topTags(qlf)
summary(decideTests(qlf))
write.table(qlf$table,'malus_insectRNA_qlf.tsv')
