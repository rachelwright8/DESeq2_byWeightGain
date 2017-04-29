library("DESeq2") 
library("arrayQualityMetrics")
library("Biobase") #for ExpressionSet
setwd("~/Documents/tsa_2014/RNASeq/byGrowth/")

#############   #############   #############   #############   #############
#############              Load data                            #############
#############   #############   #############   #############   #############

#-----------Load counts table
countdata = read.table("../allcounts_Mcav_august2015.txt",header=TRUE,row.names=1) 
head(countdata) 
length(countdata[,1])
# 62559 genes mapped

# make sample names prettier
names(countdata)
names(countdata) = sub("*.fq.trim.sam.counts","",names(countdata))
names(countdata) = sub("MC","", names(countdata))

#-----------Load design table
coldata = read.csv("../conds.csv")
head(coldata)
# sam = sample name; bank = east or west bank; geno = genotype identifier; tod = time of death in days; early = 1 for died early, 0 for did not die early; pheno = susceptible or resistant

# genotype is specified as a number but is a factor... make it so
coldata$geno = as.factor(coldata$geno)
summary(coldata)


# add growth data
growth = read.csv("../../buoyantWeights/allGrowthData.csv")
head(growth)
growth = growth[c(1,12)]

colDataGro = merge(coldata, growth, by="geno", all=T)

# get rid of NA samples
colDataGro = colDataGro[!(is.na(colDataGro$mean)),]
colDataGro

# ---- I tried using growth as a continuous measure, but the results were not meaningful. Try binning growth:

hist(colDataGro$mean, 20)

# start with three bins
colDataGro$bin3 <- cut(colDataGro$mean, c(0, 10, 15, 30), labels = 1:3)
summary(colDataGro$bin3)

# four bins
colDataGro$bin4 <- cut(colDataGro$mean, c(0, 5, 10, 15, 30), labels = 1:4)
summary(colDataGro$bin4)

# ---- remove NA for growth from col and count
countData = countdata[,names(countdata) %in% colDataGro$sam]
head(countData)

#-----------Match order of sample names and conds
colDataGro$sam
names(countData)
colDataGro = colDataGro[match(names(countData), colDataGro$sam),]

#-------check sample order
test = cbind(names(countData),as.vector(colDataGro$sam))
test
# looks good


#############   #############   #############   #############   #############
#############              Construct data object                #############
#############   #############   #############   #############   #############

dds <- DESeqDataSetFromMatrix(
  countData = countData,
  colData = colDataGro,
  design = ~ bin3)


#############   #############   #############   #############   #############
#############              Call outliers                       #############
#############   #############   #############   #############   #############

# rld = rlogTransformation(dds, blind=TRUE)
# e = ExpressionSet(assay(rld), AnnotatedDataFrame(as.data.frame(colData(rld))))
# arrayQualityMetrics(e, intgroup=c("mean"), force=T)

## No sample failed more than one test... leave them all in for now

##----Failed tests
# fail1=c("17A","17B","21B","3A","4B")
# 
# #--Remove outliers
# countdataOut=countdata[!(names(countdata) %in% fail1)]
# head(countdataOut)
# 
# coldataOut=coldata[!(coldata$sam %in% fail1),]
# head(coldataOut)
# 
# save(coldataOut, countdataOut, rld, e, file="data4DESeq2_outRemoved.Rdata")

#---------Set base mean minimum
means = apply(countdata,1,mean)
table(means>2)
# FALSE  TRUE 
# 45565 16994 

means2 = names(means[means>2])
head(means2)
length(means2)
#16994

countDataFilt = countData[row.names(countData) %in% means2,]
head(countDataFilt)

totalCountsFilt = colSums(countDataFilt)
totalCountsFilt

min(totalCountsFilt) #123396
max(totalCountsFilt) #926410
mean(totalCountsFilt) #556581.5


# #############   #############   #############   #############   #############
# #############    Reconstruct data object  (outliers removed)    #############
# #############   #############   #############   #############   #############
# 
dds3 <- DESeqDataSetFromMatrix(
   countData = countDataFilt,
   colData = colDataGro,
   design = ~ bin3)
 
dds4 <- DESeqDataSetFromMatrix(
  countData = countDataFilt,
  colData = colDataGro,
  design = ~ bin4)

#############   #############   #############   #############   #############
#############              Explore data                       #############
#############   #############   #############   #############   #############

colMeans(counts(dds3))
colMeans(counts(dds4))

apply(counts(dds3), 2, quantile, 0:10/10)
apply(counts(dds4), 2, quantile, 0:10/10)

#############   #############   #############   #############   #############
#############          DESeq                                    #############
#############   #############   #############   #############   #############

#-------------DESeq pipeline in one step: makes large DESeqDataSet
deds3 <- DESeq(dds3)
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing
# -- replacing outliers and refitting for 87 genes
# -- DESeq argument 'minReplicatesForReplace' = 7 
# -- original counts are preserved in counts(dds)
# estimating dispersions
# fitting model and testing

deds4 <- DESeq(dds4)
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing
# -- replacing outliers and refitting for 68 genes
# -- DESeq argument 'minReplicatesForReplace' = 7 
# -- original counts are preserved in counts(dds)
# estimating dispersions
# fitting model and testing


#############   #############   #############   #############   #############
#############          PCA and Heatmaps                         #############
#############   #############   #############   #############   #############

##---------------Sample distance heatmap
rld3 = rlogTransformation(dds3)
head(assay(rld3))
hist(assay(rld3))
colnames(rld3)=names(countDataFilt)
head(assay(rld3))

rld4 = rlogTransformation(dds4)
head(assay(rld4))
hist(assay(rld4))
colnames(rld4)=names(countDataFilt)
head(assay(rld4))

#############   #############   #############   #############   #############
#############          Extract results                          #############
#############   #############   #############   #############   #############

res3 <- results( deds3, independentFiltering = F)
res3
table(res3$padj < 0.05)
# FALSE  TRUE 
# 16956    38 
table( is.na(res3$padj) )
# FALSE   
# 16994    

write.csv( as.data.frame(res3), file="results_bin3Growth.csv" )

res4 <- results( deds4, independentFiltering = F)
res4
table(res4$padj < 0.05)
# FALSE  TRUE 
#16846   141 
table( is.na(res4$padj) )
# FALSE   
# 16994    

write.csv( as.data.frame(res4), file="results_bin4Growth.csv" )

save(colData, countDataFilt, dds3, dds4, deds3, deds4, res3,res4, rld3, rld4, file="results_bins.Rdata")

#load("results_bins.Rdata")

#############   #############   #############   #############   #############
#############          Diagnostics                              #############
#############   #############   #############   #############   #############

#####-------------Dispersions plot
quartz()
plotDispEsts(deds3, main="Dispersion Plot Baseline")

quartz()
plotDispEsts(deds4, main="Dispersion Plot Baseline")

####-----------MA plot
quartz()
par(mfrow=c(1,2))
plotMA(res3, ylim = c(-1, 1) ) 
plotMA(res4, ylim = c(-1, 1) ) 



#############   #############   #############   #############   #############
#############          Write results for heatmaps               #############
#############   #############   #############   #############   #############

###--------------Get pvals
head(res3)
vals = cbind(res3$pvalue, res3$padj)
head(vals)
colnames(vals) = c("pval", "padj")
length(vals[,1])
table(complete.cases(vals))

######-------------Make rlogdata and pvals table
rldd = assay(rld3)
rldpvals = cbind(rldd,vals)
head(rldpvals)
dim(rldpvals)
#16994    46
table(complete.cases(rldpvals))
#   TRUE 
#  16994 

write.csv(rldpvals, "9feb15_Mcav_bm2_growthbin3_RLDandPVALS.csv", quote=F)

head(res4)
vals = cbind(res4$pvalue, res4$padj)
head(vals)
colnames(vals) = c("pval", "padj")
length(vals[,1])
table(complete.cases(vals))

######-------------Make rlogdata and pvals table
rldd = assay(rld4)
rldpvals = cbind(rldd,vals)
head(rldpvals)
dim(rldpvals)
#16994    46
table(complete.cases(rldpvals))
#   TRUE 
#  16994 

write.csv(rldpvals, "9feb15_Mcav_bm2_growthbin4_RLDandPVALS.csv", quote=F)

#############   #############   #############   #############   #############
#############             Write results for GO                  #############
#############   #############   #############   #############   #############

head(res3)
logs = data.frame(cbind("gene" = row.names(res3),"logP" = round(-log(res3$pvalue+1e-10,10),1)))
logs$logP = as.numeric(as.character(logs$logP))
sign = rep(1,nrow(logs))
sign[res3$log2FoldChange<0] = -1  ##change to correct model
table(sign)
# -1     1 
# 8358 8636 
logs$logP = logs$logP*sign
write.table(logs,quote=F,row.names=F,file="GO_Growthbin3_logP.csv",sep=",")

head(res4)
logs = data.frame(cbind("gene" = row.names(res4),"logP" = round(-log(res4$pvalue+1e-10,10),1)))
logs$logP = as.numeric(as.character(logs$logP))
sign = rep(1,nrow(logs))
sign[res3$log2FoldChange<0] = -1  ##change to correct model
table(sign)
# -1     1 
# 8358 8636 
logs$logP = logs$logP*sign
write.table(logs,quote=F,row.names=F,file="GO_Growthbin4_logP.csv",sep=",")
