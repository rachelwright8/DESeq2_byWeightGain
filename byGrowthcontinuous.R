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
  design = ~ means)


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
means = apply(countData,1,mean)
table(means>2)
# FALSE  TRUE 
# 45427 17132 

means2 = names(means[means>2])
head(means2)
length(means2)
#17132

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
dds <- DESeqDataSetFromMatrix(
   countData = countDataFilt,
   colData = colDataGro,
   design = ~ mean)
 
#############   #############   #############   #############   #############
#############              Explore data                       #############
#############   #############   #############   #############   #############

colMeans(counts(dds))

apply(counts(dds), 2, quantile, 0:10/10)

#############   #############   #############   #############   #############
#############          DESeq                                    #############
#############   #############   #############   #############   #############

#-------------DESeq pipeline in one step: makes large DESeqDataSet
deds <- DESeq(dds)
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

#############   #############   #############   #############   #############
#############          PCA and Heatmaps                         #############
#############   #############   #############   #############   #############

##---------------Sample distance heatmap
rld = rlogTransformation(dds)
head(assay(rld))
hist(assay(rld))
colnames(rld)=names(countDataFilt)
head(assay(rld))

#############   #############   #############   #############   #############
#############          Extract results                          #############
#############   #############   #############   #############   #############

res <- results( deds, independentFiltering = F)
res
#log2 fold change (MAP): mean 

table(res$padj < 0.05)
# FALSE  TRUE 
#16912    82 
table( is.na(res$padj) )
# FALSE   
# 16994    

write.csv( as.data.frame(res), file="results_Growth.csv" )

#save(colData, countDataFilt, dds, deds, res, rld, file="results.Rdata")

load("results.Rdata")

#############   #############   #############   #############   #############
#############          Diagnostics                              #############
#############   #############   #############   #############   #############

#####-------------Dispersions plot
quartz()
plotDispEsts(deds, main="Dispersion Plot Baseline")


####-----------MA plot
quartz()
plotMA(res, ylim = c(-0.08, 0.08) ) 


#############   #############   #############   #############   #############
#############          Write results for heatmaps               #############
#############   #############   #############   #############   #############

###--------------Get pvals
head(res)
vals = cbind(res$pvalue, res$padj)
head(vals)
colnames(vals) = c("pval", "padj")
length(vals[,1])
table(complete.cases(vals))

######-------------Make rlogdata and pvals table
rldd = assay(rld)
head(rldd)
colnames(rldd)=names(countDataFilt2)
rldpvals = cbind(rldd,vals)
head(rldpvals)
dim(rldpvals)
#16994    46
table(complete.cases(rldpvals))
#   TRUE 
#  16994 

write.csv(rldpvals, "9feb15_Mcav_bm2_growth_RLDandPVALS.csv", quote=F)

#############   #############   #############   #############   #############
#############             Write results for GO                  #############
#############   #############   #############   #############   #############

head(res)
logs = data.frame(cbind("gene" = row.names(res),"logP" = round(-log(res$pvalue+1e-10,10),1)))
logs$logP = as.numeric(as.character(logs$logP))
sign = rep(1,nrow(logs))
sign[res$log2FoldChange<0] = -1  ##change to correct model
table(sign)
# -1     1 
# 8436 8558 
logs$logP = logs$logP*sign
write.table(logs,quote=F,row.names=F,file="GO_Growth_logP.csv",sep=",")



##### Trying something new: multiply p-value by LFC

#############   #############   #############   #############   #############
#############             Write results for KOG                  #############
#############   #############   #############   #############   #############

head(res)
logs = data.frame(cbind("gene" = row.names(res),"logP" = round(-log(res$pvalue+1e-10,10),1)))
logs$logP = as.numeric(as.character(logs$logP))
logs$logP = logs$logP*res$log2FoldChange
head(logs)
write.table(logs,quote=F,row.names=F,file="GO_GrowthLFC_logP.csv",sep=",")
