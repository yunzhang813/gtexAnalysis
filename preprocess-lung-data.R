## preprocess and organize GTEx lung data 
load(file="gtex-gene-counts-lung.rda")

## add subject specific data
subjtab <- read.table(file="GTEx_Data_V4_Annotations_SubjectPhenotypes_DS.txt", fill=T,
                      stringsAsFactors=FALSE, sep="\t", header=T)
subjIDs <- sapply(strsplit(colnames(lungdat),"-"),function(x) paste(x[1],x[2],sep="-"))
map <- match(subjIDs,subjtab$SUBJID)
identical(subjIDs,subjtab$SUBJID[map])
ctab <- data.frame(subjtab[map,], lungtab)

## filter genes with low counts
hist(log10(rowSums(lungdat)))
abline(v=log10(5000), col=2)
i.rm <- which(rowSums(lungdat)<5000)

### option 1: remove low counts genes and then rlog() ###
lungdat1 <- lungdat[-i.rm,]
gtab1 <- gtab[-i.rm,]
## transform data to a DESeq object
library(DESeq2)
lungDES1 <- DESeqDataSetFromMatrix(countData = lungdat1,
                                  colData = ctab,
                                  design = ~ 1)
rownames(lungDES1) <- gtab1$Name
## rlog transform
rld1 <- rlog(lungDES1, fast=TRUE)

### option 2: rlog() then remove low counts genes ###
lungDES <- DESeqDataSetFromMatrix(countData = lungdat,
                                   colData = ctab,
                                   design = ~ 1)
rownames(lungDES) <- gtab$Name
rld <- rlog(lungDES, fast=TRUE)
rld2 <- rld[-i.rm,]

### compare option 1 vs. option 2 ###
smoothScatter(x=assay(rld1), y=assay(rld2))
smoothScatter(x=(assay(rld1)+assay(rld2))/2, y=assay(rld2)-assay(rld1))

### option 3: rlog() then filt out low expressed genes based on rld ###
i.keep <- which(rowMeans(assay(rld))>5)
rld3 <- rld[i.keep,]
gtab3 <- gtab[i.keep,]

####################################################
## choose to use option 3 and save data
rld <- rld3
gtab <- gtab3
save(rld, gtab, file="rlog-lung.rda")
####################################################
