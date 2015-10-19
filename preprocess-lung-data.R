## preprocess and organize GTEx lung data 
load(file="gtex-gene-counts-lung.rda")

## filter genes with low counts
hist(log10(rowSums(lungdat)))
abline(v=log10(5000))
i.rm <- which(rowSums(lungdat)<5000)
lungdat <- lungdat[-i.rm,]
gtab <- gtab[-i.rm,]

## add subject specific data
subjtab <- read.table(file="GTEx_Data_V4_Annotations_SubjectPhenotypes_DS.txt", fill=T,
                      stringsAsFactors=FALSE, sep="\t", header=T)
subjIDs <- sapply(strsplit(colnames(lungdat),"-"),function(x) paste(x[1],x[2],sep="-"))
map <- match(subjIDs,subjtab$SUBJID)
identical(subjIDs,subjtab$SUBJID[map])
ctab <- data.frame(subjtab[map,], lungtab)

## transform data to a DESeq object
library(DESeq2)
lungDES <- DESeqDataSetFromMatrix(countData = lungdat,
                                  colData = ctab,
                                  design = ~ 1)
rownames(lungDES) <- gtab$Name

## rlog transform
rld <- rlog(lungDES, fast=TRUE)

save(rld, gtab, file="rlog-lung.rda")
