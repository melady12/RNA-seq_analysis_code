library("statmod")
library("edgeR")
library(reshape2)
library(plyr)
#args<-commandArgs(TRUE)
#inp1<-args[1]
#inp2<-args[2]
#oup<-args[3]

#setwd('/Users/outfit_zjm/Desktop/Chicken.Broiler.vs. Layer.study/new.RANseq.analysis/02.whole.samples.PCG.exp.DE.between.breed/edgeR.result/')

#tissue_list<-c('Bra', 'Hea', 'Hyp', 
#               'Kid', 'Liv', 'Lun', 
#               'Pit','SAT', 'Spl', 'VAT')
#
#tissue_list<-c('Hyp')
#for (tissue in tissue_list){
#  
#  inp<-paste(tissue,'.samples.counts.matrix.filt.txt',sep='')
#  rawdata<-read.table(file=inp,sep='\t',row.names = 1,header=T)
#  y <- DGEList(counts=rawdata, genes=row.names(rawdata))
#  y <- calcNormFactors(y,method="TMM")
  
  #inp2<-'F1.POF.forDEG.info.txt'
  #data2<-read.table(inp2,head=T,sep='\t',stringsAsFactors = F)
#  order_list<-c(rep(1,time=6),rep(2,time=6))
#  myfactors<-factor(order_list)
#  data.frame(Sample=colnames(y), myfactors)
#  design<-model.matrix(~myfactors)
#  rownames(design)<-colnames(y)
#  y<-estimateDisp(y, design, robust=TRUE)
#  fit <- glmQLFit(y,design)
#  myOP <- glmQLFTest(fit)
#  deg<-topTags(myOP, n=nrow(myOP$table), adjust.method="BH")$table
# write.table(deg, file=paste(tissue,'.PCG.edgeR.result.WB.LH.txt',sep=''), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
#
#}


gp<-read.table("Sample.gp",stringsAsFactors=F,sep="\t")
dat<-read.table("GSE107593_raw_reads_BCHRNAseq.mody.txt",stringsAsFactors=F,sep="\t",header=TRUE)
inp<-dat[,!(names(dat)%in%c('Row','Type','gene_type','chr','start','end','strand','length'))]
mdat<-melt(inp,id.var='gene_name',stringsAsFactors=F)
out<-ddply(mdat,.(gene_name,variable),summarize,val=sum(value))
dd<-dcast(out,gene_name~variable,value.var='val')
gname<-dd$gene_name
dd<-dd[,-1]
rownames(dd)<-gname

dd<-dd[,order(names(dd))]
gp<-gp[order(gp$V1),]

y <- DGEList(counts=dd, genes=row.names(dd))
y <- calcNormFactors(y,method="TMM")
myfactors<-factor(gp$V2)
data.frame(Sample=colnames(y), myfactors)
design<-model.matrix(~myfactors)
rownames(design)<-colnames(y)
y<-estimateDisp(y, design, robust=TRUE)
fit <- glmQLFit(y,design)
myOP <- glmQLFTest(fit)
deg<-topTags(myOP, n=nrow(myOP$table), adjust.method="BH")$table
write.table(deg, file="edgeR.de.genes.xls", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
