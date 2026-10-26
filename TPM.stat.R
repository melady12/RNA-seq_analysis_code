nam<-read.table("Sample.gp",stringsAsFactors=F,sep="\t")
nam<-nam[order(nam$V2),]
nam$num<-c(1:sum(nam$V2=="Inflamed"),1:sum(nam$V2=="Non-Inflamed"))
nam$ID<-paste(nam$V2,nam$num,sep="_")
dat<-read.table("GSE107593_raw_reads_BCHRNAseq.mody.txt",stringsAsFactors=F,sep="\t",header=TRUE)
for (i in 10:(dim(dat)[2])){
	vv<-dat[,i]/dat[,9]
	ww<-sum(vv)
	tt<-vv*1000000/ww
	dat[,i]<-tt
	cid<-nam$ID[nam$V1==names(dat)[i]]
	names(dat)[i]<-cid
}
write.table(dat,"whole.genes.TPM.xls",quote=F,sep="\t",row.names=F)
