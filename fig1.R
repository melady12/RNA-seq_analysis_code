dat<-read.table("edgeR.de.genes.include.genetype.xls",stringsAsFactors=F,header=TRUE)
dat$logFC<--dat$logFC
dat$P<--log(dat$FDR,10)
c1<-log(1.5,2)
c2<--log(0.001,10)
col1<-"#69C034"
col2<-"#F72130"
dat$cols<-ifelse(dat$logFC<=-c1&dat$P>=c2,col1,ifelse(dat$logFC>=c1&dat$P>=c2,col2,"gray"))
dat$cols<-adjustcolor(dat$cols,alpha.f=0.6)
pdf("Figure1.pdf",w=4,h=4)
par(xaxs='i',yaxs='i',mgp=c(1.8,0.5,0),lend=2,tcl=-0.2,mar=c(3,4,5,2))
plot(1:10,axes=F,xlim=c(-10,10),ylim=c(-0.2,17),xlab=expression(log[2]~(fold~change)),ylab=expression(-log[10]~(italic(P)~value)),type='n',main="Inflamed .vs. Non-Inflamed")

points(dat$logFC,dat$P,pch=16,col=dat$cols,cex=0.6,xpd=NA)
axis(1)
axis(2,las=2)
box(bty='l')
d1<-dat[dat$logFC<=-c1&dat$P>=c2,]
d2<-dat[dat$logFC>=c1&dat$P>=c2,]
rect(xleft=min(d1$logFC),ybottom=min(d1$P),xright=max(d1$logFC),ytop=max(d1$P),col=adjustcolor(col1,alpha.f=0.1),lty=3,lwd=0.5,border=adjustcolor(col1,alpha.f=0.6),lend=3)
rect(xleft=min(d2$logFC),ybottom=min(d2$P),xright=max(d2$logFC),ytop=max(d2$P),col=adjustcolor(col2,alpha.f=0.1),lty=3,lwd=0.5,border=adjustcolor(col2,alpha.f=0.6),lend=3)
length(unique(dat$genes[dat$FDR<=0.001&dat$logFC<=-log(1.5,2)&dat$gene_type=="protein_coding"]))
length(unique(dat$genes[dat$FDR<=0.001&dat$logFC>=log(1.5,2)&dat$gene_type=="protein_coding"]))
legend(x=-11,y=19.3,legend=c('1666 down-regulated','3216 up-regulated'),nc=2,bty='n',col=c(col1,col2),xpd=NA,cex=0.8,pch=20,pt.cex=2)
dev.off()
