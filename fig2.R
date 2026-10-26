###############################################################################
#File Name: ../../Rscript/PCA.plot.R
#Author: liyan
#Mail: liyan1380@novogene.cn 
#Created Time: 2016年12月15日 星期四 14时58分54秒
################################################################################
library(RColorBrewer)
library(ellipse)
library(reshape2)
#args<-commandArgs(TRUE)
#input<-args[1]
#group<-args[2]
#outpre<-args[3]
#logi<-args[4]
group<-"Sample.gp"
outpre<-"PCA"
logi<-FALSE
dd1<-read.table("melted.data.RNA-seq.txt",stringsAsFactors=F,sep="\t")
mdd<-acast(dd1,V1~V2,value.var="V3")
pca<-prcomp(mdd,scale=FALSE)
pc1<-summary(pca)$importance[2,1]*100
pc2<-summary(pca)$importance[2,2]*100
pc3<-summary(pca)$importance[2,3]*100
input<-data.frame(Sample=rownames(pca$rotation),Sample2=rownames(pca$rotation),x=pca$rotation[,1],y=pca$rotation[,2],z=pca$rotation[,3],stringsAsFactors=F)



op1<-paste(outpre,".PC1.and.PC2.pdf",sep='')
#op2<-paste(outpre,".PC2.and.PC3.pdf",sep='')
#op3<-paste(outpre,".PC1.and.PC3.pdf",sep='')

#pc<-read.table(input,stringsAsFactors=F)
pc<-input
names(pc)<-c('V1','V2','PC1','PC2','PC3')
ginf<-read.table(group,stringsAsFactors=F,sep="\t")
names(ginf)[2]<-'pop'
pdat<-merge(pc,ginf,by='V1')
pdat<-pdat[order(pdat$pop),]
wcols<-rep(brewer.pal(9,'Set1'),times=1000)
#cols1<-colorRampPalette(brewer.pal(9,'Set1'))(length(unique(pdat$pop)))
cols<-wcols[1:length(unique(pdat$pop))]
cols1<-adjustcolor(cols,alpha.f=0.9)
#cols2<-colorRampPalette(brewer.pal(9,'Set1'))(length(unique(pdat$pop)))
cols2<-adjustcolor(cols,alpha.f=0.2)
wpch<-rep(0:14,times=1000)
pchs<-rev(wpch)[1:length(unique(pdat$pop))]
#cols1<-rainbow(length(unique(pdat$pop)),alpha=0.8)
#cols2<-rainbow(length(unique(pdat$pop)),alpha=0.1)

tab<-table(pdat$pop)
tab<-tab[order(names(tab))]

pdat$pcol<-rep(cols1,times=tab)
pdat$ecol<-rep(cols2,times=tab)
pdat$pch<-rep(pchs,times=tab)

filt<-function(mat,probs=NULL){
	dis<-dist(mat,diag=T,upper=T)
	fil<-apply(as.matrix(dis),1,mean)
	qua<-quantile(fil,probs=probs)
	fil<=qua
}

fpca<-function(c1,c2,xlab=NULL,ylab=NULL,out=NULL,...){
	pdf(out,w=3.7,h=3.7)
	par(lend=2,mgp=c(2,0.5,0),mar=c(3.2,3.2,3,2),xpd=NA,tcl=-0.2)
	inter<-pdat[,c(c1,c2,1,6,7,8,9)]
	print(head(inter))
	gv<-unique(pdat$pop)
	inmat<-as.matrix(pdat[,c(c1,c2)])
	for (i in seq(gv)){
		inter1<-inter[inter$pop==gv[i],]
		mat<-as.matrix(inter1[,1:2])
		inter1<-inter1[filt(mat,probs=0.99),]
		co<-cov(inter1[,1:2])
		me<-apply(inter1[,1:2],2,mean)
		ell<-ellipse(x=co,centre=me,npoints = 300,level=0.98)
		inmat<-rbind(inmat,ell)
		#polygon(x=ell[,1],y=ell[,2],col=inter1$ecol[1],border=inter1$ecol[1],lwd=0.1)
		}
	xlimt<-range(inmat[,1],na.rm=T)
	ylimt<-range(inmat[,2],na.rm=T)
    
	inter<-pdat[,c(c1,c2,1,6,7,8,9)]
    plot(x=pdat[,c1],y=pdat[,c2],pch=pdat$pch,lwd=1,col=inter$pcol,xlab=xlab,ylab=ylab,bty='o',cex=0.8,cex.axis=0.8,xlim=xlimt,ylim=ylimt,axes=F)
    if(logi){text(x=inter[,1],y=inter[,2],pos=3,offset=0.3,col=inter$pcol,cex=0.2,labels=inter[,3],xpd=NA)}
    gv<-unique(pdat$pop)
    for (i in seq(gv)){
        inter1<-inter[inter$pop==gv[i],]
        mat<-as.matrix(inter1[,1:2])
        inter1<-inter1[filt(mat,probs=0.99),]
        co<-cov(inter1[,1:2])
        me<-apply(inter1[,1:2],2,mean)
        ell<-ellipse(x=co,centre=me,npoints = 300,level=0.98)
		#polygon(x=ell[,1],y=ell[,2],col=inter1$ecol[1],border=inter1$ecol[1],lwd=0.1)
        }
    leg<-unique(inter[,c('pop','pcol','pch')])
    y<-quantile(ylimt,probs=0.75)
    x<-max(xlimt)+diff(xlimt)/15
	axis(1,at=pretty(par('usr')[1:2],n=4),xpd=FALSE)
	axis(2,las=2)
	box(bty='o')
    legend("bottomright",legend=leg$pop,xpd=NA,pch=leg$pch,bty='o',col=leg$pcol,ncol=length(leg$pop)%/%8+1,cex=0.9)
    dev.off()

}
fpca(3,4,paste("PC1 (",round(pc1,2)," %)",sep=""),paste("PC2 (",round(pc2,2)," %)",sep=""),op1)
#fpca(3,5,'PC1','PC3',op3)
#fpca(4,5,'PC2','PC3',op2)
###########3dplot
#leg<-unique(pdat[,c('pop','pcol','pch')])
#pdf(op4,h=8,w=10)
#par(lend=2,mgp=c(1.9,0.5,0),xpd=NA)
#scatterplot3d(pdat[,3],pdat[,4],pdat[,5],xlab="PC1", ylab="PC3",zlab="PC2",color=pdat$pcol,pch=pdat$pch,cex.symbols=1.5,mar=c(3,3,3,10)+0.1)
#legend("right", legend=leg$pop,pch=leg$pch, col=leg$pcol,bty="n",ncol=length(leg$pop)%/%8+1 ,inset=c(-0.1,1))
#legend("right", legend=leg$pop,pch=leg$pch, col=leg$pcol,bty="n",ncol=length(leg$pop)%/%8+1 )
#dev.off()
