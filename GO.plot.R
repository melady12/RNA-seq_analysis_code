library(scales)
library(Hmisc)
dat<-read.table("GO.txt",stringsAsFactors=F,sep="\t")
vcol<-c("#2766A3","#D8137E")
#dat$cols<-ifelse(dat$V4=="Up",vcol[1],ifelse(dat$V4=="Down",vcol[2],NA))
dat$cols<-col_numeric("YlOrRd",dat$V3)(dat$V3)
dat$type<-ifelse(dat$V4=="Up",21,22)
dat$yy<-1:length(dat$V1)

dat$si<-log(dat$V3,2)



pdf("metascape.GO.pdf",w=6,h=5)
par(xaxs='i',yaxs='i',mar=c(4,16,3,3),tcl=-0.25,mgp=c(1.7,0.5,0))
plot(1:10,type='n',axes=F,xlab=expression(-Log[10]~(italic(P)~~value)),ylab="",xlim=c(5,100),ylim=rev(c(0,length(dat$V1)+1)))
grid(lty=3,col='gray',lwd=0.5)
#points(x=abs(dat$V2),y=dat$yy,pch=19,bg=dat$cols,col=dat$cols,cex=dat$si,xpd=NA)
points(x=abs(dat$V2),y=dat$yy,pch=dat$type,bg=dat$cols,col='black',cex=1.3,xpd=NA,lwd=0.1)
box(bty='l')
axis(1)
axis(2,at=dat$yy,labels=NA,xpd=NA,las=2,cex.axis=0.75)
#text(rep(0,12),dat$yy,labels=dat$V1,cex=0.75,col=c('red','red','black','black','black','black','red','red','black','black','black','black'),pos=2,xpd=NA)
text(rep(5,33),dat$yy,labels=capitalize(dat$V1),cex=0.7,col='black',pos=2,xpd=NA)

#legend(x=-3,y=-1,legend=c('Captive','Wild','Captive and Wild'),nc=3,xpd=NA,col=vcol,pch=19,bty='n',x.intersp=0,y.intersp=0,cex=1)

fleg<-function(x,y,col,size,lab,pch=19){
#par(new=TRUE,fig=c(x,x+0.5,y,y+0.5),mar=c(0,0,0,0))
#plot(1,1,pch=19,col=col,axes=F,xlab='',ylab='',bty='n',size=size)
	points(x,y,pch=pch,col='black',bg='gray',cex=size,xpd=NA,lwd=0.1)
	text(x,y,pos=4,labels=lab,xpd=NA,cex=0.8)
}

fleg(-5,-1.5,vcol[1],1.2,"down-regulated",22)
fleg(60,-1.5,vcol[2],1.2,"up-regulated",21)
#fleg(2,-0.5,vcol[3],1.2,"Captvie and wild")
#text(1,-1.2,labels="Selected regions",xpd=NA,cex=0.95)
text(-40,35.5,labels="GO biological processes",xpd=NA,cex=1)



#leg<-unique(dat[,c('V3','si')])
#leg<-leg[order(leg$V3),]
#leg$xx<-seq(-12,-5,len=length(leg$V3))
#for(i in 1:length(leg$xx)){
#	fleg(leg$xx[i],-0.5,'black',leg$si[i],"")
#	text(leg$xx[i],-0.2,pos=1,labels=leg$V3[i],xpd=NA)
#}
#text(-8,-1.3,labels="Number of genes",cex=1,xpd=NA)

par(fig=c(0.1,0.4,0.92,0.95),new=TRUE,mar=c(0,0,0,0),mgp=c(0.5,-0.1,0),cex=0.7)
imy<-seq(from=range(dat$V3,na.rm=T)[1],to=range(dat$V3,na.rm=T)[2],len=5000)
image(x=imy,y=1:2,z=t(t(imy)),col=col_numeric("YlOrRd",imy)(imy),axes=F,xlab='',ylab='',xpd=NA)
axis(1,tck=-0.15)
mtext(3,text="Number of genes",cex=0.9)
box(bty='o')



dev.off()
