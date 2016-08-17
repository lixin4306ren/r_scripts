################ calcuate sliding window methylation level from files
################
sliding.count<-function(file,window.size,step.size,min.num=20,min.mc.depth=0){
require(data.table)
require(GenomicRanges)

read.table(file,colClasses = "integer")->x
chr.name=tail(unlist(strsplit(file,'/')),n=1)
chr.name=unlist(strsplit(chr.name,'\\.'))[1]

seq(1,max(x$V1)-window.size+1,step.size)->start
seq(1+window.size-1,max(x$V1),step.size)->end

x[x$V4 > 0,]->x
x[x$V3 > min.mc.depth,]->mcg
win<-GRanges(seqnames=chr.name,IRanges(start=start,end=end))
#transfer x, y to GRange objects
cg1<-GRanges(seqnames=chr.name,IRanges(start=x$V1,end=x$V1))
mcg<-GRanges(seqnames=chr.name,IRanges(start=mcg$V1,end=mcg$V1))

as.data.frame(findOverlaps(win,cg1))->tmp1
countOverlaps(win,mcg)->mcg
countOverlaps(win,cg1)->c

which(countOverlaps(win,cg1)==0)->tmp2 #regions without Cs

data.table(id=as.character(tmp1$queryHits),cov=x$V4[tmp1$subjectHits],m=x$V3[tmp1$subjectHits])->new1
setkey(new1,id)
as.integer(new1[,sum(cov),by=id]$id)->id1
data.table(id=new1[,sum(cov),by=id]$id,start=start(win)[id1],end=end(win)[id1],m=new1[,sum(m),by=id]$V1,cov=new1[,sum(cov),by=id]$V1)->tmp1
tmp1<-rbind(data.frame(tmp1),data.frame(id=as.character(tmp2),start=start(win)[tmp2],end=end(win)[tmp2],m=rep(0,length(tmp2)),cov=rep(0,length(tmp2))))
sort(as.integer(tmp1$id),index.return = T)->ind1
tmp1[ind1$ix,]->tmp1

#data.frame(rep(chr.name,nrow(tmp1)),tmp1$start,tmp1$end,tmp1$m/tmp1$cov)->out
#paste(chr.name,".bed",sep="")->outname
#write.table(out,file=outname,row.names=F,col.names=F,sep="\t",quote=F)
tmp1$mcg=mcg
tmp1$c=c

tmp1$cov<min.num ->ind
tmp1$m_level[ind]<-NA
tmp1$m_level[!ind]<-tmp1$m[!ind]/tmp1$cov[!ind]
tmp1
}


################################# calcuate sliding window methylation level for one chromsome from objects
#################################

sliding.count2<-function(x,chr.name="chr",window.size=10000,step.size=10000,min.num=20,min.mc.depth=0){
library(data.table)
library(GenomicRanges)

seq(1,max(x$V1)-window.size+1,step.size)->start
seq(1+window.size-1,max(x$V1),step.size)->end

x[x$V4 > 0,]->x
x[x$V3 > min.mc.depth,]->mcg
win<-GRanges(seqnames=chr.name,IRanges(start=start,end=end))
#transfer x, y to GRange objects
cg1<-GRanges(seqnames=chr.name,IRanges(start=x$V1,end=x$V1))
mcg<-GRanges(seqnames=chr.name,IRanges(start=mcg$V1,end=mcg$V1))

as.data.frame(findOverlaps(win,cg1))->tmp1
countOverlaps(win,mcg)->mcg
countOverlaps(win,cg1)->c

which(countOverlaps(win,cg1)==0)->tmp2 #regions without Cs

data.table(id=as.character(tmp1$queryHits),cov=x$V4[tmp1$subjectHits],m=x$V3[tmp1$subjectHits])->new1
setkey(new1,id)
as.integer(new1[,sum(cov),by=id]$id)->id1
data.table(id=new1[,sum(cov),by=id]$id,start=start(win)[id1],end=end(win)[id1],m=new1[,sum(m),by=id]$V1,cov=new1[,sum(cov),by=id]$V1)->tmp1
tmp1<-rbind(data.frame(tmp1),data.frame(id=as.character(tmp2),start=start(win)[tmp2],end=end(win)[tmp2],m=rep(0,length(tmp2)),cov=rep(0,length(tmp2))))
sort(as.integer(tmp1$id),index.return = T)->ind1
tmp1[ind1$ix,]->tmp1

tmp1$mcg=mcg
tmp1$c=c
tmp1[tmp1$cov>min.num,]->tmp1
tmp1
}

################################ plot chromorsome level methylation level
################################
plotCG<-function(all,chr.name="chr"){
window.size<-(all[[1]]$end[1]-all[[1]]$start[1]+1[1])
all[[1]]$start+window.size/2->plot.x
all[[1]]$c/window.size->plot.y
plot(plot.x,plot.y,pch=".",cex=0.1,col="black",ylab="CpG density",xlab=chr.name,type="n")
lines(smooth.spline(plot.x[which(plot.y != "NaN")], plot.y[which(plot.y != "NaN")], spar=0.3))
}
################################
plotmanyCG<-function(all,chr.name="chr"){
window.size<-(all[[1]]$end[1]-all[[1]]$start[1]+1[1])
all[[1]]$start+window.size/2->plot.x
all[[1]]$c/window.size->plot.y
plot(plot.x,plot.y,pch=".",cex=0.1,col="blue",ylab="CpG density",xlab=chr.name,type="n")

num<-length(all)
#rep(c("black","red"),each=3)->col
c("black","black","black","red","blue","green")->col
for(i in 1:num){
window.size<-(all[[i]]$end[1]-all[[i]]$start[1]+1[1])
all[[i]]$start+window.size/2->plot.x
all[[i]]$c/window.size->plot.y
lines(smooth.spline(plot.x[which(plot.y != "NaN")], plot.y[which(plot.y != "NaN")], spar=0.3),col=col[i])
}
}

####################### plot meth.slding according to region range
plot.meth.sliding<-function(x,region,win.base=5000,step.base=5000,up.base=200000,down.base=200000,win.percent=0.01,step.percent=0.01,up.percent=1,down.percent=1,min.num=20,min.mc.depth=0,frame=TRUE,main="DNAm",col="black",append=FALSE,dist="base",x.axis=TRUE,ylim=c(0,1)){
library(data.table)
library(GenomicRanges)
small.num=0.00001
	f<-function(x,a,b,up,down){
		width(x)*a->win
		width(x)*b->step
		region.start<-start(x)-width(x)*up
		region.end<-end(x)+width(x)*down
		round(seq(region.start,region.end-win+1+small.num,step))->start
		round(seq(region.start+win-1,region.end+small.num,step))->end
		GRanges(seqnames=seqnames(x),range=IRanges(start=start,end=end))
	}

x[x[,5] > 0,]->cg.x
x[x[,4] > min.mc.depth,]->mcg.x
cg1<-GRanges(seqnames=cg.x[,1],IRanges(start=cg.x[,2],end=cg.x[,2]))
mcg1<-GRanges(seqnames=mcg.x[,1],IRanges(start=mcg.x[,2],end=mcg.x[,2]))

if(dist=="percent"){
win<-f(region[1],win.percent,step.percent,up.percent,down.percent)
for(i in 2:length(region)){
	tmp.win<-f(region[i],win.percent,step.percent,up.percent,down.percent)
	c(win,tmp.win)->win
}
}
else if(dist=="base"){
win<-generate.sliding.window(region,win.percent,step.base,up.base,down.base)
}
print("starting")

as.data.frame(findOverlaps(win,cg1))->tmp1
countOverlaps(win,mcg1)->mcg
countOverlaps(win,cg1)->c

which(countOverlaps(win,cg1)==0)->tmp2 #regions without Cs

data.table(id=as.character(tmp1$queryHits),cov=cg.x[tmp1$subjectHits,5],m=cg.x[tmp1$subjectHits,4])->new1
setkey(new1,id)
as.integer(new1[,sum(cov),by=id]$id)->id1
data.table(id=new1[,sum(cov),by=id]$id,start=start(win)[id1],end=end(win)[id1],m=new1[,sum(m),by=id]$V1,cov=new1[,sum(cov),by=id]$V1)->tmp1
tmp1<-rbind(data.frame(tmp1),data.frame(id=as.character(tmp2),start=start(win)[tmp2],end=end(win)[tmp2],m=rep(0,length(tmp2)),cov=rep(0,length(tmp2))))
sort(as.integer(tmp1$id),index.return = T)->ind1
tmp1[ind1$ix,]->tmp1

tmp1$mcg=mcg
tmp1$c=c
tmp1$cov<min.num ->ind
tmp1$cov[ind]<-0
tmp1$m[ind]<-0
tmp1$meth=tmp1$m/tmp1$cov
matrix(tmp1$meth,nrow=length(region))->tbls

if(dist=="percent"){
	pos1<-1
	pos2<-up.percent/step.percent+1
	pos3<-ncol(tbls)-down.percent/step.percent
	pos4<-ncol(tbls)
axis=c(-up.percent,"0 %","100 %",down.percent)
}
else if(dist=="base"){
	pos1<-1
	pos2<-(ncol(tbls)-1/step.percent)/2+1
	pos3<-ncol(tbls)-(ncol(tbls)-1/step.percent)/2
	pos4<-ncol(tbls)

axis=c(-up.base,"0 %","100 %",down.base)
}

if (append==FALSE){
	if(frame==TRUE){
	plot(colMeans(tbls,na.rm=T),ylim=ylim,ylab="methylation level",xlab="",xaxt="n",type="l",main=main,col=col)
	}else{
	plot(colMeans(tbls,na.rm=T),ylim=ylim,ylab="methylation level",xlab="",xaxt="n",yaxt="n",type="l",main=main,col=col,axes=F)
	}
	axis(2)
}else{
	points(colMeans(tbls,na.rm=T),type="l",col=col)

}
if(x.axis==TRUE){axis(1,c(pos1,pos2,pos3,pos4),axis)}
list(tbls=tbls,x.axis.pos=c(pos1,pos2,pos3,pos4),x.axis=axis)->plot.out
return(plot.out)
}

################ for bsseq object
plot.meth.sliding.bsseq<-function(BS.BM,region,win.base=5000,step.base=5000,up.base=200000,down.base=200000,win.percent=0.01,step.percent=0.01,up.percent=1,down.percent=1,min.num=20,min.mc.depth=0,frame=TRUE,main="DNAm",col="black",append=FALSE,dist="base",x.axis=TRUE,ylim=c(0,1)){
require(data.table)
require(GenomicRanges)
require(bsseq)
small.num=0.00001
	f<-function(x,a,b,up,down){
		width(x)*a->win
		width(x)*b->step
		region.start<-start(x)-width(x)*up
		region.end<-end(x)+width(x)*down
		round(seq(region.start,region.end-win+1+small.num,step))->start
		round(seq(region.start+win-1,region.end+small.num,step))->end
		GRanges(seqnames=seqnames(x),range=IRanges(start=start,end=end))
	}

rowSums(getCoverage(BS.BM,type="Cov"))->Cov
rowSums(getCoverage(BS.BM,type="M"))->M

Cov[Cov > 0]->cg.x.cov
M[Cov>0]->cg.x.m
M[M > min.mc.depth]->mcg.x

cg1<- BS.BM@rowRanges[Cov>0]
mcg1<- BS.BM@rowRanges[M>min.mc.depth]

if(dist=="percent"){
win<-f(region[1],win.percent,step.percent,up.percent,down.percent)
for(i in 2:length(region)){
	tmp.win<-f(region[i],win.percent,step.percent,up.percent,down.percent)
	c(win,tmp.win)->win
}
}
else if(dist=="base"){
win<-generate.sliding.window(region,win.percent,step.base,up.base,down.base)
}
print("starting")

as.data.frame(findOverlaps(win,cg1))->tmp1
countOverlaps(win,mcg1)->mcg
countOverlaps(win,cg1)->c

which(countOverlaps(win,cg1)==0)->tmp2 #regions without Cs

data.table(id=as.character(tmp1$queryHits),cov=cg.x.cov[tmp1$subjectHits],m=cg.x.m[tmp1$subjectHits])->new1
setkey(new1,id)
as.integer(new1[,sum(cov),by=id]$id)->id1
data.table(id=new1[,sum(cov),by=id]$id,start=start(win)[id1],end=end(win)[id1],m=new1[,sum(m),by=id]$V1,cov=new1[,sum(cov),by=id]$V1)->tmp1
tmp1<-rbind(data.frame(tmp1),data.frame(id=as.character(tmp2),start=start(win)[tmp2],end=end(win)[tmp2],m=rep(0,length(tmp2)),cov=rep(0,length(tmp2))))
sort(as.integer(tmp1$id),index.return = T)->ind1
tmp1[ind1$ix,]->tmp1

tmp1$mcg=mcg
tmp1$c=c
tmp1$cov<min.num ->ind
tmp1$cov[ind]<-0
tmp1$m[ind]<-0
tmp1$meth=tmp1$m/tmp1$cov
matrix(tmp1$meth,nrow=length(region))->tbls

if(dist=="percent"){
	pos1<-1
	pos2<-up.percent/step.percent+1
	pos3<-ncol(tbls)-down.percent/step.percent
	pos4<-ncol(tbls)
axis=c(-up.percent,"0 %","100 %",down.percent)
}
else if(dist=="base"){
	pos1<-1
	pos2<-(ncol(tbls)-1/step.percent)/2+1
	pos3<-ncol(tbls)-(ncol(tbls)-1/step.percent)/2
	pos4<-ncol(tbls)

axis=c(-up.base,"0 %","100 %",down.base)
}

if (append==FALSE){
	if(frame==TRUE){
	plot(colMeans(tbls,na.rm=T),ylim=ylim,ylab="methylation level",xlab="",xaxt="n",type="l",main=main,col=col)
	}else{
	plot(colMeans(tbls,na.rm=T),ylim=ylim,ylab="methylation level",xlab="",xaxt="n",yaxt="n",type="l",main=main,col=col,axes=F)
	}
	axis(2)
}else{
	points(colMeans(tbls,na.rm=T),type="l",col=col)

}
if(x.axis==TRUE){axis(1,c(pos1,pos2,pos3,pos4),axis)}
list(tbls=tbls,x.axis.pos=c(pos1,pos2,pos3,pos4),x.axis=axis)->plot.out
return(plot.out)
}

####################### only calculate meth.sliding according to region range
calculate.meth.sliding<-function(x,region,dist="base",win.base=5000,step.base=5000,up.base=200000,down.base=200000,win.percent=0.01,step.percent=0.01,up.percent=1,down.percent=1,min.num=20,min.mc.depth=0,normalize=FALSE,type="average_meth_by_CpG"){
library(data.table)
library(GenomicRanges)
x[x[,5] > 0,]->cg.x
x[x[,4] > min.mc.depth,]->mcg.x

cg1<-GRanges(seqnames=cg.x[,1],IRanges(start=cg.x[,2],end=cg.x[,2]))
mcg1<-GRanges(seqnames=mcg.x[,1],IRanges(start=mcg.x[,2],end=mcg.x[,2]))

if(dist=="base"){
win<-generate.sliding.window(region,win.percent,step.base,up.base,down.base)
}

if(win.percent==1){
win<-generate.sliding.window.center(region,step.base,up.base,from.start=TRUE)
}
print("starting")
as.data.frame(findOverlaps(win,cg1))->tmp1
countOverlaps(win,mcg1)->mcg
countOverlaps(win,cg1)->c
which(countOverlaps(win,cg1)==0)->tmp2 #regions without Cs
data.table(id=as.character(tmp1$queryHits),cov=cg.x[tmp1$subjectHits,5],m=cg.x[tmp1$subjectHits,4])->new1
setkey(new1,id)
as.integer(new1[,sum(cov),by=id]$id)->id1
data.table(id=new1[,sum(cov),by=id]$id,start=start(win)[id1],end=end(win)[id1],m=new1[,sum(m),by=id]$V1,cov=new1[,sum(cov),by=id]$V1)->tmp1
tmp1<-rbind(data.frame(tmp1),data.frame(id=as.character(tmp2),start=start(win)[tmp2],end=end(win)[tmp2],m=rep(0,length(tmp2)),cov=rep(0,length(tmp2))))
sort(as.integer(tmp1$id),index.return = T)->ind1
tmp1[ind1$ix,]->tmp1

tmp1$mcg=mcg
tmp1$c=c
tmp1$cov<min.num ->ind
tmp1$cov[ind]<-0
tmp1$m[ind]<-0
tmp1$meth=tmp1$m/tmp1$cov
#tmp1$meth[ind]<-0

if(normalize==TRUE){
if(type=="average_meth_by_base")
{tmp1$meth*tmp1$c/width(win)->tmp}
else if(type=="total"){
tmp1$meth*tmp1$c->tmp
}
else if(type=="average_meth_by_CpG"){
tmp1$meth->tmp
}
matrix(tmp,nrow=length(region))->tbls
return(tbls)
}else{
matrix(tmp1$meth,nrow=length(region))->tbls
return(tbls)
}
}

####################### only calculate meth.sliding according to region range and smoothed value

calculate.smooth.meth.sliding<-function(meth,cpg.all,region,cpg.with.value=NULL,dist="base",win.base=5000,step.base=5000,up.base=200000,down.base=200000,win.percent=0.01,step.percent=0.01,up.percent=1,down.percent=1,normalize=T,type="average_meth_by_CpG"){
library(data.table)
library(GenomicRanges)
if(!is.null(cpg.with.value)){
cg1<-cpg.with.value}
else{
cg1<-cpg.all
}

stopifnot(length(meth)==length(cg1))

all.cpg<-cpg.all
if(dist=="base"){
win<-generate.sliding.window(region,win.percent,step.base,up.base,down.base)
}

if(win.percent==1){
win<-generate.sliding.window.center(region,step.base,up.base,from.start=TRUE)
}
print("starting")
as.data.frame(findOverlaps(win,cg1))->tmp1
countOverlaps(win,all.cpg)->c
countOverlaps(win,cg1)->cc
tmp2<-c(which(countOverlaps(win,all.cpg)==0),which(countOverlaps(win,all.cpg)>0 & countOverlaps(win,cg1)==0)) #regions without Cs


data.table(id=as.character(tmp1$queryHits),meth=meth[tmp1$subjectHits])->new1
setkey(new1,id)
as.integer(new1[,mean(meth,na.rm=T),by=id]$id)->id1
data.table(id=new1[,mean(meth,na.rm=T),by=id]$id,start=start(win)[id1],end=end(win)[id1],meth=new1[,mean(meth,na.rm=T),by=id]$V1)->tmp1
tmp1<-rbind(data.frame(tmp1),data.frame(id=as.character(tmp2),start=start(win)[tmp2],end=end(win)[tmp2],meth=rep(NA,length(tmp2))))
sort(as.integer(tmp1$id),index.return = T)->ind1
tmp1[ind1$ix,]->tmp1

tmp1$c=c
tmp1$meth<-as.numeric(tmp1$meth)
if(normalize==TRUE){
	if(type=="average_meth_by_base")
		{tmp1$meth*cc/width(win)->tmp}
	else if(type=="total"){
		tmp1$meth*cc->tmp
	}
	else if(type=="average_meth_by_CpG"){
		tmp1$meth*cc/tmp1$c->tmp
	}
	else if(type=="density"){
		cc/tmp1$c->tmp
	}
	else if(type=="None"){
		tmp1$meth->tmp
	}

matrix(tmp,nrow=length(region))->tbls
return(tbls)
}else{
matrix(tmp1$meth,nrow=length(region))->tbls
return(tbls)
}
}
######################## calculate meth level for specific regions
########################
meth.count.region<-function(x,range,type="cg",min.num=20,min.mc.depth=0){
require(data.table)
require(GenomicRanges)

win<-range
if(is.null(names(win))){names(win)<-1:length(win)}
cg2<-GRanges(seqnames=x[,1],IRanges(start=x[,2],end=x[,2]))

if(type=="cg"){
x[,5] > 0 & (x[,3]==1|x[,3]==4)->ind
x[ind,]->x
x[x[,4] > min.mc.depth,]->mcg
}
else{
x[,5] > 0 & !(x[,3]==1|x[,3]==4)->ind
x[ind,]->x
x[x[,4] > min.mc.depth,]->mcg
}



#transfer x, y to GRange objects
cg1<-GRanges(seqnames=x[,1],IRanges(start=x[,2],end=x[,2]))
mcg<-GRanges(seqnames=mcg[,1],IRanges(start=mcg[,2],end=mcg[,2]))

as.data.frame(findOverlaps(win,cg1))->tmp1
#countOverlaps(win,mcg)->mcg

which(countOverlaps(win,cg1)==0)->tmp2 #regions without sequenced Cs

data.table(id=names(win)[tmp1$queryHits],cov=x[,5][tmp1$subjectHits],m=x[,4][tmp1$subjectHits])->new1
setkey(new1,id)
new1[,sum(cov),by=id]$id->id1
data.table(id=new1[,sum(cov),by=id]$id,chr=as.vector(seqnames(win[id1])),start=start(win[id1]),end=end(win[id1]),m=new1[,sum(m),by=id]$V1,cov=new1[,sum(cov),by=id]$V1)->tmp1
tmp1<-rbind(data.frame(tmp1),data.frame(id=names(win)[tmp2],chr=as.vector(seqnames(win[tmp2])),start=start(win)[tmp2],end=end(win)[tmp2],m=rep(0,length(tmp2)),cov=rep(0,length(tmp2))))
order(tmp1$chr,tmp1$start)->ind1
tmp1[ind1,]->tmp1

tmp1$cov<min.num ->ind

tmp1$m_level[ind]<-NA
tmp1$m_level[!ind]<-tmp1$m[!ind]/tmp1$cov[!ind]

tmp1.range=GRanges(seqnames=tmp1$chr,range=IRanges(start=tmp1$start,end=tmp1$end))
countOverlaps(tmp1.range,cg2)->total.c
countOverlaps(tmp1.range,cg1)->covered.c

tmp1$total.c<-total.c
tmp1$covered.c<-covered.c
return(tmp1)
}


#################输入文件为bsseq object
meth.count.region.for.bsseq<-function(BS,range,min.cov=5,min.cpg.num=1){
require(data.table)
require(GenomicRanges)
require(bsseq)
win<-range
cpg.range<-BS@rowRanges

f.tmp<-function(BS.tmp,win,min.cov,min.cpg.num){
if(is.null(names(win))){names(win)<-1:length(win)}
cov<-getCoverage(BS.tmp,type="Cov")
m<-getCoverage(BS.tmp,type="M")
as.data.frame(findOverlaps(win,cpg.range))->tmp1
which(countOverlaps(win,cpg.range)==0)->tmp2 #regions without sequenced Cs
data.table(id=names(win)[tmp1$queryHits],cov=cov[tmp1$subjectHits],m=m[tmp1$subjectHits])->new1
setkey(new1,id)
new1[,sum(cov),by=id]$id->id1
data.table(id=new1[,sum(cov),by=id]$id,chr=as.vector(seqnames(win[id1])),start=start(win[id1]),end=end(win[id1]),m=new1[,sum(m),by=id]$V1,cov=new1[,sum(cov),by=id]$V1)->tmp1
tmp1<-rbind(data.frame(tmp1),data.frame(id=names(win)[tmp2],chr=as.vector(seqnames(win[tmp2])),start=start(win)[tmp2],end=end(win)[tmp2],m=rep(0,length(tmp2)),cov=rep(0,length(tmp2))))
order(tmp1$chr,tmp1$start)->ind1
tmp1[ind1,]->tmp1
tmp1.range=GRanges(seqnames=tmp1$chr,range=IRanges(start=tmp1$start,end=tmp1$end))
countOverlaps(tmp1.range,cpg.range)->total.c
countOverlaps(tmp1.range,cpg.range)->covered.c
tmp1$total.c<-total.c
tmp1$covered.c<-covered.c

tmp1$cov<min.cov | tmp1$covered.c<min.cpg.num ->ind
tmp1$m_level[ind]<-NA
tmp1$m_level[!ind]<-tmp1$m[!ind]/tmp1$cov[!ind]
return(tmp1)
}

out<-f.tmp(BS[,1],win,min.cov,min.cpg.num)[,c(2,3,4,9)]
print (1)
if(length(sampleNames(BS))>1){    
for(i in 2:length(sampleNames(BS))){
print(i)
BS.tmp<-BS[,i]
f.tmp(BS.tmp,win,min.cov,min.cpg.num)->out.tmp
cbind(out,out.tmp[,9])->out
}
}
return(out)

}


####################输入文件为每个cpg位点的meth level，例如450k数据或者用getmeth先计算好meth level
meth.count.region.array<-function(x,range,min.sites=1,data.type="BED"){
print("sort region first, warnning.............................")
library(data.table)
library(GenomicRanges)

win<-range
#transfer x, y to GRange objects

if(data.type=="BED"){
cg1<-GRanges(seqnames=x[,1],IRanges(start=x[,2],end=x[,2]))
}
else{
cg1<-x
}

which(countOverlaps(win,cg1)==0) ->tmp2 #regions without 450k cpg
as.data.frame(findOverlaps(win,cg1))->tmp1

data.table(id=names(win)[tmp1$queryHits],meth=x[,4][tmp1$subjectHits])->new1
setkey(new1,id)
new1[,mean(meth,na.rm=T),by=id]$id->id1
tmp1<-data.table(id=new1[,mean(meth,na.rm=T),by=id]$id,chr=as.vector(seqnames(win[id1])),start=start(win[id1]),end=end(win[id1]),meth=new1[,mean(meth,na.rm=T),by=id]$V1)
tmp1<-rbind(data.frame(tmp1),data.frame(id=names(win)[tmp2],chr=as.vector(seqnames(win[tmp2])),start=start(win)[tmp2],end=end(win)[tmp2],meth=rep("NA",length(tmp2))))
data.frame(tmp1)->tmp1
order(tmp1$chr,tmp1$start)->ind1
tmp1[ind1,]->tmp1
tmp1
}
