sliding.count<-function(file,window.size=10000,step.size=10000,min.num=20,min.mc.depth=0){
library(data.table)
library(GenomicRanges)

read.table(file,colClasses = "integer")->x
x[x$V4 > 0,]->x
x[x$V3 > min.mc.depth,]->mcg
chr.name=tail(unlist(strsplit(file,'/')),n=1)
chr.name=unlist(strsplit(chr.name,'\\.'))[1]
#generate sliding window pos
seq(1,max(x$V1)-window.size+1,step.size)->start
seq(1+window.size-1,max(x$V1),step.size)->end
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
tmp1[tmp1$cov>min.num,]->tmp1
tmp1
}


hypo.region<-function(file1,file2,diff=0.2,merge.dis=20000){
chr.name=tail(unlist(strsplit(file1,'/')),n=1)
chr.name=unlist(strsplit(chr.name,'\\.'))[1]
sliding.count(file1)->s1
sliding.count(file2)->s2
library(GenomicRanges)
intersect(s1$id,s2$id)->ind
s1[s1$id %in% ind,]->s1
s2[s2$id %in% ind,]->s2
s1$m/s1$cov-s2$m/s2$cov >diff ->ind
GRanges(seqnames=chr.name,ranges=IRanges(start=s1[ind,]$start,end=s1[ind,]$end))->gr
gaps(gr)->gap
gap[width(gap)<merge.dis]->gap
reduce(c(gr,gap))->gr
library(rtracklayer)
outname<-paste(chr.name,".hypo.bed",sep="")
export.bed(gr,outname)
}
