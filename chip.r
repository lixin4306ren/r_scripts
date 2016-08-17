#####################
annotatePeaks<-function(peaks,genedb="~/amber3/no_back_up/data/R_data/human.gencodeV19",peak.type="BED",by_gene=TRUE,min.distance=2000,pro.len=500,pro=TRUE,subchr=FALSE,head=FALSE){
library(GenomicFeatures)
library(bumphunter)
library(hash)

if(peak.type=="BED"){
	if(head==TRUE){
		read.table(peaks,header=T)->peak
	}
	else{
		read.table(peaks)->peak
	}
	if(subchr==TRUE){
		GRanges(seqname=sub("chr","",peak[,1]),range=IRanges(start=peak[,2],end=peak[,3]))->peak.range
	}
	else{
		GRanges(seqname=peak[,1],range=IRanges(start=peak[,2],end=peak[,3]))->peak.range
	}
}
else if(peak.type=="RANGE"){
peaks->peak.range
}


names(peak.range)<-1:length(peak.range)

loadDb(genedb)->db
if(by_gene==TRUE){
genes(db)->gene
}
else{
transcriptsBy(db,by="gene")->gene
unlist(gene)->gene
}

annotateNearest(peak.range,gene)->anno
anno$peak_id<-as.numeric(names(peak.range))
abs(anno$dist)<=min.distance -> ind
anno[ind,]->anno
names(gene)[anno$subjectHits]->anno$gene

if(pro==TRUE){
if(by_gene==TRUE){
promoters(gene, upstream=pro.len, downstream=pro.len)->pro
}
else{
promoters(db, upstream=pro.len, downstream=pro.len)->pro
}
if(by_gene==TRUE){
hash(keys=mcols(gene)$gene_id,values=names(gene))->table
values(table,mcols(pro)$gene_id)->names(pro)
}
else{
hash(keys=mcols(gene)$tx_name,values=names(gene))->table
values(table,mcols(pro)$tx_name)->names(pro)
}

annotateNearest(peak.range,pro)->anno.pro
anno.pro$peak_id<-as.numeric(names(peak.range))
abs(anno.pro$dist)==0 & !is.na(anno.pro$dist)-> ind
anno.pro[ind,]->anno.pro
names(pro)[anno.pro$subjectHits]->anno.pro$gene
list(all=anno,promoter=anno.pro)->out
}
else{
list(all=anno)->out
}
return(out)
}


#####################################
annotatePeaksfromgff<-function(peaks,genedb,min.distance=2000,pro=TRUE,subchr=TRUE,header=FALSE){
if(header==TRUE){read.table(peaks,header=T)->peak}else{read.table(peaks)->peak}
library(GenomicFeatures)
library(bumphunter)
if(subchr==TRUE){GRanges(seqname=sub("chr","",peak[,1]),range=IRanges(start=peak[,2],end=peak[,3]))->peak.range}else{GRanges(seqname=peak[,1],range=IRanges(start=peak[,2],end=peak[,3]))->peak.range}
read.table(genedb,sep="\t")->gene.gff

f2<-function(x){sub("Name=","",unlist(strsplit(as.vector(x),";"))[2])}
f<-function(x){sub("ID=","",unlist(strsplit(as.vector(x),";"))[1])}
sapply(as.vector(gene.gff[,9]),f)->geneid
sapply(as.vector(gene.gff[,9]),f2)->genename

GRanges(seqnames=gene.gff[,1],range=IRanges(start=gene.gff[,4],end=gene.gff[,5]),strand=gene.gff[,7])->gene
names(gene)<-geneid

annotateNearest(peak.range,gene)->anno
abs(anno$dist)<=min.distance -> ind
anno[ind,]->anno
names(gene)[anno$matchIndex]->anno$geneid
names(genename)[anno$matchIndex]->anno$genename

if(pro==TRUE){
flank(gene,min.distance,both=T)->pro
annotateNearest(peak.range,pro)->anno.pro
abs(anno.pro$dist)==0 -> ind
anno.pro[ind,]->anno.pro
names(pro)[anno.pro$matchIndex]->anno.pro$geneid
names(genename)[anno.pro$matchIndex]->anno.pro$genename

list(all=anno,promoter=anno.pro)->out
}
else{
list(all=anno)->out
}
return(out)
}

#####################################

annotatePeaks2<-function(peaks,genedb="~xinli/amber3/no_back_up/data/R_data/human.GRCh37.p8.gene.annotation",min.distance=2000,pro=T,subchr=TRUE){
peaks->peak
peak[,1]!="chrM" ->ind
peak[ind,]->peak
library(GenomicFeatures)
library(bumphunter)
library(hash)
if(subchr==TRUE){
GRanges(seqname=sub("chr","",peak[,1]),range=IRanges(start=peak[,2],end=peak[,3]))->peak.range
}
else{
GRanges(seqname=peak[,1],range=IRanges(start=peak[,2],end=peak[,3]))->peak.range
}
loadDb(genedb)->db
transcriptsBy(db,by="gene")->gene
unlist(gene)->gene

annotateNearest(peak.range,gene)->anno
#peak$pval->anno$pval
#peak$padj->anno$padj
names(gene)[anno$matchIndex]->anno$gene
cbind(peak,anno)->anno
abs(anno$dist)< min.distance -> ind
anno[ind,]->anno

if(pro==T){
promoters(db, upstream=2000, downstream=2000)->pro
hash(keys=mcols(gene)$tx_name,values=names(gene))->table
values(table,mcols(pro)$tx_name)->names(pro)
annotateNearest(peak.range,pro)->anno.pro
#peak$pval->anno.pro$pval
#peak$padj->anno.pro$padj
names(pro)[anno.pro$matchIndex]->anno.pro$gene
#peak$nearestGene<-anno.pro$gene
cbind(peak,anno.pro)->anno.pro

abs(anno.pro$dist)==0 -> ind
anno.pro[ind,]->anno.pro
#peak[ind,]->peak.out
list(all=anno,promoter=anno.pro)->out
}
else{
list(all=anno)->out
}
return(out)
}

#####################################

get_common_peak<-function(peak1,peak2,return.range=FALSE,min.overlap=0){
library(GenomicRanges)
read.table(peak1)->peak1
read.table(peak2)->peak2
GRanges(seqname=peak1$V1,range=IRanges(start=peak1$V2,end=peak1$V3))->peak1.range
GRanges(seqname=peak2$V1,range=IRanges(start=peak2$V2,end=peak2$V3))->peak2.range

peak1.range %over% peak2.range ->ind1
peak2.range %over% peak1.range ->ind2
peak1.range[ind1]->peak1.range
peak2.range[ind2]->peak2.range

c(peak1.range,peak2.range)-> both.peak.range
reduce(both.peak.range)-> both.peak.range
if(return.range==FALSE){
data.frame(chr=as.vector(seqnames(both.peak.range)),start=start(both.peak.range),end=end(both.peak.range))->both.peak
return(both.peak)
}
else{
return(both.peak.range)
}
}

#####################################

get_common_peak3<-function(peak1,peak2,return.range=FALSE,min.overlap=0.5){
library(GenomicRanges)
bedtools.coveragebed.range(bed1=peak2,bed2=peak1,mcol=F)->cov1
bedtools.coveragebed.range(bed1=peak1,bed2=peak2,mcol=F)->cov2

cov1[,ncol(cov1)]>min.overlap ->ind1
cov2[,ncol(cov2)]>min.overlap ->ind2
cov1[ind1,c(1,2,3)]->peak1
cov2[ind2,c(1,2,3)]->peak2

GRanges(seqname=peak1$V1,range=IRanges(start=peak1$V2,end=peak1$V3))->peak1.range
GRanges(seqname=peak2$V1,range=IRanges(start=peak2$V2,end=peak2$V3))->peak2.range

peak1.range %over% peak2.range ->ind1
peak2.range %over% peak1.range ->ind2
peak1.range[ind1]->peak1.range
peak2.range[ind2]->peak2.range

c(peak1.range,peak2.range)-> both.peak.range
reduce(both.peak.range)-> both.peak.range
if(return.range==FALSE){
data.frame(chr=as.vector(seqnames(both.peak.range)),start=start(both.peak.range),end=end(both.peak.range))->both.peak
return(both.peak)
}
else{
return(both.peak.range)
}
}

##################################### inputs are grange objects

get_common_peak2<-function(peak1,peak2,return.range=FALSE){
library(GenomicRanges)
peak1->peak1.range
peak2->peak2.range

peak1.range %over% peak2.range ->ind1
peak2.range %over% peak1.range ->ind2
peak1.range[ind1]->peak1.range
peak2.range[ind2]->peak2.range


c(peak1.range,peak2.range)-> both.peak.range
reduce(both.peak.range)-> both.peak.range
if(return.range==FALSE){
data.frame(chr=as.vector(seqnames(both.peak.range)),start=start(both.peak.range),end=end(both.peak.range))->both.peak
return(both.peak)
}
else{
return(both.peak.range)
}
}

#####################################

get_all_peak<-function(peak1,peak2,return.range=FALSE,min.overlap=0){
library(GenomicRanges)
read.table(peak1)->peak1
read.table(peak2)->peak2
GRanges(seqname=peak1$V1,range=IRanges(start=peak1$V2,end=peak1$V3))->peak1.range
GRanges(seqname=peak2$V1,range=IRanges(start=peak2$V2,end=peak2$V3))->peak2.range

reduce(c(peak1.range,peak2.range))-> all.peak.range
if(return.range==FALSE){
data.frame(chr=as.vector(seqnames(all.peak.range)),start=start(all.peak.range),end=end(all.peak.range))->all.peak
return(all.peak)
}
else{
return(all.peak.range)
}
}

#####################################
get_unique_peak<-function(peak1,peak2){
library(GenomicRanges)
read.table(peak1)->peak1
read.table(peak2)->peak2

GRanges(seqname=peak1$V1,range=IRanges(start=peak1$V2,end=peak1$V3))->peak1.range
GRanges(seqname=peak2$V1,range=IRanges(start=peak2$V2,end=peak2$V3))->peak2.range

peak1.range %over% peak2.range ->ind1
peak1.range[!ind1]->peak1.range

data.frame(chr=as.vector(seqnames(peak1.range)),start=start(peak1.range),end=end(peak1.range))->unique.peak
return(unique.peak)
}

#####################################

peak.pattern.dist<-function(block1,block2){
library(GenomicRanges)
read.table(block1)->block1
read.table(block2)->block2
GRanges(seqnames=block1[,1],ranges=IRanges(start=block1[,2],end=block1[,3]))->block1
GRanges(seqnames=block2[,1],ranges=IRanges(start=block2[,2],end=block2[,3]))->block2

dist<-sum(width(intersect(block1,block2)))/sum(width(union(block1,block2)))
dist
}

#####################################

merge_peaks<-function(peak=peakfile,min.dis=10000,remove.small=NULL,remove.small0=0,return.range=FALSE){
######### first remove region less than remove.small0
library(GenomicFeatures)
read.table(peak)->x
GRanges(seqnames=x[,1],range=IRanges(start=x[,2],end=x[,3]))->y
width(y)>=remove.small0 ->ind
y[ind]->y
gaps(y)->gap
width(gap)< min.dis ->ind
reduce(c(y,gap[ind]))->z

######### after merge by mis.dis, remove merged region less than remove.small
if(!is.null(remove.small)){
width(z)>=remove.small -> ind
z[ind]->z
}
if(return.range==FALSE){
data.frame(chr=as.vector(seqnames(z)),start=start(z),end=end(z))->out
return(out)
}
else{
return(z)
}

}

#####################################

merge_peaks2<-function(peak=peakobject,min.dis=10000,remove.small0=0,remove.small=NULL,return.range=FALSE){
library(GenomicFeatures)
peak->x
GRanges(seqnames=x[,1],range=IRanges(start=x[,2],end=x[,3]))->y
width(y)>=remove.small0 ->ind
y[ind]->y
gaps(y)->gap
width(gap)< min.dis ->ind
reduce(c(y,gap[ind]))->z
if(!is.null(remove.small)){
width(z)>=remove.small -> ind
z[ind]->z
}
if(return.range==FALSE){
data.frame(chr=as.vector(seqnames(z)),start=start(z),end=end(z))->out
return(out)
}
else{
return(z)
}
}

#                                                                              min.dis=10000,remove.small0=0,remove.small=NULL,remove.small.for.diff=NULLremove.small2=0,min.dis.for.diff=10000, remove.small0.for.diff=0,remove.small.each.replicate=NULL,name1=NULL,name2=NULL,name3=NULL){
#####################################

diff_peaks<-function(peak1=NULL,peak2=NULL,remove.peak1=NULL,remove.peak2=NULL,min.dis=10000,remove.small0=0,remove.small=NULL,remove.small.for.diff=NULL,remove.small2=0,min.dis.for.diff=10000,remove.small0.for.diff=0,name1=NULL,name2=NULL,name3=NULL){
if(!is.null(remove.small)){merge_peaks(peak1,min.dis,remove.small=remove.small,return.range=T,remove.small0=remove.small0)->merge1}
else{merge_peaks(peak1,min.dis,return.range=T,remove.small0=remove.small0)->merge1}
if(!is.null(remove.small)){merge_peaks(peak2,min.dis,remove.small=remove.small,return.range=T,remove.small0=remove.small0)->merge2}
else{merge_peaks(peak2,min.dis,return.range=T,remove.small0=remove.small0)->merge2}


if(!is.null(remove.peak1) | !is.null(remove.peak2)){
read.table(remove.peak1)->tmp
GRanges(seqnames=tmp[,1],range=IRanges(start=tmp[,2],end=tmp[,3]))->tmp.range
setdiff(merge1,tmp.range)->merge1
}
data.frame(chr=as.vector(seqnames(merge1)),start=start(merge1),end=end(merge1))->tmp
file.name<-paste(name1,name3,"bed",sep=".")
write.table(tmp,file=file.name,quote=F,col.names=F,row.names=F,sep="\t")

if(!is.null(remove.peak1) | !is.null(remove.peak2)){
read.table(remove.peak2)->tmp
GRanges(seqnames=tmp[,1],range=IRanges(start=tmp[,2],end=tmp[,3]))->tmp.range
setdiff(merge2,tmp.range)->merge2
}

data.frame(chr=as.vector(seqnames(merge2)),start=start(merge2),end=end(merge2))->tmp
file.name<-paste(name2,name3,"bed",sep=".")
write.table(tmp,file=file.name,quote=F,col.names=F,row.names=F,sep="\t")

merge_peaks(peak1,min.dis.for.diff,return.range=T,remove.small0=remove.small0.for.diff)->merge1.for.diff
merge_peaks(peak2,min.dis.for.diff,return.range=T,remove.small0=remove.small0.for.diff)->merge2.for.diff

######output all domain
data.frame(chr=as.vector(seqnames(merge1.for.diff)),start=start(merge1.for.diff),end=end(merge1.for.diff))->tmp
file.name<-paste(name1,name3,"all.bed",sep=".")
write.table(tmp,file=file.name,quote=F,col.names=F,row.names=F,sep="\t")
data.frame(chr=as.vector(seqnames(merge2.for.diff)),start=start(merge2.for.diff),end=end(merge2.for.diff))->tmp
file.name<-paste(name2,name3,"all.bed",sep=".")
write.table(tmp,file=file.name,quote=F,col.names=F,row.names=F,sep="\t")

setdiff(merge1,merge2.for.diff)->peak.up
setdiff(merge2,merge1.for.diff)->peak.down
width(peak.up)>=remove.small2 ->ind1
width(peak.down)>=remove.small2 ->ind2
peak.up[ind1]->peak.up
peak.down[ind2]->peak.down
data.frame(chr=as.vector(seqnames(peak.up)),start=start(peak.up),end=end(peak.up),direction="Up")->out1
data.frame(chr=as.vector(seqnames(peak.down)),start=start(peak.down),end=end(peak.down),direction="Down")->out2
rbind(out1,out2)->out
file.name<-paste(name1,name2,name3,"diff.bed",sep=".")
write.table(out,file=file.name,quote=F,col.names=F,row.names=F,sep="\t")

}
#####################################

diff_peaks_replicate<-function(peak1=NULL,peak2=NULL,peak3=NULL,peak4=NULL,min.dis=10000,min.dis.for.diff=10000,remove.small0=0,remove.small0.for.diff=0,remove.small.for.diff=NULL,remove.small=NULL,remove.small.each.replicate=NULL,remove.small2=0,name1=NULL,name2=NULL,name3=NULL){
if(!is.null(remove.small)){merge_peaks(peak1,min.dis,remove.small=remove.small.each.replicate,return.range=T,remove.small0=remove.small0)->merge1}
else{merge_peaks(peak1,min.dis,return.range=T,remove.small0=remove.small0)->merge1}
if(!is.null(remove.small)){merge_peaks(peak2,min.dis,remove.small=remove.small.each.replicate,return.range=T,remove.small0=remove.small0)->merge2}
else{merge_peaks(peak2,min.dis,return.range=T,remove.small0=remove.small0)->merge2}
if(!is.null(remove.small)){merge_peaks(peak3,min.dis,remove.small=remove.small.each.replicate,return.range=T,remove.small0=remove.small0)->merge3}
else{merge_peaks(peak3,min.dis,return.range=T,remove.small0=remove.small0)->merge3}
if(!is.null(remove.small)){merge_peaks(peak4,min.dis,remove.small=remove.small.each.replicate,return.range=T,remove.small0=remove.small0)->merge4}
else{merge_peaks(peak4,min.dis,return.range=T,remove.small0=remove.small0)->merge4}

get_common_peak3(merge1,merge2,return.range=T)->merge.peak1
width(merge.peak1)>=remove.small ->ind1
merge.peak1[ind1]->merge.peak1
data.frame(chr=as.vector(seqnames(merge.peak1)),start=start(merge.peak1),end=end(merge.peak1))->tmp
file.name<-paste(name1,name3,"large.bed",sep=".")
write.table(tmp,file=file.name,quote=F,col.names=F,row.names=F,sep="\t")
get_common_peak3(merge3,merge4,return.range=T)->merge.peak2
width(merge.peak2)>=remove.small ->ind2
merge.peak2[ind2]->merge.peak2
data.frame(chr=as.vector(seqnames(merge.peak2)),start=start(merge.peak2),end=end(merge.peak2))->tmp
file.name<-paste(name2,name3,"large.bed",sep=".")
write.table(tmp,file=file.name,quote=F,col.names=F,row.names=F,sep="\t")


merge_peaks(peak1,min.dis=min.dis.for.diff,return.range=T,remove.small0=remove.small0.for.diff)->merge1.for.diff
merge_peaks(peak2,min.dis=min.dis.for.diff,return.range=T,remove.small0=remove.small0.for.diff)->merge2.for.diff
merge_peaks(peak3,min.dis=min.dis.for.diff,return.range=T,remove.small0=remove.small0.for.diff)->merge3.for.diff
merge_peaks(peak4,min.dis=min.dis.for.diff,return.range=T,remove.small0=remove.small0.for.diff)->merge4.for.diff

get_common_peak3(merge1.for.diff,merge2.for.diff,return.range=T)->merge.peak1.for.diff
width(merge.peak1.for.diff)>=remove.small.for.diff ->ind1
merge.peak1.for.diff[ind1]->merge.peak1.for.diff

get_common_peak3(merge3.for.diff,merge4.for.diff,return.range=T)->merge.peak2.for.diff
width(merge.peak2.for.diff)>=remove.small.for.diff ->ind2
merge.peak2.for.diff[ind2]->merge.peak2.for.diff


######output all domain
data.frame(chr=as.vector(seqnames(merge.peak1.for.diff)),start=start(merge.peak1.for.diff),end=end(merge.peak1.for.diff))->tmp
file.name<-paste(name1,name3,"all.bed",sep=".")
write.table(tmp,file=file.name,quote=F,col.names=F,row.names=F,sep="\t")
data.frame(chr=as.vector(seqnames(merge.peak2.for.diff)),start=start(merge.peak2.for.diff),end=end(merge.peak2.for.diff))->tmp
file.name<-paste(name2,name3,"all.bed",sep=".")
write.table(tmp,file=file.name,quote=F,col.names=F,row.names=F,sep="\t")


setdiff(merge.peak1,reduce(c(merge.peak2,merge.peak2.for.diff)))->peak.up
setdiff(merge.peak2,reduce(c(merge.peak1,merge.peak1.for.diff)))->peak.down
#setdiff(merge.peak1,merge.peak2)->peak.up
#setdiff(merge.peak2,merge.peak1)->peak.down

width(peak.up)>=remove.small2 ->ind1
width(peak.down)>=remove.small2 ->ind2
peak.up[ind1]->peak.up
peak.down[ind2]->peak.down
data.frame(chr=as.vector(seqnames(peak.up)),start=start(peak.up),end=end(peak.up),direction="Up")->out1
data.frame(chr=as.vector(seqnames(peak.down)),start=start(peak.down),end=end(peak.down),direction="Down")->out2
rbind(out1,out2)->out
file.name<-paste(name1,name2,name3,"diff.bed",sep=".")
write.table(out,file=file.name,quote=F,col.names=F,row.names=F,sep="\t")

}
#####################################

plot_chip_seq_replicates<-function(tmp.list,region=NULL){
as.vector(read.table(tmp.list)$V1)->list
data.frame(V1=read.table(list[1])$V4)->x
read.table(list[1])->y
GRanges(seqnames=y[,1],range=IRanges(start=y[,2],end=y[,3]))->all.region
if(is.null(region)){
all.region %over% all.region -> ind
}
else{
all.region %over% region -> ind
}

#filename<-paste(name,"pdf",sep=".")
#pdf(filename)
name<-basename(unlist(strsplit(list[1],"\\."))[1])
for(i in 2:(length(list))){
	read.table(list[i])->tmp
	cbind(x,tmp[,4])->x
	tmp.name<-basename(unlist(strsplit(list[i],"\\."))[1])
	append(name,tmp.name)->name
}
x[ind,]->x
na.omit(x)->x
colnames(x)<-name
return(as.dist(1-cor(x)))
#plot(hclust(as.dist(1-cor(x))),main=title)
#dev.off
}

#####################################


plot.chip.comp<-function(no.strand=TRUE,anno.type="BED",anno,anno2=NULL,chip1,chip2,input1,input2,up=200000,down=200000,flanking.freq=5000,up.percent=100,down.percent=100,freq=1,flanking.dist="base",dist="percent",cap.q = 0.99,s.width=5000,frame=TRUE,main=NULL,legend=NULL,paste.chr=FALSE,ylim=c(0,1),gapped=FALSE,append=FALSE,text.col=c("black","red"),max.cvg=NULL,min.cvg=NULL,x.axis=TRUE,y.axis=TRUE,logratio=FALSE){
	library(Repitools)
	library(GenomicRanges)
	if(gapped==FALSE){BAMpair2GRanges(chip1)->x1}else{BAM2GRanges(chip1)->x1}
	if(paste.chr==TRUE){paste("chr",seqlevels(x1),sep="")->seqlevels(x1)}
	if(gapped==FALSE){BAMpair2GRanges(chip2)->x2}else{BAM2GRanges(chip2)->x2}
	if(paste.chr==TRUE){paste("chr",seqlevels(x2),sep="")->seqlevels(x2)}
	if(gapped==FALSE){BAMpair2GRanges(input1)->x3}else{BAM2GRanges(input1)->x3}
	if(paste.chr==TRUE){paste("chr",seqlevels(x3),sep="")->seqlevels(x3)}
	if(gapped==FALSE){BAMpair2GRanges(input2)->x4}else{BAM2GRanges(input2)->x4}
	if(paste.chr==TRUE){paste("chr",seqlevels(x4),sep="")->seqlevels(x4)}
	
	if(anno.type=="BED"){
	read.table(anno)->y
		if(no.strand==TRUE){
			GRanges(seqnames=y[,1],range=IRanges(start=y[,2],end=y[,3]),strand="+")->region
		}
		else{
			GRanges(seqnames=y[,1],range=IRanges(start=y[,2],end=y[,3]),strand=y[,6])->region
		}
	}
	else if(anno.type=="RANGE"){
		if(no.strand==TRUE){
			strand(anno)<-"+"
			anno->region
		}
		else{
			anno->region
		}
	}
	
	if(flanking.dist=="base"){
	if(freq!=100){fs1.middle=featureScores(x1,region,up=0,down=100,freq=freq,verbose=TRUE,dist=dist,s.width=s.width)}
	fs1.left=featureScores(x1,region,up=up,down=0,freq=flanking.freq,verbose=TRUE,dist=flanking.dist,s.width=s.width)
	if(freq!=100){fs3.middle=featureScores(x3,region,up=0,down=100,freq=freq,verbose=TRUE,dist=dist,s.width=s.width)}
	fs3.left=featureScores(x3,region,up=up,down=0,freq=flanking.freq,verbose=TRUE,dist=flanking.dist,s.width=s.width)
	if(no.strand==TRUE){
		strand(region)<-"-"
	}
	else if(no.strand==FALSE){
		strand(region)=="+" -> ind1
		strand(region)=="-" -> ind2
		strand(region[ind1])<- "-"
		strand(region[ind2])<- "+"
	}
	fs1.right=featureScores(x1,region,up=down,down=0,freq=flanking.freq,verbose=TRUE,dist=flanking.dist,s.width=s.width)
	fs3.right=featureScores(x3,region,up=down,down=0,freq=flanking.freq,verbose=TRUE,dist=flanking.dist,s.width=s.width)
	}
	else if(flanking.dist=="percent"){
	fs1=featureScores(x1,region,up=up.percent,down=down.percent,freq=freq,verbose=TRUE,dist=dist,s.width=s.width)
	fs3=featureScores(x3,region,up=up.percent,down=down.percent,freq=freq,verbose=TRUE,dist=dist,s.width=s.width)
	}

	if(flanking.dist=="base"){
	if(is.null(anno2)){
		if(freq!=100){fs2.middle=featureScores(x2,region,up=0,down=100,freq=freq,verbose=TRUE,dist=dist,s.width=s.width)}
		fs2.left=featureScores(x2,region,up=up,down=0,freq=flanking.freq,verbose=TRUE,dist=flanking.dist,s.width=s.width)
		if(freq!=100){fs4.middle=featureScores(x4,region,up=0,down=100,freq=freq,verbose=TRUE,dist=dist,s.width=s.width)}
		fs4.left=featureScores(x4,region,up=up,down=0,freq=flanking.freq,verbose=TRUE,dist=flanking.dist,s.width=s.width)
	if(no.strand==TRUE){
		strand(region)<-"-"
	}
	else if(no.strand==FALSE){
		strand(region)=="+" -> ind1
		strand(region)=="-" -> ind2
		strand(region[ind1])<- "-"
		strand(region[ind2])<- "+"
	}
	fs2.right=featureScores(x2,region,up=down,down=0,freq=flanking.freq,verbose=TRUE,dist=flanking.dist,s.width=s.width)
	fs4.right=featureScores(x4,region,up=down,down=0,freq=flanking.freq,verbose=TRUE,dist=flanking.dist,s.width=s.width)
	}
	else{

	if(anno.type=="BED"){
	read.table(anno2)->y2
		if(no.strand==TRUE){
			GRanges(seqnames=y2[,1],range=IRanges(start=y2[,2],end=y2[,3]),strand="+")->region2
		}
		else{
			GRanges(seqnames=y2[,1],range=IRanges(start=y2[,2],end=y2[,3]),strand=y2[,6])->region2
		}
	}
	else if(anno.type=="RANGE"){
		if(no.strand==TRUE){
			strand(anno2)<-"+"
			anno2->region2
		}
		else{
			anno2->region2
		}
	}
		
		if(freq!=100){fs2.middle=featureScores(x2,region2,up=0,down=100,freq=freq,verbose=TRUE,dist=dist,s.width=s.width)}
		fs2.left=featureScores(x2,region2,up=up,down=0,freq=flanking.freq,verbose=TRUE,dist=flanking.dist,s.width=s.width)
		if(freq!=100){fs4.middle=featureScores(x4,region2,up=0,down=100,freq=freq,verbose=TRUE,dist=dist,s.width=s.width)}
		fs4.left=featureScores(x4,region2,up=up,down=0,freq=flanking.freq,verbose=TRUE,dist=flanking.dist,s.width=s.width)
		strand(region2)<-"-"
		fs2.right=featureScores(x2,region2,up=down,down=0,freq=flanking.freq,verbose=TRUE,dist=flanking.dist,s.width=s.width)
		fs4.right=featureScores(x4,region2,up=down,down=0,freq=flanking.freq,verbose=TRUE,dist=flanking.dist,s.width=s.width)
	}
	}
	else if(flanking.dist=="percent"){
	if(is.null(anno2)){
		read.table(anno)->y
		GRanges(seqnames=y[,1],range=IRanges(start=y[,2],end=y[,3]),strand="+")->region
		fs2=featureScores(x2,region,up=up.percent,down=down.percent,freq=freq,verbose=TRUE,dist=dist,s.width=s.width)
		fs4=featureScores(x4,region,up=up.percent,down=down.percent,freq=freq,verbose=TRUE,dist=dist,s.width=s.width)
	}
	else{
		read.table(anno2)->y2
		GRanges(seqnames=y2[,1],range=IRanges(start=y2[,2],end=y2[,3]),strand="+")->region2
		fs2=featureScores(x2,region2,up=up.percent,down=down.percent,freq=freq,verbose=TRUE,dist=dist,s.width=s.width)
		fs4=featureScores(x4,region2,up=up.percent,down=down.percent,freq=freq,verbose=TRUE,dist=dist,s.width=s.width)
	}
	}
	
if(flanking.dist=="base"){	
if(freq!=100){cbind(Repitools::tables(fs1.left)[[1]],Repitools::tables(fs1.middle)[[1]],t(apply(Repitools::tables(fs1.right)[[1]],1,rev)))->tmp1}
else{cbind(Repitools::tables(fs1.left)[[1]],t(apply(Repitools::tables(fs1.right)[[1]],1,rev)))->tmp1}
if(freq!=100){cbind(Repitools::tables(fs2.left)[[1]],Repitools::tables(fs2.middle)[[1]],t(apply(Repitools::tables(fs2.right)[[1]],1,rev)))->tmp2}
else{cbind(Repitools::tables(fs2.left)[[1]],t(apply(Repitools::tables(fs2.right)[[1]],1,rev)))->tmp2}
if(freq!=100){cbind(Repitools::tables(fs3.left)[[1]],Repitools::tables(fs3.middle)[[1]],t(apply(Repitools::tables(fs3.right)[[1]],1,rev)))->tmp3}
else{cbind(Repitools::tables(fs3.left)[[1]],t(apply(Repitools::tables(fs3.right)[[1]],1,rev)))->tmp3}
if(freq!=100){cbind(Repitools::tables(fs4.left)[[1]],Repitools::tables(fs4.middle)[[1]],t(apply(Repitools::tables(fs4.right)[[1]],1,rev)))->tmp4}
else{cbind(Repitools::tables(fs4.left)[[1]],t(apply(Repitools::tables(fs4.right)[[1]],1,rev)))->tmp4}
}
else if(flanking.dist=="percent"){
Repitools::tables(fs1)[[1]]->tmp1
Repitools::tables(fs2)[[1]]->tmp2
Repitools::tables(fs3)[[1]]->tmp3
Repitools::tables(fs4)[[1]]->tmp4
}

if(logratio==FALSE){
tmp1-tmp3->out1
tmp2-tmp4->out2
}
else{
log2((tmp1+1e-10)/(tmp3+1e-10))->out1
log2((tmp2+1e-10)/(tmp4+1e-10))->out2
}

if(is.null(max.cvg)){
max.cvg <- quantile(c(unlist(out1),unlist(out2)), cap.q, na.rm = TRUE)
min.cvg <- quantile(c(unlist(out1),unlist(out2)), 1-cap.q, na.rm = TRUE)
}
else{
max.cvg<-max.cvg
min.cvg<-min.cvg
}
print(max.cvg)
print(min.cvg)
list(max=max.cvg,min=min.cvg)->cvg.cutoff
    f.tmp<-function(x)
    {
        x[x > max.cvg] = max.cvg
        x[x < min.cvg] = min.cvg
        (x-min(x)) / (max(x)-min(x))
    }

    f.tmp2<-function(x){
        x[x > max.cvg] = max.cvg
        x[x < min.cvg] = min.cvg
	x
    }
    
    if(logratio==FALSE){
    f.tmp(out1)->tbls1
    f.tmp(out2)->tbls2
    }
    else{
    f.tmp2(out1)->tbls1
    f.tmp2(out2)->tbls2
    }
	x.lab<-names(colMeans(tbls1,na.rm=T))
	pos1<-1
	which(x.lab=="0 %")->pos2
	which(x.lab=="100 %")->pos3
	pos4<-length(x.lab)

#max(colMeans(tbls1),colMeans(tbls2))
if(append==FALSE){
	if(frame==TRUE){
	plot(colMeans(tbls1,na.rm=T),ylim=ylim,ylab="Normalized read density",xlab="",xaxt="n",yaxt="n",type="l",main=main,col=text.col[1])
	points(colMeans(tbls2,na.rm=T),type="l",col=text.col[2])
	}else{
	plot(colMeans(tbls1,na.rm=T),ylim=ylim,ylab="Normalized read density",xlab="",xaxt="n",yaxt="n",type="l",axes=F,main=main,col=text.col[1])
	points(colMeans(tbls2,na.rm=T),type="l",col=text.col[2])
	}
	if (y.axis==TRUE){axis(2)}

}
else{
	points(colMeans(tbls1,na.rm=T),type="l",col=text.col[1])
	points(colMeans(tbls2,na.rm=T),type="l",col=text.col[2])
}

if (x.axis==TRUE){axis(1,c(pos1,pos2,pos3,pos4),x.lab[c(pos1,pos2,pos3,pos4)])}

	if(!is.null(legend)){
		legend("topright", legend = legend,  text.col = text.col)
	}

list(cvg.cutoff=cvg.cutoff,tbls1=tbls1,tbls2=tbls2,x.axis.pos=c(pos1,pos2,pos3,pos4),x.axis=x.lab[c(pos1,pos2,pos3,pos4)])->plot.out
return(plot.out)
}

plot.chip.comp.around.TSS<-function(anno,anno2=NULL,chip1,chip2,input1,input2,up=10000,down=10000,freq=200,cap.q = 0.99,s.width=100,frame=TRUE,main=NULL,legend=NULL,paste.chr=FALSE,ylim=c(0,1),gapped=FALSE,append=FALSE,text.col=c("black","red"),max.cvg=NULL,min.cvg=NULL,x.axis=TRUE,y.axis=TRUE,logratio=FALSE,input.type="BAM",anno.type="BED"){
	require(Repitools)
	require(GenomicRanges)
	if(input.type=="BAM"){
	if(gapped==FALSE){BAMpair2GRanges(chip1)->x1}else{BAM2GRanges(chip1)->x1}
	if(paste.chr==TRUE){paste("chr",seqlevels(x1),sep="")->seqlevels(x1)}
	if(gapped==FALSE){BAMpair2GRanges(chip2)->x2}else{BAM2GRanges(chip2)->x2}
	if(paste.chr==TRUE){paste("chr",seqlevels(x2),sep="")->seqlevels(x2)}
	if(gapped==FALSE){BAMpair2GRanges(input1)->x3}else{BAM2GRanges(input1)->x3}
	if(paste.chr==TRUE){paste("chr",seqlevels(x3),sep="")->seqlevels(x3)}
	if(gapped==FALSE){BAMpair2GRanges(input2)->x4}else{BAM2GRanges(input2)->x4}
	if(paste.chr==TRUE){paste("chr",seqlevels(x4),sep="")->seqlevels(x4)}
	}
	else if(input.type=="BED"){
	read.table(chip1)->tmp
	GRanges(tmp[,1],range=IRanges(start=tmp[,2],end=tmp[,3]))->x1
	read.table(chip2)->tmp
	GRanges(tmp[,1],range=IRanges(start=tmp[,2],end=tmp[,3]))->x2
	read.table(input1)->tmp
	GRanges(tmp[,1],range=IRanges(start=tmp[,2],end=tmp[,3]))->x3
	read.table(input2)->tmp
	GRanges(tmp[,1],range=IRanges(start=tmp[,2],end=tmp[,3]))->x4
	}
	
		if(anno.type=="BED"){
		read.table(anno)->y
		GRanges(seqnames=y[,1],range=IRanges(start=y[,2],end=y[,3]),strand=y[,6])->region
		}
		else if(anno.type=="RANGE"){
			anno->region
		}

		fs1=featureScores(x1,region,up=up,down=down,freq=freq,verbose=TRUE,s.width=s.width)
		fs3=featureScores(x3,region,up=up,down=down,freq=freq,verbose=TRUE,s.width=s.width)


	if(is.null(anno2)){
		if(anno.type=="BED"){
		read.table(anno)->y
		GRanges(seqnames=y[,1],range=IRanges(start=y[,2],end=y[,3]),strand=y[,6])->region
		}
		else if(anno.type=="RANGE"){
			anno->region
		}
		fs2=featureScores(x2,region,up=up,down=down,freq=freq,verbose=TRUE,s.width=s.width)
		fs4=featureScores(x4,region,up=up,down=down,freq=freq,verbose=TRUE,s.width=s.width)
	}
	else{
		if(anno.type=="BED"){
		read.table(anno2)->y2
		GRanges(seqnames=y2[,1],range=IRanges(start=y2[,2],end=y2[,3]),strand=y2[,6])->region2
		}
		else if(anno.type=="RANGE"){
			anno2->region2
		}
		fs2=featureScores(x2,region2,up=up,down=down,freq=freq,verbose=TRUE,s.width=s.width)
		fs4=featureScores(x4,region2,up=up,down=down,freq=freq,verbose=TRUE,s.width=s.width)
	}

	
tables(fs1)[[1]]->tmp1
tables(fs2)[[1]]->tmp2
tables(fs3)[[1]]->tmp3
tables(fs4)[[1]]->tmp4



if(logratio==FALSE){
tmp1-tmp3->out1
tmp2-tmp4->out2
}
else{
log2((tmp1+1e-10)/(tmp3+1e-10))->out1
log2((tmp2+1e-10)/(tmp4+1e-10))->out2
}

if(is.null(max.cvg)){
max.cvg <- quantile(c(unlist(out1),unlist(out2)), cap.q, na.rm = TRUE)
min.cvg <- quantile(c(unlist(out1),unlist(out2)), 1-cap.q, na.rm = TRUE)
}
else{
max.cvg<-max.cvg
min.cvg<-min.cvg
}
print(max.cvg)
print(min.cvg)
list(max=max.cvg,min=min.cvg)->cvg.cutoff
    f.tmp<-function(x)
    {
        x[x > max.cvg] = max.cvg
        x[x < min.cvg] = min.cvg
        (x-min(x)) / (max(x)-min(x))
    }

    f.tmp2<-function(x){
        x[x > max.cvg] = max.cvg
        x[x < min.cvg] = min.cvg
	x
    }
    
    if(logratio==FALSE){
    f.tmp(out1)->tbls1
    f.tmp(out2)->tbls2
    }
    else{
    f.tmp2(out1)->tbls1
    f.tmp2(out2)->tbls2
    }
	x.lab<-names(colMeans(tbls1,na.rm=T))
	pos1<-1
	which(x.lab=="0")->pos2
	x.lab[pos2]<-"TSS"
	pos3<-length(x.lab)

#max(colMeans(tbls1),colMeans(tbls2))
if(append==FALSE){
	if(frame==TRUE){
	plot(colMeans(tbls1,na.rm=T),ylim=ylim,ylab="Normalized read density",xlab="",xaxt="n",yaxt="n",type="l",main=main)
	points(colMeans(tbls2,na.rm=T),type="l",col="red")
	}else{
	plot(colMeans(tbls1,na.rm=T),ylim=ylim,ylab="Normalized read density",xlab="",xaxt="n",yaxt="n",type="l",axes=F,main=main)
	points(colMeans(tbls2,na.rm=T),type="l",col="red")
	}
	if (y.axis==TRUE){axis(2)}

}
else{
	points(colMeans(tbls1,na.rm=T),type="l")
	points(colMeans(tbls2,na.rm=T),type="l",col="red")
}

if (x.axis==TRUE){axis(1,c(pos1,pos2,pos3),x.lab[c(pos1,pos2,pos3)])}

	if(!is.null(legend)){
		legend("topright", legend = legend,  text.col = text.col)
	}

list(cvg.cutoff=cvg.cutoff,tbls1=tbls1,tbls2=tbls2,x.axis.pos=c(pos1,pos2,pos3),x.axis=x.lab[c(pos1,pos2,pos3)])->plot.out
return(plot.out)
}



plot.chip.single.around.TSS<-function(no.strand=TRUE,anno,chip1,input1,up=10000,down=10000,freq=200,cap.q = 0.99,s.width=100,frame=TRUE,main=NULL,legend=NULL,paste.chr=FALSE,ylim=c(0,1),gapped=FALSE,append=FALSE,text.col=c("black"),max.cvg=NULL,min.cvg=NULL,x.axis=TRUE,y.axis=TRUE,logratio=FALSE,input.type="BAM",anno.type="BED"){
	library(Repitools)
	library(GenomicRanges)
	if(input.type=="BAM"){
	if(gapped==FALSE){BAMpair2GRanges(chip1)->x1}else{BAM2GRanges(chip1)->x1}
	if(paste.chr==TRUE){paste("chr",seqlevels(x1),sep="")->seqlevels(x1)}

	if(gapped==FALSE){BAMpair2GRanges(input1)->x3}else{BAM2GRanges(input1)->x3}
	if(paste.chr==TRUE){paste("chr",seqlevels(x3),sep="")->seqlevels(x3)}

	}
	else if(input.type=="BED"){
	read.table(chip1)->tmp
	GRanges(tmp[,1],range=IRanges(start=tmp[,2],end=tmp[,3]))->x1
	read.table(input1)->tmp
	GRanges(tmp[,1],range=IRanges(start=tmp[,2],end=tmp[,3]))->x3
	}
	
		if(anno.type=="BED"){
		read.table(anno)->y
		GRanges(seqnames=y[,1],range=IRanges(start=y[,2],end=y[,3]),strand=y[,6])->region
		}
		else if(anno.type=="RANGE"){
			if(no.strand==TRUE){
				strand(anno)<-"+"
				anno->region
			}
			else{
				anno->region
		}	}

		fs1=featureScores(x1,region,up=up,down=down,freq=freq,verbose=TRUE,s.width=s.width)
		fs3=featureScores(x3,region,up=up,down=down,freq=freq,verbose=TRUE,s.width=s.width)

tables(fs1)[[1]]->tmp1
tables(fs3)[[1]]->tmp3

if(logratio==FALSE){
tmp1-tmp3->out1
}
else{
log2((tmp1+1e-10)/(tmp3+1e-10))->out1
}

if(is.null(max.cvg)){
max.cvg <- quantile(unlist(out1), cap.q, na.rm = TRUE)
min.cvg <- quantile(unlist(out1), 1-cap.q, na.rm = TRUE)
}
else{
max.cvg<-max.cvg
min.cvg<-min.cvg
}
print(max.cvg)
print(min.cvg)
list(max=max.cvg,min=min.cvg)->cvg.cutoff
    f.tmp<-function(x)
    {
        x[x > max.cvg] = max.cvg
        x[x < min.cvg] = min.cvg
        (x-min(x)) / (max(x)-min(x))
    }

    f.tmp2<-function(x){
        x[x > max.cvg] = max.cvg
        x[x < min.cvg] = min.cvg
	x
    }
    
    if(logratio==FALSE){
    f.tmp(out1)->tbls1
    }
    else{
    f.tmp2(out1)->tbls1
    }
	x.lab<-names(colMeans(tbls1,na.rm=T))
	pos1<-1
	which(x.lab=="0")->pos2
	x.lab[pos2]<-"TSS"
	pos3<-length(x.lab)

if(append==FALSE){
	if(frame==TRUE){
	plot(colMeans(tbls1,na.rm=T),ylim=ylim,ylab="Normalized read density",xlab="",xaxt="n",yaxt="n",type="l",main=main)
	}else{
	plot(colMeans(tbls1,na.rm=T),ylim=ylim,ylab="Normalized read density",xlab="",xaxt="n",yaxt="n",type="l",axes=F,main=main)
	}
	if (y.axis==TRUE){axis(2)}

}
else{
	points(colMeans(tbls1,na.rm=T),type="l")
}

if (x.axis==TRUE){axis(1,c(pos1,pos2,pos3),x.lab[c(pos1,pos2,pos3)])}

	if(!is.null(legend)){
		legend("topright", legend = legend,  text.col = text.col)
	}

list(cvg.cutoff=cvg.cutoff,tbls1=tbls1,x.axis.pos=c(pos1,pos2,pos3),x.axis=x.lab[c(pos1,pos2,pos3)])->plot.out
return(plot.out)
}


plot.chip.single.no.input.around.TSS<-function(anno,chip1,up=10000,down=10000,freq=200,cap.q = 0.99,s.width=100,frame=TRUE,main=NULL,legend=NULL,paste.chr=FALSE,ylim=c(0,1),gapped=FALSE,append=FALSE,text.col=c("black"),max.cvg=NULL,min.cvg=NULL,x.axis=TRUE,y.axis=TRUE,input.type="BAM",anno.type="BED"){
	library(Repitools)
	library(GenomicRanges)
	if(input.type=="BAM"){
	if(gapped==FALSE){BAMpair2GRanges(chip1)->x1}else{BAM2GRanges(chip1)->x1}
	if(paste.chr==TRUE){paste("chr",seqlevels(x1),sep="")->seqlevels(x1)}
	}
	else if(input.type=="BED"){
	read.table(chip1)->tmp
	GRanges(tmp[,1],range=IRanges(start=tmp[,2],end=tmp[,3]))->x1
	read.table(input1)->tmp
	GRanges(tmp[,1],range=IRanges(start=tmp[,2],end=tmp[,3]))->x3
	}
	
		if(anno.type=="BED"){
		read.table(anno)->y
		GRanges(seqnames=y[,1],range=IRanges(start=y[,2],end=y[,3]),strand=y[,6])->region
		}
		else if(anno.type=="RANGE"){
			anno->region
		}

		fs1=featureScores(x1,region,up=up,down=down,freq=freq,verbose=TRUE,s.width=s.width)

tables(fs1)[[1]]->tmp1

tmp1->out1

if(is.null(max.cvg)){
max.cvg <- quantile(c(unlist(out1)), cap.q, na.rm = TRUE)
}
else{
max.cvg<-max.cvg
}
print(max.cvg)

    f.tmp<-function(x)
    {
        x[x > max.cvg] = max.cvg
        x / max.cvg
    }

f.tmp(out1)->tbls1

print(max.cvg)

list(max=max.cvg)->cvg.cutoff



	x.lab<-names(colMeans(tbls1,na.rm=T))
	pos1<-1
	which(x.lab=="0")->pos2
	x.lab[pos2]<-"TSS"
	pos3<-length(x.lab)

if(append==FALSE){
	if(frame==TRUE){
	plot(colMeans(tbls1,na.rm=T),ylim=ylim,ylab="Normalized read density",xlab="",xaxt="n",yaxt="n",type="l",main=main)
	}else{
	plot(colMeans(tbls1,na.rm=T),ylim=ylim,ylab="Normalized read density",xlab="",xaxt="n",yaxt="n",type="l",axes=F,main=main)
	}
	if (y.axis==TRUE){axis(2)}

}
else{
	points(colMeans(tbls1,na.rm=T),type="l")
}

if (x.axis==TRUE){axis(1,c(pos1,pos2,pos3),x.lab[c(pos1,pos2,pos3)])}

	if(!is.null(legend)){
		legend("topright", legend = legend,  text.col = text.col)
	}

list(cvg.cutoff=cvg.cutoff,tbls1=tbls1,x.axis.pos=c(pos1,pos2,pos3),x.axis=x.lab[c(pos1,pos2,pos3)])->plot.out
return(plot.out)
}


############# using strand information
plot.chip.gene.single<-function(anno,chip1,input1,up=200000,down=200000,flanking.freq=5000,up.percent=100,down.percent=100,freq=1,flanking.dist="base",dist="percent",cap.q = 0.99,s.width=5000,frame=TRUE,main=NULL,legend=NULL,paste.chr=FALSE,ylim=c(0,1),gapped=FALSE,append=FALSE,text.col=c("black"),max.cvg=NULL,min.cvg=NULL,x.axis=TRUE,y.axis=TRUE,logratio=FALSE){
	require(Repitools)
	require(GenomicRanges)
	if(gapped==FALSE){BAMpair2GRanges(chip1)->x1}else{BAM2GRanges(chip1)->x1}
	if(paste.chr==TRUE){paste("chr",seqlevels(x1),sep="")->seqlevels(x1)}
	if(gapped==FALSE){BAMpair2GRanges(input1)->x3}else{BAM2GRanges(input1)->x3}
	if(paste.chr==TRUE){paste("chr",seqlevels(x3),sep="")->seqlevels(x3)}


	read.table(anno)->y
	GRanges(seqnames=y[,1],range=IRanges(start=y[,2],end=y[,3]),strand=y[,6])->region
	if(flanking.dist=="base"){
	fs1.middle=featureScores(x1,region,up=0,down=100,freq=freq,verbose=TRUE,dist=dist,s.width=s.width)
	fs1.left=featureScores(x1,region,up=up,down=0,freq=flanking.freq,verbose=TRUE,dist=flanking.dist,s.width=s.width)
	fs3.middle=featureScores(x3,region,up=0,down=100,freq=freq,verbose=TRUE,dist=dist,s.width=s.width)
	fs3.left=featureScores(x3,region,up=up,down=0,freq=flanking.freq,verbose=TRUE,dist=flanking.dist,s.width=s.width)
	strand(region)=="-" ->ind.minus
	strand(region)=="+" ->ind.plus
	strand(region)[ind.minus]<-"+"
	strand(region)[ind.plus]<-"-"
	fs1.right=featureScores(x1,region,up=down,down=0,freq=flanking.freq,verbose=TRUE,dist=flanking.dist,s.width=s.width)
	fs3.right=featureScores(x3,region,up=down,down=0,freq=flanking.freq,verbose=TRUE,dist=flanking.dist,s.width=s.width)
	}
	else if(flanking.dist=="percent"){
	fs1=featureScores(x1,region,up=up.percent,down=down.percent,freq=freq,verbose=TRUE,dist=dist,s.width=s.width)
	fs3=featureScores(x3,region,up=up.percent,down=down.percent,freq=freq,verbose=TRUE,dist=dist,s.width=s.width)
	}

	

if(flanking.dist=="base"){
cbind(tables(fs1.left)[[1]],tables(fs1.middle)[[1]],t(apply(tables(fs1.right)[[1]],1,rev)))->tmp1
cbind(tables(fs3.left)[[1]],tables(fs3.middle)[[1]],t(apply(tables(fs3.right)[[1]],1,rev)))->tmp3
}
else if(flanking.dist=="percent"){
tables(fs1)[[1]]->tmp1
tables(fs3)[[1]]->tmp3
}

if(logratio==FALSE){
tmp1-tmp3->out1
}
else{
log2((tmp1+1e-10)/(tmp3+1e-10))->out1
}

if(is.null(max.cvg)){
max.cvg <- quantile(c(unlist(out1)), cap.q, na.rm = TRUE)
min.cvg <- quantile(c(unlist(out1)), 1-cap.q, na.rm = TRUE)
}
else{
max.cvg<-max.cvg
min.cvg<-min.cvg
}
print(max.cvg)
print(min.cvg)
list(max=max.cvg,min=min.cvg)->cvg.cutoff
    f.tmp<-function(x)
    {
        x[x > max.cvg] = max.cvg
        x[x < min.cvg] = min.cvg
        (x-min(x)) / (max(x)-min(x))
    }

    f.tmp2<-function(x){
        x[x > max.cvg] = max.cvg
        x[x < min.cvg] = min.cvg
	x
    }

    if(logratio==FALSE){
    f.tmp(out1)->tbls1
    }
    else{
    f.tmp2(out1)->tbls1
    }
	x.lab<-names(colMeans(tbls1,na.rm=T))
	pos1<-1
	which(x.lab=="0 %")->pos2
	which(x.lab=="100 %")->pos3
	pos4<-length(x.lab)

if(append==FALSE){
	if(frame==TRUE){
	plot(colMeans(tbls1,na.rm=T),ylim=ylim,ylab="Normalized read density",xlab="",xaxt="n",yaxt="n",type="l",main=main,col=text.col[1])
	}else{
	plot(colMeans(tbls1,na.rm=T),ylim=ylim,ylab="Normalized read density",xlab="",xaxt="n",yaxt="n",type="l",axes=F,main=main,col=text.col[1])
	}
	if (y.axis==TRUE){axis(2)}

}
else{
	points(colMeans(tbls1,na.rm=T),type="l",text.col=text.col[1])
}

if (x.axis==TRUE){axis(1,c(pos1,pos2,pos3,pos4),x.lab[c(pos1,pos2,pos3,pos4)])}

	if(!is.null(legend)){
		legend("topright", legend = legend,  text.col = text.col)
	}

list(cvg.cutoff=cvg.cutoff,tbls1=tbls1,x.axis.pos=c(pos1,pos2,pos3,pos4),x.axis=x.lab[c(pos1,pos2,pos3,pos4)])->plot.out
return(plot.out)
}
#####################################



############# using strand information
plot.chip.gene.comp<-function(anno,anno2=NULL,chip1,chip2,input1,input2,up=200000,down=200000,flanking.freq=5000,up.percent=100,down.percent=100,freq=1,flanking.dist="base",dist="percent",cap.q = 0.99,s.width=5000,frame=TRUE,main=NULL,legend=NULL,paste.chr=FALSE,ylim=c(0,1),gapped=FALSE,append=FALSE,text.col=c("black","red"),max.cvg=NULL,min.cvg=NULL,x.axis=TRUE,y.axis=TRUE,logratio=FALSE){
	require(Repitools)
	require(GenomicRanges)
	if(gapped==FALSE){BAMpair2GRanges(chip1)->x1}else{BAM2GRanges(chip1)->x1}
	if(paste.chr==TRUE){paste("chr",seqlevels(x1),sep="")->seqlevels(x1)}
	if(gapped==FALSE){BAMpair2GRanges(chip2)->x2}else{BAM2GRanges(chip2)->x2}
	if(paste.chr==TRUE){paste("chr",seqlevels(x2),sep="")->seqlevels(x2)}
	if(gapped==FALSE){BAMpair2GRanges(input1)->x3}else{BAM2GRanges(input1)->x3}
	if(paste.chr==TRUE){paste("chr",seqlevels(x3),sep="")->seqlevels(x3)}
	if(gapped==FALSE){BAMpair2GRanges(input2)->x4}else{BAM2GRanges(input2)->x4}
	if(paste.chr==TRUE){paste("chr",seqlevels(x4),sep="")->seqlevels(x4)}
	
	
	read.table(anno)->y
	GRanges(seqnames=y[,1],range=IRanges(start=y[,2],end=y[,3]),strand=y[,6])->region
	if(flanking.dist=="base"){
	fs1.middle=featureScores(x1,region,up=0,down=100,freq=freq,verbose=TRUE,dist=dist,s.width=s.width)
	fs1.left=featureScores(x1,region,up=up,down=0,freq=flanking.freq,verbose=TRUE,dist=flanking.dist,s.width=s.width)
	fs3.middle=featureScores(x3,region,up=0,down=100,freq=freq,verbose=TRUE,dist=dist,s.width=s.width)
	fs3.left=featureScores(x3,region,up=up,down=0,freq=flanking.freq,verbose=TRUE,dist=flanking.dist,s.width=s.width)
	strand(region)=="-" ->ind.minus
	strand(region)=="+" ->ind.plus
	strand(region)[ind.minus]<-"+"
	strand(region)[ind.plus]<-"-"
	fs1.right=featureScores(x1,region,up=down,down=0,freq=flanking.freq,verbose=TRUE,dist=flanking.dist,s.width=s.width)
	fs3.right=featureScores(x3,region,up=down,down=0,freq=flanking.freq,verbose=TRUE,dist=flanking.dist,s.width=s.width)
	}
	else if(flanking.dist=="percent"){
	fs1=featureScores(x1,region,up=up.percent,down=down.percent,freq=freq,verbose=TRUE,dist=dist,s.width=s.width)
	fs3=featureScores(x3,region,up=up.percent,down=down.percent,freq=freq,verbose=TRUE,dist=dist,s.width=s.width)
	}

	if(flanking.dist=="base"){
	if(is.null(anno2)){
		read.table(anno)->y
		GRanges(seqnames=y[,1],range=IRanges(start=y[,2],end=y[,3]),strand=y[,6])->region
		fs2.middle=featureScores(x2,region,up=0,down=100,freq=freq,verbose=TRUE,dist=dist,s.width=s.width)
		fs2.left=featureScores(x2,region,up=up,down=0,freq=flanking.freq,verbose=TRUE,dist=flanking.dist,s.width=s.width)
		fs4.middle=featureScores(x4,region,up=0,down=100,freq=freq,verbose=TRUE,dist=dist,s.width=s.width)
		fs4.left=featureScores(x4,region,up=up,down=0,freq=flanking.freq,verbose=TRUE,dist=flanking.dist,s.width=s.width)
		strand(region)=="-" ->ind.minus
		strand(region)=="+" ->ind.plus
		strand(region)[ind.minus]<-"+"
		strand(region)[ind.plus]<-"-"
		fs2.right=featureScores(x2,region,up=down,down=0,freq=flanking.freq,verbose=TRUE,dist=flanking.dist,s.width=s.width)
		fs4.right=featureScores(x4,region,up=down,down=0,freq=flanking.freq,verbose=TRUE,dist=flanking.dist,s.width=s.width)
	}
	else{
		read.table(anno2)->y2
		GRanges(seqnames=y2[,1],range=IRanges(start=y2[,2],end=y2[,3]),strand=y2[,6])->region2
		fs2.middle=featureScores(x2,region2,up=0,down=100,freq=freq,verbose=TRUE,dist=dist,s.width=s.width)
		fs2.left=featureScores(x2,region2,up=up,down=0,freq=flanking.freq,verbose=TRUE,dist=flanking.dist,s.width=s.width)
		fs4.middle=featureScores(x4,region2,up=0,down=100,freq=freq,verbose=TRUE,dist=dist,s.width=s.width)
		fs4.left=featureScores(x4,region2,up=up,down=0,freq=flanking.freq,verbose=TRUE,dist=flanking.dist,s.width=s.width)
		strand(region2)=="-" ->ind.minus
		strand(region2)=="+" ->ind.plus
		strand(region2)[ind.minus]<-"+"
		strand(region2)[ind.plus]<-"-"
		fs2.right=featureScores(x2,region2,up=down,down=0,freq=flanking.freq,verbose=TRUE,dist=flanking.dist,s.width=s.width)
		fs4.right=featureScores(x4,region2,up=down,down=0,freq=flanking.freq,verbose=TRUE,dist=flanking.dist,s.width=s.width)
	}
	}
	else if(flanking.dist=="percent"){
	if(is.null(anno2)){
		read.table(anno)->y
		GRanges(seqnames=y[,1],range=IRanges(start=y[,2],end=y[,3]),strand=y[,6])->region
		fs2=featureScores(x2,region,up=up.percent,down=down.percent,freq=freq,verbose=TRUE,dist=dist,s.width=s.width)
		fs4=featureScores(x4,region,up=up.percent,down=down.percent,freq=freq,verbose=TRUE,dist=dist,s.width=s.width)
	}
	else{
		read.table(anno2)->y2
		GRanges(seqnames=y2[,1],range=IRanges(start=y2[,2],end=y2[,3]),strand=y2[,6])->region2
		fs2=featureScores(x2,region2,up=up.percent,down=down.percent,freq=freq,verbose=TRUE,dist=dist,s.width=s.width)
		fs4=featureScores(x4,region2,up=up.percent,down=down.percent,freq=freq,verbose=TRUE,dist=dist,s.width=s.width)
	}
	}
	
if(flanking.dist=="base"){	
cbind(tables(fs1.left)[[1]],tables(fs1.middle)[[1]],t(apply(tables(fs1.right)[[1]],1,rev)))->tmp1
cbind(tables(fs2.left)[[1]],tables(fs2.middle)[[1]],t(apply(tables(fs2.right)[[1]],1,rev)))->tmp2
cbind(tables(fs3.left)[[1]],tables(fs3.middle)[[1]],t(apply(tables(fs3.right)[[1]],1,rev)))->tmp3
cbind(tables(fs4.left)[[1]],tables(fs4.middle)[[1]],t(apply(tables(fs4.right)[[1]],1,rev)))->tmp4
}
else if(flanking.dist=="percent"){
tables(fs1)[[1]]->tmp1
tables(fs2)[[1]]->tmp2
tables(fs3)[[1]]->tmp3
tables(fs4)[[1]]->tmp4
}

if(logratio==FALSE){
tmp1-tmp3->out1
tmp2-tmp4->out2
}
else{
log2((tmp1+1e-10)/(tmp3+1e-10))->out1
log2((tmp2+1e-10)/(tmp4+1e-10))->out2
}

if(is.null(max.cvg)){
max.cvg <- quantile(c(unlist(out1),unlist(out2)), cap.q, na.rm = TRUE)
min.cvg <- quantile(c(unlist(out1),unlist(out2)), 1-cap.q, na.rm = TRUE)
}
else{
max.cvg<-max.cvg
min.cvg<-min.cvg
}
print(max.cvg)
print(min.cvg)
list(max=max.cvg,min=min.cvg)->cvg.cutoff
    f.tmp<-function(x)
    {
        x[x > max.cvg] = max.cvg
        x[x < min.cvg] = min.cvg
        (x-min(x)) / (max(x)-min(x))
    }

    f.tmp2<-function(x){
        x[x > max.cvg] = max.cvg
        x[x < min.cvg] = min.cvg
	x
    }
    
    if(logratio==FALSE){
    f.tmp(out1)->tbls1
    f.tmp(out2)->tbls2
    }
    else{
    f.tmp2(out1)->tbls1
    f.tmp2(out2)->tbls2
    }
	x.lab<-names(colMeans(tbls1,na.rm=T))
	pos1<-1
	which(x.lab=="0 %")->pos2
	which(x.lab=="100 %")->pos3
	pos4<-length(x.lab)

#max(colMeans(tbls1),colMeans(tbls2))
if(append==FALSE){
	if(frame==TRUE){
	plot(colMeans(tbls1,na.rm=T),ylim=ylim,ylab="Normalized read density",xlab="",xaxt="n",yaxt="n",type="l",main=main,col=text.col[1])
	points(colMeans(tbls2,na.rm=T),type="l",col=text.col[2])
	}else{
	plot(colMeans(tbls1,na.rm=T),ylim=ylim,ylab="Normalized read density",xlab="",xaxt="n",yaxt="n",type="l",axes=F,main=main,col=text.col[1])
	points(colMeans(tbls2,na.rm=T),type="l",col=text.col[2])
	}
	if (y.axis==TRUE){axis(2)}

}
else{
	points(colMeans(tbls1,na.rm=T),type="l",text.col=text.col[1])
	points(colMeans(tbls2,na.rm=T),type="l",col=text.col[2])
}

if (x.axis==TRUE){axis(1,c(pos1,pos2,pos3,pos4),x.lab[c(pos1,pos2,pos3,pos4)])}

	if(!is.null(legend)){
		legend("topright", legend = legend,  text.col = text.col)
	}

list(cvg.cutoff=cvg.cutoff,tbls1=tbls1,tbls2=tbls2,x.axis.pos=c(pos1,pos2,pos3,pos4),x.axis=x.lab[c(pos1,pos2,pos3,pos4)])->plot.out
return(plot.out)
}
#####################################

plot.chip.single<-function(no.strand=TRUE,anno,chip1,input1,up=200000,down=200000,flanking.freq=5000,freq=1,flanking.dist="base",dist="percent",max.cvg=NULL,min.cvg=NULL,s.width=5000,frame=TRUE,main=NULL,legend=NULL,paste.chr=FALSE,ylim=c(0,1),gapped=FALSE,append=FALSE,col="black",text.col="black",cap.q=0.99,input.type="BAM",anno.type="RANGE",logratio=FALSE){
	library(Repitools)
	library(GenomicRanges)
	if(input.type=="BAM"){
		if(gapped==FALSE){BAMpair2GRanges(chip1)->x1}else{BAM2GRanges(chip1)->x1}
		if(paste.chr==TRUE){paste("chr",seqlevels(x1),sep="")->seqlevels(x1)}

		if(gapped==FALSE){BAMpair2GRanges(input1)->x3}else{BAM2GRanges(input1)->x3}
		if(paste.chr==TRUE){paste("chr",seqlevels(x3),sep="")->seqlevels(x3)}
	}
	else if(input.type=="BED"){
		read.table(chip1)->tmp
		GRanges(tmp[,1],range=IRanges(start=tmp[,2],end=tmp[,3]))->x1
		read.table(input1)->tmp
		GRanges(tmp[,1],range=IRanges(start=tmp[,2],end=tmp[,3]))->x3
	}

	if(anno.type=="BED"){
	read.table(anno)->y
		if(no.strand==TRUE){
			GRanges(seqnames=y[,1],range=IRanges(start=y[,2],end=y[,3]),strand="+")->region
		}
		else{
			GRanges(seqnames=y[,1],range=IRanges(start=y[,2],end=y[,3]),strand=y[,6])->region
		}
	}
	else if(anno.type=="RANGE"){
		if(no.strand==TRUE){
			strand(anno)<-"+"
			anno->region
		}
		else{
			anno->region
		}
	}

	fs1.middle=featureScores(x1,region,up=0,down=100,freq=freq,verbose=TRUE,dist=dist,s.width=s.width)
	fs1.left=featureScores(x1,region,up=up,down=0,freq=flanking.freq,verbose=TRUE,dist=flanking.dist,s.width=s.width)
	fs3.middle=featureScores(x3,region,up=0,down=100,freq=freq,verbose=TRUE,dist=dist,s.width=s.width)
	fs3.left=featureScores(x3,region,up=up,down=0,freq=flanking.freq,verbose=TRUE,dist=flanking.dist,s.width=s.width)
	if(no.strand==TRUE){
		strand(region)<-"-"
	}
	else if(no.strand==FALSE){
		strand(region)=="+" -> ind1
		strand(region)=="-" -> ind2
		strand(region[ind1])<- "-"
		strand(region[ind2])<- "+"
	}
	fs1.right=featureScores(x1,region,up=down,down=0,freq=flanking.freq,verbose=TRUE,dist=flanking.dist,s.width=s.width)
	fs3.right=featureScores(x3,region,up=down,down=0,freq=flanking.freq,verbose=TRUE,dist=flanking.dist,s.width=s.width)

cbind(tables(fs1.left)[[1]],tables(fs1.middle)[[1]],t(apply(tables(fs1.right)[[1]],1,rev)))->tmp1
cbind(tables(fs3.left)[[1]],tables(fs3.middle)[[1]],t(apply(tables(fs3.right)[[1]],1,rev)))->tmp3

if(logratio==FALSE){
tmp1-tmp3->out1
}
else{
log2((tmp1+1e-10)/(tmp3+1e-10))->out1
}

if(is.null(max.cvg)){
max.cvg <- quantile(c(unlist(out1)), cap.q, na.rm = TRUE)
min.cvg <- quantile(c(unlist(out1)), 1-cap.q, na.rm = TRUE)
}
else{
max.cvg<-max.cvg
min.cvg<-min.cvg
}
print(max.cvg)
print(min.cvg)

    f.tmp<-function(x)
    {
        x[x > max.cvg] = max.cvg
        x[x < min.cvg] = min.cvg
        (x-min(x)) / (max(x)-min(x))
    }

    f.tmp2<-function(x){
        x[x > max.cvg] = max.cvg
        x[x < min.cvg] = min.cvg
	x
    }

    if(logratio==FALSE){
    f.tmp(out1)->tbls1
    }
    else{
    f.tmp2(out1)->tbls1
    }
    
if(append==FALSE){
	if(frame==TRUE){
	plot(colMeans(tbls1),ylim=ylim,ylab="Normalized read density",xlab="relative position",xaxt="n",type="l",main=main,col=col)
	x.lab<-names(colMeans(tbls1))
	pos1<-1
	which(x.lab=="0 %")->pos2
	which(x.lab=="100 %")->pos3
	pos4<-length(x.lab)
	axis(1,c(pos1,pos2,pos3,pos4),x.lab[c(pos1,pos2,pos3,pos4)])
	}else{
	plot(colMeans(tbls1),ylim=ylim,ylab="",xlab="",xaxt="n",yaxt="n",type="l",main=main,col=col)
	}
}else{
	points(colMeans(tbls1),type="l",main=main,col=col)	
}
	if(!is.null(legend)){
		legend("topright", legend = legend,  text.col = text.col)
	}
	
	x.lab<-names(colMeans(tbls1))
	pos1<-1
	which(x.lab=="0 %")->pos2
	which(x.lab=="100 %")->pos3
	pos4<-length(x.lab)

list(max=max.cvg,min=min.cvg)->cvg.cutoff
list(cvg.cutoff=cvg.cutoff,tbls1=tbls1,x.axis.pos=c(pos1,pos2,pos3,pos4),x.axis=x.lab[c(pos1,pos2,pos3,pos4)])->plot.out
return(plot.out)
}


###################### only calculate chip enrichment according region
##################################### simple function
cal.rpkm.chip<-function(read,range,input="bam",pair=F){
require(GenomicRanges)
require(Repitools)
if(input=="bam"){
	if(pair==F){BAM2GRanges(read)->chip}
	else if(pair==T){BAMpair2GRanges(read)->chip}
}
else if(input=="bed"){
read.table(read)->tmp
GRanges(seqname=tmp$V1,range=IRanges(start=tmp$V2,end=tmp$V3))->chip
}

countOverlaps(range,chip)->tmp
tmp/(length(chip)/1000000)/(width(range)/1000)->rpkm
mcols(range)$rpkm<-rpkm
return(range)
}

###################################
calculate.chip.single<-function(anno,chip1,input1,max.cvg=NULL,min.cvg=NULL,paste.chr=FALSE,gapped=FALSE,cap.q=0.99,input.type="BAM",anno.type="RANGE",logratio=FALSE){
	require(GenomicRanges)
	require(Repitools)
	if(input.type=="BAM"){
		if(gapped==FALSE){BAMpair2GRanges(chip1)->x1}else{BAM2GRanges(chip1)->x1}
		if(paste.chr==TRUE){paste("chr",seqlevels(x1),sep="")->seqlevels(x1)}

		if(gapped==FALSE){BAMpair2GRanges(input1)->x3}else{BAM2GRanges(input1)->x3}
		if(paste.chr==TRUE){paste("chr",seqlevels(x3),sep="")->seqlevels(x3)}
	}
	else if(input.type=="BED"){
		read.table(chip1)->tmp
		GRanges(tmp[,1],range=IRanges(start=tmp[,2],end=tmp[,3]))->x1
		read.table(input1)->tmp
		GRanges(tmp[,1],range=IRanges(start=tmp[,2],end=tmp[,3]))->x3
	}

	if(anno.type=="BED"){
	read.table(anno)->y
		if(no.strand==TRUE){
			GRanges(seqnames=y[,1],range=IRanges(start=y[,2],end=y[,3]),strand="+")->region
		}
		else{
			GRanges(seqnames=y[,1],range=IRanges(start=y[,2],end=y[,3]),strand=y[,6])->region
		}
	}
	else if(anno.type=="RANGE"){
		anno->region
	}

	countOverlaps(region,x1)/width(region)/length(x1)*1000000->tmp1
	countOverlaps(region,x3)/width(region)/length(x3)*1000000->tmp3

if(logratio==FALSE){
tmp1-tmp3->out1
}
else{
log2((tmp1+1e-10)/(tmp3+1e-10))->out1
}

if(is.null(max.cvg)){
max.cvg <- quantile(c(unlist(out1)), cap.q, na.rm = TRUE)
min.cvg <- quantile(c(unlist(out1)), 1-cap.q, na.rm = TRUE)
}
else{
max.cvg<-max.cvg
min.cvg<-min.cvg
}
print(max.cvg)
print(min.cvg)

    f.tmp<-function(x)
    {
        x[x > max.cvg] = max.cvg
        x[x < min.cvg] = min.cvg
        (x-min(x)) / (max(x)-min(x))
    }

    f.tmp2<-function(x){
        x[x > max.cvg] = max.cvg
        x[x < min.cvg] = min.cvg
	x
    }

    if(logratio==FALSE){
    f.tmp(out1)->tbls1
    }
    else{
    f.tmp2(out1)->tbls1
    }
    

list(max=max.cvg,min=min.cvg)->cvg.cutoff
list(cvg.cutoff=cvg.cutoff,tbls1=tbls1)->plot.out
return(plot.out)
}

#############################
calculate.chip.single.no.input<-function(anno,chip1,max.cvg=NULL,min.cvg=NULL,paste.chr=FALSE,gapped=FALSE,cap.q=0.99,input.type="BAM",anno.type="RANGE"){
	require(GenomicRanges)
	require(Repitools)
	if(input.type=="BAM"){
		if(gapped==FALSE){BAMpair2GRanges(chip1)->x1}else{BAM2GRanges(chip1)->x1}
		if(paste.chr==TRUE){paste("chr",seqlevels(x1),sep="")->seqlevels(x1)}
	}
	else if(input.type=="BED"){
		read.table(chip1)->tmp
		GRanges(tmp[,1],range=IRanges(start=tmp[,2],end=tmp[,3]))->x1
	}

	if(anno.type=="BED"){
	read.table(anno)->y
		if(no.strand==TRUE){
			GRanges(seqnames=y[,1],range=IRanges(start=y[,2],end=y[,3]),strand="+")->region
		}
		else{
			GRanges(seqnames=y[,1],range=IRanges(start=y[,2],end=y[,3]),strand=y[,6])->region
		}
	}
	else if(anno.type=="RANGE"){
		anno->region
	}

	countOverlaps(region,x1)/width(region)/length(x1)*1000000->tmp1

tmp1->out1

if(is.null(max.cvg)){
max.cvg <- quantile(c(unlist(out1)), cap.q, na.rm = TRUE)
min.cvg <- quantile(c(unlist(out1)), 1-cap.q, na.rm = TRUE)
}
else{
max.cvg<-max.cvg
min.cvg<-min.cvg
}
print(max.cvg)
print(min.cvg)

    f.tmp2<-function(x){
        x[x > max.cvg] = max.cvg
        x[x < min.cvg] = min.cvg
	x
    }

    f.tmp2(out1)->tbls1

list(max=max.cvg,min=min.cvg)->cvg.cutoff
list(cvg.cutoff=cvg.cutoff,tbls1=tbls1)->plot.out
return(plot.out)
}

BAMpair2GRanges<-function(path, what = c("rname", "strand", "pos", "isize"),
             flag = scanBamFlag(isUnmappedQuery = FALSE, isDuplicate = FALSE,isPaired=TRUE,isProperPair=TRUE,hasUnmappedMate=FALSE,isMinusStrand=FALSE,isMateMinusStrand=TRUE),verbose = TRUE)
{
    require(Rsamtools)
    if(length(path) > 1)
	stop("This method is only for one BAM file. See ?BAM2GRangesList.")

    if(verbose == TRUE)
	cat("Reading BAM file ", basename(path), ".\n", sep = '')
    filters <- ScanBamParam(what = what, flag = flag)
    bam <- scanBam(path, param = filters)[[1]]
    anno <- !names(bam) %in% c("rname", "strand", "pos", "isize")
    if(verbose == TRUE) cat("Creating GRanges from BAM datalist.\n")
    gr <- GRanges(bam$rname, IRanges(start=bam$pos, end = bam$pos+bam$isize-1),
                  bam$strand, bam[anno])
}

#####################################

plot.chip.comp.no.input<-function(no.strand=TRUE,anno.type="BED",anno,anno2=NULL,chip1,chip2,up=200000,down=200000,flanking.freq=5000,freq=1,flanking.dist="base",dist="percent",cap.q = 0.99,s.width=5000,frame=TRUE,append=FALSE,main=NULL,legend=NULL,paste.chr=FALSE,ylim=c(0,1),gapped=FALSE,text.col=c("black","red"),max.cvg=NULL){
	library(Repitools)
	library(GenomicRanges)
	paste(chip1,".rda",sep="")->file
	if(file.exists(file)){
		print("loading bam file")
		print (file)
		load(file);bam->x1
	}
	else{
		if(gapped==FALSE){BAMpair2GRanges(chip1)->x1}else{BAM2GRanges(chip1)->x1}
		if(paste.chr==TRUE){paste("chr",seqlevels(x1),sep="")->seqlevels(x1)}
		x1->bam
		print("save bam...")
		save(bam,file=file)
	}

	paste(chip2,".rda",sep="")->file

	if(file.exists(file)){
		print("loading bam")
		print (file)
		load(file);bam->x2
	}
	else{
		if(gapped==FALSE){BAMpair2GRanges(chip2)->x2}else{BAM2GRanges(chip2)->x2}
		if(paste.chr==TRUE){paste("chr",seqlevels(x2),sep="")->seqlevels(x2)}
		x2->bam
		print("save bam...")
		save(bam,file=file)
	}
	
	if(anno.type=="BED"){
	read.table(anno)->y
		if(no.strand==TRUE){
			GRanges(seqnames=y[,1],range=IRanges(start=y[,2],end=y[,3]),strand="+")->region
		}
		else{
			GRanges(seqnames=y[,1],range=IRanges(start=y[,2],end=y[,3]),strand=y[,6])->region
		}
	}
	else if(anno.type=="RANGE"){
		if(no.strand==TRUE){
			strand(anno)<-"+"
			anno->region
		}
		else{
			anno->region
		}
	}

	if(flanking.dist=="base"){
	fs1.middle=featureScores(x1,region,up=0,down=100,freq=freq,verbose=TRUE,dist=dist,s.width=s.width)
	fs1.left=featureScores(x1,region,up=up,down=0,freq=flanking.freq,verbose=TRUE,dist=flanking.dist,s.width=s.width)
		if(no.strand==TRUE){
		strand(region)<-"-"
	}
	else if(no.strand==FALSE){
		strand(region)=="+" -> ind1
		strand(region)=="-" -> ind2
		strand(region[ind1])<- "-"
		strand(region[ind2])<- "+"
	}

	fs1.right=featureScores(x1,region,up=down,down=0,freq=flanking.freq,verbose=TRUE,dist=flanking.dist,s.width=s.width)
	}
	else if(flanking.dist=="percent"){
	fs1=featureScores(x1,region,up=100,down=200,freq=freq,verbose=TRUE,dist=dist,s.width=s.width)
	fs3=featureScores(x3,region,up=100,down=200,freq=freq,verbose=TRUE,dist=dist,s.width=s.width)
	}

	if(flanking.dist=="base"){
	if(is.null(anno2)){
		fs2.middle=featureScores(x2,region,up=0,down=100,freq=freq,verbose=TRUE,dist=dist,s.width=s.width)
		fs2.left=featureScores(x2,region,up=up,down=0,freq=flanking.freq,verbose=TRUE,dist=flanking.dist,s.width=s.width)
	if(no.strand==TRUE){
		strand(region)<-"-"
	}
	else if(no.strand==FALSE){
		strand(region)=="+" -> ind1
		strand(region)=="-" -> ind2
		strand(region[ind1])<- "-"
		strand(region[ind2])<- "+"
	}
	fs2.right=featureScores(x2,region,up=down,down=0,freq=flanking.freq,verbose=TRUE,dist=flanking.dist,s.width=s.width)
	}
	else{
	#read.table(anno2)->y2
	#GRanges(seqnames=y2[,1],range=IRanges(start=y2[,2],end=y2[,3]),strand="+")->region2

	if(anno.type=="BED"){
	read.table(anno2)->y2
		if(no.strand==TRUE){
			GRanges(seqnames=y2[,1],range=IRanges(start=y2[,2],end=y2[,3]),strand="+")->region2
		}
		else{
			GRanges(seqnames=y2[,1],range=IRanges(start=y2[,2],end=y2[,3]),strand=y2[,6])->region2
		}
	}
	else if(anno.type=="RANGE"){
		if(no.strand==TRUE){
			strand(anno)<-"+"
			anno2->region2
		}
		else{
			anno2->region2
		}
	}

		fs2.middle=featureScores(x2,region2,up=0,down=100,freq=freq,verbose=TRUE,dist=dist,s.width=s.width)
		fs2.left=featureScores(x2,region2,up=up,down=0,freq=flanking.freq,verbose=TRUE,dist=flanking.dist,s.width=s.width)
		strand(region2)<-"-"
		fs2.right=featureScores(x2,region2,up=down,down=0,freq=flanking.freq,verbose=TRUE,dist=flanking.dist,s.width=s.width)
	}
	}
	else if(flanking.dist=="percent"){
	if(is.null(anno2)){
		read.table(anno)->y
		GRanges(seqnames=y[,1],range=IRanges(start=y[,2],end=y[,3]),strand="+")->region
		fs2=featureScores(x2,region,up=100,down=200,freq=freq,verbose=TRUE,dist=dist,s.width=s.width)
	}
	else{
		read.table(anno2)->y2
		GRanges(seqnames=y2[,1],range=IRanges(start=y2[,2],end=y2[,3]),strand="+")->region2
		fs2=featureScores(x2,region2,up=100,down=200,freq=freq,verbose=TRUE,dist=dist,s.width=s.width)
	}
	}
	
if(flanking.dist=="base"){	
cbind(tables(fs1.left)[[1]],tables(fs1.middle)[[1]],t(apply(tables(fs1.right)[[1]],1,rev)))->tmp1
cbind(tables(fs2.left)[[1]],tables(fs2.middle)[[1]],t(apply(tables(fs2.right)[[1]],1,rev)))->tmp2
}
else if(flanking.dist=="percent"){
tables(fs1)[[1]]->tmp1
tables(fs2)[[1]]->tmp2
}
tmp1->out1
tmp2->out2
if(is.null(max.cvg)){
max.cvg <- quantile(c(unlist(out1),unlist(out2)), cap.q, na.rm = TRUE)
}
print(max.cvg)

    f.tmp<-function(x)
    {
        x[x > max.cvg] = max.cvg
        x /max.cvg
    }
    f.tmp(out1)->tbls1
    f.tmp(out2)->tbls2

	x.lab<-names(colMeans(tbls1))
	pos1<-1
	which(x.lab=="0 %")->pos2
	which(x.lab=="100 %")->pos3
	pos4<-length(x.lab)

#max(colMeans(tbls1),colMeans(tbls2))
if(append==FALSE){
	if(frame==TRUE){
	plot(colMeans(tbls1),ylim=ylim,ylab="Normalized read density",xlab="",xaxt="n",type="l",main=main,col=text.col[1])
	axis(1,c(pos1,pos2,pos3,pos4),x.lab[c(pos1,pos2,pos3,pos4)])
	points(colMeans(tbls2),type="l",col=text.col[2])
	}else{
	plot(colMeans(tbls1),ylim=ylim,ylab="",xlab="",xaxt="n",yaxt="n",type="l",main=main,axes=F,col=text.col[1])
	axis(1,c(pos1,pos2,pos3,pos4),x.lab[c(pos1,pos2,pos3,pos4)])
	axis(2)
	points(colMeans(tbls2),type="l",col=text.col[2])
	}
}else{
	points(colMeans(tbls1),type="l",col=text.col[1])
	points(colMeans(tbls2),type="l",col=text.col[2])	
}


	if(!is.null(legend)){
		legend("topright", legend = legend,  text.col = text.col)
	}
list(max=max.cvg,tbls1=tbls1,tbls2=tbls2)->out
return(out)
}

#####################################

plot.chip.single.no.input<-function(no.strand=False,anno.type ="BED",anno,chip1,input1,up=200000,down=200000,flanking.freq=5000,freq=1,flanking.dist="base",dist="percent",max.cvg=NULL,min.cvg=NULL,s.width=5000,frame=TRUE,main=NULL,legend=NULL,paste.chr=FALSE,ylim=c(0,1),gapped=FALSE,append=FALSE,col="black",text.col="black",cap.q=0.99){
	library(Repitools)
	library(GenomicRanges)

	paste(chip1,".rda",sep="")->file
	if(file.exists(file)){
		print("loading bam")
		load(file);bam->x1
	}
	else{
		if(gapped==FALSE){BAMpair2GRanges(chip1)->x1}else{BAM2GRanges(chip1)->x1}
		if(paste.chr==TRUE){paste("chr",seqlevels(x1),sep="")->seqlevels(x1)}
		x1->bam
		print("save bam...")
		save(bam,file=file)
	}
	
	if(anno.type=="BED"){
	read.table(anno)->y
		if(no.strand==TRUE){
			GRanges(seqnames=y[,1],range=IRanges(start=y[,2],end=y[,3]),strand="+")->region
		}
		else{
			GRanges(seqnames=y[,1],range=IRanges(start=y[,2],end=y[,3]),strand=y[,6])->region
		}
	}
	else if(anno.type=="RANGE"){
		if(no.strand==TRUE){
			strand(anno)<-"+"
			anno->region
		}
		else{
			anno->region
		}
	}
	
	fs1.middle=featureScores(x1,region,up=0,down=100,freq=freq,verbose=TRUE,dist=dist,s.width=s.width)
	fs1.left=featureScores(x1,region,up=up,down=0,freq=flanking.freq,verbose=TRUE,dist=flanking.dist,s.width=s.width)
		if(no.strand==TRUE){
		strand(region)<-"-"
	}
	else if(no.strand==FALSE){
		strand(region)=="+" -> ind1
		strand(region)=="-" -> ind2
		strand(region[ind1])<- "-"
		strand(region[ind2])<- "+"
	}
	fs1.right=featureScores(x1,region,up=down,down=0,freq=flanking.freq,verbose=TRUE,dist=flanking.dist,s.width=s.width)

cbind(tables(fs1.left)[[1]],tables(fs1.middle)[[1]],t(apply(tables(fs1.right)[[1]],1,rev)))->tmp1
tmp1->out1

if(is.null(max.cvg)){
max.cvg <- quantile(c(unlist(out1)), cap.q, na.rm = TRUE)
}
else{
max.cvg<-max.cvg
}
print(max.cvg)

    f.tmp<-function(x)
    {
        x[x > max.cvg] = max.cvg
        x / max.cvg
    }
    f.tmp(out1)->tbls1
	
if(append==FALSE){
	if(frame==TRUE){
	plot(colMeans(tbls1),ylim=ylim,ylab="Normalized read density",xlab="relative position",xaxt="n",type="l",main=main,col=col)
	x.lab<-names(colMeans(tbls1))
	pos1<-1
	which(x.lab=="0 %")->pos2
	which(x.lab=="100 %")->pos3
	pos4<-length(x.lab)
	axis(1,c(pos1,pos2,pos3,pos4),x.lab[c(pos1,pos2,pos3,pos4)])
	}else{
	plot(colMeans(tbls1),ylim=ylim,ylab="",xlab="",xaxt="n",yaxt="n",type="l",main=main,col=col)
	}
}else{
	points(colMeans(tbls1),type="l",main=main,col=col)	
}
	if(!is.null(legend)){
		legend("topright", legend = legend,  text.col = text.col)
	}

}


hist.plot<-function(x,breaks=10,ylim=c(0,0.3),append=FALSE,col="black",ylab="frequency",xlab="",max.cvg=NULL,min.cvg=NULL,cap.q=0.99,spar=NULL,type="l",xlim=NULL){
	na.omit(x)->x

	if(is.null(max.cvg)){
	quantile(x,cap.q)->max1
	quantile(x,1-cap.q)->max2
	}
	else{
		max1<-max.cvg
		max2<-min.cvg
	}
	x[x>max1]<-max1
	x[x<max2]<-max2

	cat(range(x))
	hist(x,breaks=breaks,plot=F)->out
	#print(out$mids)
	#print(out$counts/sum(out$counts))
	if(append==FALSE){
		plot(smooth.spline(out$mids,out$counts/sum(out$counts),spar=spar),type=type,ylim=ylim,ylab=ylab,xlab=xlab, xlim=xlim)
		points(out$mids,out$counts/sum(out$counts),type="p")
	}
	else{
		points(smooth.spline(out$mids,out$counts/sum(out$counts),spar=spar),type=type,col=col)
		points(out$mids,out$counts/sum(out$counts),type="p",col=col)
	}
}

density.plot<-function(x,y,cap.q=0.99,max.cvg=NULL,min.cvg=NULL,single.tail=FALSE,sample1=NULL,sample2=NULL,xlab=NULL,append=FALSE,lty=c(1,1),col=c("black","red"),text.col=c("black","red"),main=NULL){
	na.omit(x)->x
	na.omit(y)->y
	require(sm)
	c(x,y)->data
	factor(c(rep(1,length(x)),rep(2,length(y))),labels=c(sample1,sample2))->fa

if(single.tail==FALSE){
	if(is.null(max.cvg)){
	quantile(data,cap.q)->max1
	quantile(data,1-cap.q)->max2
	}
	else{
		max1<-max.cvg
		max2<-min.cvg
	}
	data[data>max1]<-max1
	data[data<max2]<-max2
}
else{
	if(is.null(max.cvg)){
	quantile(data,cap.q)->max1
	}
	else{
		max1<-max.cvg
	}
	data[data>max1]<-max1
}
	cat(range(data))

	if(append==FALSE){
		sm.density.compare(data, fa, xlab=xlab,col=col,lty=lty,xlim=c(max1,max2))
	
	}
	else{
	}
	legend=c(sample1,sample2)
	legend("topright",legend=legend,text.col=text.col)
	title(main)
}


#####################################

bedtools.coveragebed<-function(functionstring="coverageBed",bed1,bed2,opt.string="")
{
  #create temp files
  #a.file=tempfile()
  #b.file=tempfile()
  out   =tempfile()
  options(scipen =99) # not to use scientific notation when writing out
 
  #write.table(bed1,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
  #write.table(bed2,file=b.file,quote=F,sep="\t",col.names=F,row.names=F)
 
  # create the command string and call the command using system()
  command=paste(functionstring,"-a",bed1,"-b",bed2,opt.string,">",out,sep=" ")
  cat(command,"\n")
  try(system(command))
 
  res=read.table(out,header=F)
  #unlink(a.file);
  #unlink(b.file);
  unlink(out)
  return(res)
}

#####################################

bedtools.coveragebed.range<-function(functionstring="coverageBed",bed1,bed2,opt.string="",mcol=TRUE)
{
  #create temp files
  a.file=tempfile()
  b.file=tempfile()
  out   =tempfile()
  options(scipen =99) # not to use scientific notation when writing out
 data.frame(chr=as.character(seqnames(bed1)),start=start(bed1),end=end(bed1))->a.out
 if(mcol==TRUE){
 data.frame(chr=as.character(seqnames(bed2)),start=start(bed2),end=end(bed2),gene=mcols(bed2)$gene)->b.out
 }
 else{
 data.frame(chr=as.character(seqnames(bed2)),start=start(bed2),end=end(bed2))->b.out
 }

write.table(a.out,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
write.table(b.out,file=b.file,quote=F,sep="\t",col.names=F,row.names=F)
 
  # create the command string and call the command using system()
  command=paste(functionstring,"-a",a.file,"-b",b.file,opt.string,">",out,sep=" ")
  cat(command,"\n")
  try(system(command))
 
  res=read.table(out,header=F)
  return(res)
}

