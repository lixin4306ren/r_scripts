merge.count<-function(file,condition.sep="_",condition.col=1,col.name=NULL,tissue=NULL,has.condtion.col=1){
if(has.condtion.col>1){
as.vector(read.table(file)[,has.condtion.col])->group
}
as.vector(read.table(file)$V1)->file
if(!is.null(tissue)){file[grepl(pattern=tissue,file)]->file}
out<-NULL
name_col<-NULL
condition<-NULL
read.table(file[1],row.names=1)->tmp

grepl("ENS",rownames(tmp))->ind
rownames(tmp)[ind]->name
as.data.frame(tmp[ind,])->tmp
rownames(tmp)<-name
out<-tmp

#name_col<-unlist(strsplit(file[1],"\\."))[1]
name_col<-unlist(file[1])

if(has.condtion.col==1){
condition<-unlist(strsplit(file[1],condition.sep))[condition.col]
}
else{
condition<-group[1]
}
for (i in 2:length(file)){
read.table(file[i],row.names=1)->tmp

rownames(tmp)[ind]->name
if(!identical(rownames(out),name)){stop("wrong")}
grepl("ENS",rownames(tmp))->ind
as.data.frame(tmp[ind,])->tmp
rownames(tmp)<-name
out<-cbind(out,tmp)
#tmp<-unlist(strsplit(file[i],"\\."))[1]
tmp<-unlist(file[i])
append(name_col,tmp)->name_col
if(has.condtion.col==1){
tmp<-unlist(strsplit(file[i],condition.sep))[condition.col]
}
else{
tmp<-group[i]
}
append(condition,tmp)->condition

}

if(is.null(col.name)){
colnames(out)<-name_col
}else{
colnames(out)<-col.name
}

library(DESeq)
cds = newCountDataSet(out, condition )
return(cds)
}


merge.count2<-function(file){
as.vector(read.table(file)$V1)->file
out<-NULL
name_col<-NULL
condition<-NULL
read.table(file[1],row.names=1)->tmp

grepl("ENSG",rownames(tmp))->ind
rownames(tmp)[ind]->name
as.data.frame(tmp[ind,])->tmp
rownames(tmp)<-name
out<-tmp

name_col<-unlist(strsplit(file[1],"\\."))[1]
condition<-unlist(strsplit(file[1],"_"))[1]

for (i in 2:length(file)){
read.table(file[i],row.names=1)->tmp

rownames(tmp)[ind]->name
if(!identical(rownames(out),name)){stop("wrong")}
grepl("ENSG",rownames(tmp))->ind
as.data.frame(tmp[ind,])->tmp
rownames(tmp)<-name
out<-cbind(out,tmp)
tmp<-unlist(strsplit(file[i],"\\."))[1]
append(name_col,tmp)->name_col

tmp<-unlist(strsplit(file[i],"_"))[1]
append(condition,tmp)->condition

}
colnames(out)<-name_col

return(out)
}


exp_cluster_plot<-function(sample.list=NULL,method="ward"){
merge.count2(sample.list)->exp
colnames(exp)->name
source("~/r_scripts/functions.R")
rowSums(exp)>10->ind
exp[ind,]->exp
hclust(as.dist(1-cor(exp,method="spearman")),method=method)->h

pdf(file="whole_exp_cluster.pdf")
myplclust(h,lab=name,lab.col=as.integer(as.factor(name)),main="",xlab="")
dev.off()
}

cal.rpkm<-function(file,gene="/home/jhmi/xinli/amber3/no_back_up/data/gene/gene.len.for.rpkm"){
as.vector(read.table(file)$V1)->file
out<-NULL
name_col<-NULL
condition<-NULL
read.table(file[1],row.names=1)->tmp

grepl("ENSG",rownames(tmp))->ind
rownames(tmp)[ind]->name
as.data.frame(tmp[ind,])->tmp
rownames(tmp)<-name
out<-tmp

name_col<-unlist(strsplit(file[1],"\\."))[1]
condition<-unlist(strsplit(file[1],"_"))[1]

for (i in 2:length(file)){
read.table(file[i],row.names=1)->tmp

rownames(tmp)[ind]->name
if(!identical(rownames(out),name)){stop("wrong")}
grepl("ENSG",rownames(tmp))->ind
as.data.frame(tmp[ind,])->tmp
rownames(tmp)<-name
out<-cbind(out,tmp)
tmp<-unlist(strsplit(file[i],"\\."))[1]
append(name_col,tmp)->name_col

tmp<-unlist(strsplit(file[i],"_"))[1]
append(condition,tmp)->condition

}
colnames(out)<-name_col
read.table(gene)->gene.len
match(rownames(out),gene.len[,1])->ind
gene.len[ind,]->gene.len2
out2<-out/(gene.len2[,2]/1000)
out3<-t(t(out2)/colSums(out))*1000000
return(out3)
}

plot.exp.distribution<-function(rpkm,gene.list,append=FALSE,col="black"){
(names(rpkm) %in% gene.list)->ind
rpkm[ind]->rpkm
if(append==FALSE){
plot(density(log10(rpkm+0.001)),col=col)
}
else{
points(density(log10(rpkm+0.001)),col=col,type="l")
}
}

count.read.exon<-function(bam,exon,type="single"){
require(GenomicRanges)
if(type=="single"){
readGappedAlignments(bam)->x
}
else if(type=="pair"){
readGappedAlignmentPairs(bam)->x
}
countOverlaps(exon,x)->count
data.frame(names(exon),count)->out
return(out)
}



