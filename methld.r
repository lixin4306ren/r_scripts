methblock<-function(meth,pos,num=3,file="meth_block.bed",log=FALSE){

chr<-unique(pos$chr)
for(i in 1:length(chr)){
pos[pos$chr==chr[i],]->meth.info
sort(meth.info$pos,index.return=T)$ix->ind
meth.info[ind,]->meth.info
#which(rownames(meth) %in% rownames(meth.info)) ->ind
match(rownames(meth.info),rownames(meth)) ->ind
meth[ind,]->chr.meth

if(log==TRUE){
ilogit=function(x) 1/(1+exp(-x))
ilogit(chr.meth)->chr.meth
cat ("log TRUE")
}

meth.info$position->tmp

methblock2(chr.meth,chr=chr[i])->block
#head(block,100)
#if(1)stop("exit")
block$cgstart<-block$start
block$cgend<-block$end
tmp[block$start]->block$start
tmp[block$end]->block$end
block$len<-block$end-block$start+1
block[block$num>=num,]->block

write.table(file=file,x=block,quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)
}
block
} 

######################################

methblock.whole<-function(meth,pos,num=3,file="meth_block.rda",log=FALSE){

chr<-unique(pos$chr)

  pos[order(pos$chr,pos$position),]->pos
  match(rownames(meth),rownames(pos)) ->ind
  meth[ind,]->meth
  
  if(log==TRUE){
  ilogit=function(x) 1/(1+exp(-x))
  ilogit(meth)->meth
  cat ("log TRUE")
  }
  
  pos$position->tmp.pos
  pos$chr->tmp.chr

methblock2(meth,chr="genome")->block
block$cgstart<-block$start
block$cgend<-block$end
tmp.pos[block$start]->block$start
tmp.pos[block$end]->block$end
tmp.chr[block$cgstart]->block$chr
block$len<-block$end-block$start+1
block[block$num>=num,]->block
save(x=block,file=file)
} 

left.com<-function(a,b){
 b1 <- matrix(apply(b,1,function(x) paste(b[1,],b[2,],sep=" ")),,2)
 a1 <- matrix(apply(a,1,function(x) paste(a[1,],a[2,],sep=" ")),,2)
 m<-b1[,1] %in% a1[,1]
 b[,!m]
}



############### basic methblock finder
methblock2<-function(data,a=0.4,b=0.5,chr=chr){
n=nrow(data);
start<-vector();
end<-vector();
pair<-vector();
index<-1;
tag<-0;
total.pair=0;
n.sig=0;
i=1

while (i<=n){
	start[index]<-i
        for(j in (i+1):n){
	   total.pair=0
	   
	   n.sig=0
          if(j-i>1){
		left<-left.com(combn(i:(j-1),2),combn(i:j,2));
		for(k in 1:length(left[1,])) {
			total.pair<-total.pair+1
			if(sum((data[left[,k][1],] & data[left[,k][2],]),na.rm=T)==0){next}
			if(sum(!(is.na(data[left[,k][1],])) & !(is.na(data[left[,k][2],])))<5){total.pair<-total.pair-1;next;}
			r<-cor(data[left[,k][1],],data[left[,k][2],],use="complete.obs")^2;
			if(is.na(r)){next}else{
			if (r>=a){n.sig<-n.sig+1}}
		}
	  }
	  else{
		com=combn(i:j,2);
		for(k in 1:length(com[1,])) {
			total.pair<-total.pair+1
			if(sum((data[com[,k][1],] & data[com[,k][2],]),na.rm=T)==0){next}
			if(sum(!(is.na(data[com[,k][1],])) & !(is.na(data[com[,k][2],])))<5){next;}
			r<-cor(data[com[,k][1],],data[com[,k][2],],use="complete.obs")^2;
			if(is.na(r)){next}else{
			if (r>=a){n.sig<-n.sig+1}}
		}
	  }



          if(total.pair > 0 & n.sig/total.pair>=b){
                end[index]<-j;tag<-0;
          }
	  else if(total.pair==0){end[index]<-j;tag<-0}
          else{
                end[index]<-j-1;
		pair[index]<-total.pair;
                index<-index+1;
                total.pair=0;
                n.sig=0;
                tag<-1;
                i<-j-1;
                break;
          }
        }
        i<-i+1;
	
        if(j==n & tag==0){pair[index]<-total.pair;break;}
        else if(j==n & tag==1){start[index]<-j;end[index]<-j;pair[index]<-total.pair;break;}
}
data.frame(chr=paste("chr",chr,sep=""),start=start,end=end,num=end-start+1,pair=pair)->block
block[block$end-block$start>0,]

}

methblock3<-function(data,chr=chr){
n=nrow(data)
block<-NULL
i<-1
while (i<=n-1){
	if(sum((data[i,] & data[i+1,]),na.rm=T)==0){next}
	r<-cor(data[i,],data[i+1,],use="complete.obs")^2;
        
        i<-i+1;
data.frame(chr=paste("chr",chr,sep=""),start=i,end=i+1,r=r)->tmp
rbind(block,tmp)->block

}
return(block)
}






clear.single<-function(g,meth.info,dist=100){
out<-vector()
n=length(g)
for (i in 1:n){
   if(length(g[[i]])==1){
	out[i]<-FALSE
   }
   else{
	nn<-length(g[[i]])
	meth.info$position[match(names(g[[i]]),rownames(meth.info))] -> tmp.pos
	for (j in 1:nn){
	 if(sum(abs(tmp.pos-tmp.pos[j]) < dist)<2){
		g[[i]][j]<-NA
	 }
	}
	na.omit(g[[i]])->g[[i]]
	if(length(g[[i]])>=2){
		out[i]<-TRUE
	}
	else{out[i]<-FALSE}
   }
}
g[out]->g
return (g)
}

meth.pattern.dist<-function(block1,block2){
library(GenomicRanges)
GRanges(seqnames=block1[,1],ranges=IRanges(start=block1[,7],end=block1[,8]))->block1
block1[seqnames(block1)!="chrX" & seqnames(block1)!="chrY"] ->block1
GRanges(seqnames=block2[,1],ranges=IRanges(start=block2[,7],end=block2[,8]))->block2
block2[seqnames(block2)!="chrX" & seqnames(block2)!="chrY"] ->block2

dist<-sum(width(intersect(block1,block2)))/sum(width(union(block1,block2)))
dist
}

block.interact.intra<-function(block,meth,chr,a=0.4,b=0.5){
  block<-block[block[,1]==paste("chr",chr,sep=""),]
  ilogit=function(x) 1/(1+exp(-x))
  ilogit(meth)->meth
  pos[pos$chr==chr,]->meth.info
  sort(meth.info$pos,index.return=T)$ix->ind
  meth.info[ind,]->meth.info
  match(rownames(meth.info),rownames(meth)) ->ind
  meth[ind,]->chr.meth

n<-nrow(block)
com<-combn(n,2)
out<-data.frame()
index=1
for(i in 1:dim(com)[2]){
	start1=block[com[1,i],7]
	end1=block[com[1,i],8]
	start2=block[com[2,i],7]
	end2=block[com[2,i],8]
	num.sig=0
	total=0

	for(j in start1:end1){
		for(k in start2:end2){
			abs(cor(chr.meth[j,],chr.meth[k,]))^2->r
			total<-total+1
			#print(r)
			if (r>=a) {num.sig<-num.sig+1}
		}
	}
	if(num.sig/total>=b){
		out[index,1]<-chr
		out[index,2]<-com[1,i]
		out[index,3]<-com[2,i]
		out[index,4]<-start1
		out[index,5]<-end1
		out[index,6]<-start2
		out[index,7]<-end2
		out[index,8]<-block[com[1,i],4]
		out[index,9]<-block[com[2,i],4]
		out[index,10]<-block[com[1,i],2]
		out[index,11]<-block[com[1,i],3]
		out[index,12]<-block[com[2,i],2]
		out[index,13]<-block[com[2,i],3]
		out[index,14]<-min(abs(start1-end2),abs(start2-end1))-1
		index<-index+1
	}
	
}

return (out)
}


block.interact.inter<-function(block,meth,chr,a=0.4,b=0.5){
  ilogit=function(x) 1/(1+exp(-x))
  ilogit(meth)->meth

  block1<-block[block[,1]==paste("chr",chr,sep=""),]


  pos[pos$chr==chr,]->meth.info
  sort(meth.info$pos,index.return=T)$ix->ind
  meth.info[ind,]->meth.info
  match(rownames(meth.info),rownames(meth)) ->ind
  meth[ind,]->chr1.meth

  pos[pos$chr!=chr,]->meth.info

chrn<-unique(meth.info$chr)
out<-data.frame()
index=1


for(i.chr in 1:length(chrn)){
  pos[pos$chr==chrn[i.chr],]->meth.info
  sort(meth.info$pos,index.return=T)$ix->ind
  meth.info[ind,]->meth.info
  match(rownames(meth.info),rownames(meth)) ->ind
  meth[ind,]->chr2.meth
  block2<-block[block[,1]==paste("chr",chrn[i.chr],sep=""),]


for(i in 1:nrow(block1)){
	for (ii in 1:nrow(block2)){
	start1=block1[i,7]
	end1=block1[i,8]
	start2=block2[ii,7]
	end2=block2[ii,8]
	num.sig=0
	total=0
	
	for(j in start1:end1){
		for(k in start2:end2){
			abs(cor(chr1.meth[j,],chr2.meth[k,]))^2->r
			total<-total+1
			if (r>=a) {num.sig<-num.sig+1}
		}
	}
	if(num.sig/total>=b){
		out[index,1]<-chr
		out[index,2]<-start1
		out[index,3]<-end1
		out[index,4]<-block1[i,4]
		out[index,5]<-block1[i,2]
		out[index,6]<-block1[i,3]
		out[index,7]<-chrn[i.chr]
		out[index,8]<-start2
		out[index,9]<-end2
		out[index,10]<-block2[ii,4]
		out[index,11]<-block2[ii,2]
		out[index,12]<-block2[ii,3]
		index<-index+1
	}

	
	}
}


}
retur(out)
}


block.interact.pair<-function(block1,block2,meth,a=0.4,b=0.5){
out<-data.frame()
index=1

	start1=block[1,6]
	end1=block[1,7]
	start2=block[1,6]
	end2=block[1,7]
	num.sig=0
	total=0

	for(j in start1:end1){
		for(k in start2:end2){
			abs(cor(meth[j,],meth[k,]))^2->r
			total<-total+1
			if (r>=a) {num.sig<-num.sig+1}
		}
	}
	if(num.sig/total>=b){
		out[index,1]<-com[1,i]
		out[index,2]<-com[2,i]
		out[index,3]<-start1
		out[index,4]<-end1
		out[index,5]<-start2
		out[index,6]<-end2
		out[index,7]<-block[com[1,i],4]
		out[index,8]<-block[com[2,i],4]
		out[index,9]<-block[com[1,i],2]
		out[index,10]<-block[com[1,i],3]
		out[index,11]<-block[com[2,i],2]
		out[index,12]<-block[com[2,i],3]

		index<-index+1
	}
	
out
}

bootstrap<-function(meth,pos,num=3,chr,n=100,cutoff=0.5,log=FALSE){
library(GenomicRanges)
as.matrix(meth)->meth

pos[pos$chr==chr,]->meth.info
sort(meth.info$pos,index.return=T)$ix->ind
meth.info[ind,]->meth.info
#which(rownames(meth) %in% rownames(meth.info)) ->ind
match(rownames(meth.info),rownames(meth)) ->ind
meth[ind,]->chr.meth

if(log==TRUE){
ilogit=function(x) 1/(1+exp(-x))
ilogit(chr.meth)->chr.meth
cat ("log TRUE")
}
meth.info$position->tmp
methblock2(chr.meth,chr=chr)->ref.block
ref.block[ref.block$num>=num,]->ref.block

if(nrow(ref.block)==0){
return(NULL)
}
else{
#print(ref.block)
ref.range<-GRanges(seqnames=chr,ranges=IRanges(start=ref.block$start,end=ref.block$end))
ind<-vector(length=length(ref.range))
if (cutoff>0 & n>0){
for (i in 1:n){
#print(i)
sample.name<-colnames(chr.meth)
chr.meth[,sample(sample.name,length(sample.name),replace=T)]->tmp.meth
methblock2(tmp.meth,chr=chr)->tmp.block
if(nrow(tmp.block)==0){
	tmp.range<-GRanges(seqnames="0",ranges=IRanges(start=0,end=0))
}else{
	tmp.range<-GRanges(seqnames=chr,ranges=IRanges(start=tmp.block$start,end=tmp.block$end))
	ind<-ind+(ref.range %in% tmp.range)
}
}

ref.block$bootstrap<-ind/n

which(ind>=cutoff*n)->ind
ref.block[ind,]->ref.block
}
else{
ref.block$bootstrap<-0
}

ref.block$cgstart<-ref.block$start
ref.block$cgend<-ref.block$end
tmp[ref.block$start]->ref.block$start
tmp[ref.block$end]->ref.block$end
ref.block$len<-ref.block$end-ref.block$start+1

return(ref.block)
#print(ref.block)
}
}



######################
###################### bootstrap only for VMRs
bootstrap.vmr<-function(vmr.list,meth,pos,num=3,n=1000,log=FALSE){
block<-data.frame()
for(i in 1:nrow(vmr.list)){
print(i)
sub("chr","",vmr.list[i,1])->chr
pos$chr==chr -> ind
pos[ind,]->tmp.pos
rownames(tmp.pos[sort(which(tmp.pos$position>=vmr.list[i,2] & tmp.pos$position<=vmr.list[i,3])),]) -> tmp.vmr.cg.name
beta[tmp.vmr.cg.name,]->tmp.vmr.meth
tmp.pos[sort(which(tmp.pos$position>=vmr.list[i,2] & tmp.pos$position<=vmr.list[i,3])),]->tmp.vmr.pos
bootstrap(meth=tmp.vmr.meth,pos=tmp.vmr.pos,chr=chr,n=n,log=log)->tmp.block
rbind(block,tmp.block)->block
}
return(block)
}


###################### for across tissue
bootstrap2<-function(meth,pos,tissue,meanmeth,num=3,chr,n=1000,cutoff=0.95,log=FALSE){
library(GenomicRanges)

order(pos$chr,pos$position)->ind
pos[ind,]->pos
meth[ind,]->meth
meanmeth[ind,]->meanmeth

  pos[pos$chr==chr,]->meth.info
  match(rownames(meth.info),rownames(meth)) ->ind
  meth[ind,]->meth
  meth->chr.meth
  meanmeth[ind,]->ref.meth

  
  
  if(log==TRUE){
  ilogit=function(x) 1/(1+exp(-x))
  ilogit(chr.meth)->chr.meth
  ilogit(ref.meth)->ref.meth
  cat ("log TRUE")
  }
meth.info$position->tmp.pos
methblock2(ref.meth,chr=chr)->ref.block
ref.block[ref.block$num>=num,]->ref.block
ref.range<-GRanges(seqnames=chr,ranges=IRanges(start=ref.block$start,end=ref.block$end))

unique(tissue)->tissue.vector

ind.boot<-vector(length=length(ref.range))
if (cutoff>0){

for (i in 1:n){
out<-NULL
print (i)
for (j in 1:length(tissue.vector)) {
sample(which(tissue==tissue.vector[j]),1,replace=T)->ind
meth[,ind]->tmp
matrix(tmp,length(tmp),1,dimnames=list(rownames(tmp),tissue.vector[j]))->tmp
cbind(out,tmp)->out
}
methblock2(out,chr=chr)->tmp.block
tmp.range<-GRanges(seqnames=chr,ranges=IRanges(start=tmp.block$start,end=tmp.block$end))
ind.boot<-ind.boot+(ref.range %in% tmp.range) 
#print (ind)
}

ref.block$bootstrap<-ind.boot/n
which(ind.boot>=cutoff*n)->ind.boot
ref.block[ind.boot,]->ref.block
}
else{
ref.block$bootstrap<-0
}

ref.block$cgstart<-ref.block$start
ref.block$cgend<-ref.block$end
tmp.pos[ref.block$start]->ref.block$start
tmp.pos[ref.block$end]->ref.block$end
ref.block$len<-ref.block$end-ref.block$start+1

ref.block
}


bootstrap.seq<-function(meth,pos,num=3,chr,n=1000,cutoff=0.95){
library(GenomicRanges)
load(meth)
load(pos)
get("meth.tmp")->chr.meth
get("pos.tmp")->meth.info



meth.info$position->tmp
methblock2(chr.meth,chr=chr)->ref.block
ref.block[ref.block$num>=num,]->ref.block
ref.range<-GRanges(seqnames=chr,ranges=IRanges(start=ref.block$start,end=ref.block$end))


ind<-vector(length=length(ref.range))
if (cutoff>0){
for (i in 1:n){
print (i)
sample.name<-colnames(chr.meth)
chr.meth[,sample(sample.name,length(sample.name),replace=T)]->tmp.meth
methblock2(tmp.meth,chr=chr)->tmp.block
tmp.range<-GRanges(seqnames=chr,ranges=IRanges(start=tmp.block$start,end=tmp.block$end))
ind<-ind+(ref.range %in% tmp.range) 
}

ref.block$bootstrap<-ind/n

which(ind>=cutoff*n)->ind
ref.block[ind,]->ref.block
}
else{
ref.block$bootstrap<-0
}

ref.block$cgstart<-ref.block$start
ref.block$cgend<-ref.block$end
tmp[ref.block$start]->ref.block$start
tmp[ref.block$end]->ref.block$end
ref.block$len<-ref.block$end-ref.block$start+1

ref.block
}

###############################
############################### plot
plot.meth.block<-function(block=NULL,meth=NULL,pos=NULL){
 pdf("block.heatmap.pdf")
 nrow(block)->n
#cat(n,"\n")
 i<-1
#cat(i,"\n")
 for(i in 1:n){
start<-block[i,6]-50
 if(start<1){start=1}
 end<-block[i,7]+50
 if(end>max(block[,7])){end=max(block[,7])}
 chr<-sub("chr","",block[i,1])
 pos[pos$chr==chr,]->meth.info
 match(rownames(meth.info),rownames(meth)) ->ind
 meth[ind,]->chr.meth
 cor(t(chr.meth[start:end,]),method="pearson")->r
 r^2->r
 r[lower.tri(r)]<-NA
 library(gplots)
 ilogit=function(x) 1/(1+exp(-x))
tmp1<-block[i,6]-start+1;tmp2<-block[i,7]-start+1
heatmap.2(r,Rowv=FALSE,Colv=FALSE,dendrogram="none",breaks=c(0,0.1,0.5,1),col=rev(heat.colors(3)),trace="none",add.expr={lines(c(tmp1,tmp1),c(0,tmp2));lines(c(tmp2,tmp2),c(0,tmp1));points(1:(end-start+1),rowMeans(ilogit(chr.meth[start:end,]))*10,type="l",col="red")})
 }
 dev.off()
}


chr.meth.r.plot<-function(meth,pos,chr,file=file,log=FALSE){


  pos[pos$chr==chr,]->meth.info
  sort(meth.info$pos,index.return=T)$ix->ind
  meth.info[ind,]->meth.info
  #which(rownames(meth) %in% rownames(meth.info)) ->ind
  match(rownames(meth.info),rownames(meth)) ->ind
  meth[ind,]->chr.meth
  
  if(log==TRUE){
  ilogit=function(x) 1/(1+exp(-x))
  ilogit(chr.meth)->chr.meth
  cat ("log TRUE")
  }
  


methblock3(chr.meth,chr=chr)->block
meth.info$position->tmp

block$cgstart<-block$start
block$cgend<-block$end
tmp[block$start]->block$start
tmp[block$end]->block$end
write.table(file=file,x=block,quote=FALSE,row.names=FALSE,col.names=FALSE)
return (block)
}




##################### split meth and pos file for large dataset
#####################
split.meth.pos<-function(file,len=20000,split=TRUE){
load(file)
unlist(strsplit(file,"\\."))->str
chr=str[1]
tissue=str[2]
get(paste(chr,tissue,sep="."))->x
x[,4]==1 ->ind
x[,4]==4 ->ind2

x[ind,]->cg.plus
x[ind2,]->cg.minus
if(!identical(as.integer(cg.minus[,2]-1),as.integer(cg.plus[,2]))){stop("wrong")}

seq(5,dim(x)[2],3)->ind
seq(6,dim(x)[2],3)->ind2

cg.plus[,ind]+cg.minus[,ind]->m
cg.plus[,ind2]+cg.minus[,ind2]->c
c <5 ->ind
c[ind]<-0
meth<-as.matrix(m/c)
is.infinite(meth)->ind
meth[ind]<-NaN

is.na(meth)->tmp
rowSums(tmp) >= 10 ->ind

cg.plus[ind,1:2]->pos
meth[ind,]->meth
names(pos)<-c("chr","position")



seq(1,nrow(pos),by=len)->start
end<-start+len-1
end[end > nrow(pos)]<-nrow(pos)

if (split==TRUE){
for(i in 1:length(start)){
	file<-paste(chr,tissue,"meth",i,"rda",sep=".")
	meth.tmp<-meth[start[i]:end[i],]
	save(meth.tmp,file=file)
	file<-paste(chr,tissue,"pos",i,"rda",sep=".")
	pos.tmp<-pos[start[i]:end[i],]
	save(pos.tmp,file=file)
	}
}
else{
	file<-paste(chr,tissue,"meth","rda",sep=".")
	save(meth,file=file)
	file<-paste(chr,tissue,"pos","rda",sep=".")
	save(pos,file=file)
}

}



#########################
split.nonCG.meth.pos<-function(file,len=20000,split=TRUE){
load(file)
unlist(strsplit(file,"\\."))->str
chr=str[1]
tissue=str[2]
get(paste(chr,tissue,sep="."))->x
x[,4]!=1 & x[,4]!=4 ->ind

x[ind,]->noncg

seq(5,dim(x)[2],3)->ind
seq(6,dim(x)[2],3)->ind2

noncg[,ind]->m
noncg[,ind2]->c
c <5 ->ind  
c[ind]<-0
meth<-as.matrix(m/c)
is.infinite(meth)->ind
meth[ind]<-NaN

is.na(meth)->tmp
rowSums(tmp) >= 10 ->ind

noncg[ind,1:2]->pos
meth[ind,]->meth
names(pos)<-c("chr","position")



seq(1,nrow(pos),by=len)->start
end<-start+len-1
end[end > nrow(pos)]<-nrow(pos)

if (split==TRUE){
for(i in 1:length(start)){
	file<-paste(chr,tissue,"nonCG.meth",i,"rda",sep=".")
	meth.tmp<-meth[start[i]:end[i],]
	save(meth.tmp,file=file)
	file<-paste(chr,tissue,"nonCG.pos",i,"rda",sep=".")
	pos.tmp<-pos[start[i]:end[i],]
	save(pos.tmp,file=file)
	}
}
else{
	file<-paste(chr,tissue,"nonCG.meth","rda",sep=".")
	save(meth,file=file)
	file<-paste(chr,tissue,"nonCG.pos","rda",sep=".")
	save(pos,file=file)
}

}










heatmap.3<-function (x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
    distfun = dist, hclustfun = hclust, dendrogram = c("both",
        "row", "column", "none"), symm = FALSE, scale = c("none",
        "row", "column"), na.rm = TRUE, revC = identical(Colv,
        "Rowv"), add.expr, breaks, symbreaks = min(x < 0, na.rm = TRUE) ||
        scale != "none", col = "heat.colors", colsep, rowsep,
    sepcolor = "white", sepwidth = c(0.05, 0.05), cellnote, notecex = 1,
    notecol = "cyan", na.color = par("bg"), trace = c("column",
        "row", "both", "none"), tracecol = "cyan", hline = median(breaks),
    vline = median(breaks), linecol = tracecol, margins = c(5,
        5), ColSideColors, RowSideColors, cexRow = 0.2 + 1/log10(nr),
    cexCol = 0.2 + 1/log10(nc), labRow = NULL, labCol = NULL,
    key = TRUE, keysize = 1.5, density.info = c("histogram",
        "density", "none"), denscol = tracecol, symkey = min(x <
        0, na.rm = TRUE) || symbreaks, densadj = 0.25, main = NULL,
    xlab = NULL, ylab = NULL, lmat = NULL, lhei = NULL, lwid = NULL,
    ...)
{
    scale01 <- function(x, low = min(x), high = max(x)) {
        x <- (x - low)/(high - low)
        x
    }
    retval <- list()
    scale <- if (symm && missing(scale))
        "none"
    else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (length(col) == 1 && is.character(col))
        col <- get(col, mode = "function")
    if (!missing(breaks) && (scale != "none"))
        warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
    if (is.null(Rowv) || is.na(Rowv))
        Rowv <- FALSE
    if (is.null(Colv) || is.na(Colv))
        Colv <- FALSE
    else if (Colv == "Rowv" && !isTRUE(Rowv))
        Colv <- FALSE
    if (length(di <- dim(x)) != 2 || !is.numeric(x))
        stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1)
        stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2)
        stop("`margins' must be a numeric vector of length 2")
    if (missing(cellnote))
        cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
        if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
            c("both", "row"))) {
            if (is.logical(Colv) && (Colv))
                dendrogram <- "column"
            else dedrogram <- "none"
            warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting row dendogram.")
        }
    }
    if (!inherits(Colv, "dendrogram")) {
        if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
            c("both", "column"))) {
            if (is.logical(Rowv) && (Rowv))
                dendrogram <- "row"
            else dendrogram <- "none"
            warning("Discrepancy: Colv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting column dendogram.")
        }
    }
    if (inherits(Rowv, "dendrogram")) {
        ddr <- Rowv
        rowInd <- order.dendrogram(ddr)
    }
    else if (is.integer(Rowv)) {
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd))
            stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
        Rowv <- rowMeans(x, na.rm = na.rm)
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd))
            stop("row dendrogram ordering gave index of wrong length")
    }
    else {
        rowInd <- nr:1
    }
    if (inherits(Colv, "dendrogram")) {
        ddc <- Colv
        colInd <- order.dendrogram(ddc)
    }
    else if (identical(Colv, "Rowv")) {
        if (nr != nc)
            stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
        if (exists("ddr")) {
            ddc <- ddr
            colInd <- order.dendrogram(ddc)
        }
        else colInd <- rowInd
    }
    else if (is.integer(Colv)) {
        hcc <- hclustfun(distfun(if (symm)
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd))
            stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
        Colv <- colMeans(x, na.rm = na.rm)
        hcc <- hclustfun(distfun(if (symm)
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd))
            stop("column dendrogram ordering gave index of wrong length")
    }
    else {
        colInd <- 1:nc
    }
    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()
    x <- x[rowInd, colInd]
    x.unscaled <- x
    cellnote <- cellnote[rowInd, colInd]
    if (is.null(labRow))
        labRow <- if (is.null(rownames(x)))
            (1:nr)[rowInd]
        else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol))
        labCol <- if (is.null(colnames(x)))
            (1:nc)[colInd]
        else colnames(x)
    else labCol <- labCol[colInd]
    if (scale == "row") {
        retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
        x <- sweep(x, 1, rm)
        retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
        x <- sweep(x, 2, rm)
        retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) <
        1) {
        if (missing(col) || is.function(col))
            breaks <- 16
        else breaks <- length(col) + 1
    }
    if (length(breaks) == 1) {
        if (!symbreaks)
            breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                length = breaks)
        else {
            extreme <- max(abs(x), na.rm = TRUE)
            breaks <- seq(-extreme, extreme, length = breaks)
        }
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (class(col) == "function")
        col <- col(ncol)
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks
#    if (missing(lhei) || is.null(lhei))
#        lhei <- c(keysize, 4)
#    if (missing(lwid) || is.null(lwid))
#        lwid <- c(keysize, 4)
#    if (missing(lmat) || is.null(lmat)) {
#        lmat <- rbind(4:3, 2:1)
#        if (!missing(ColSideColors)) {
#            if (!is.character(ColSideColors) || length(ColSideColors) !=
#                nc)
#                stop("'ColSideColors' must be a character vector of length ncol(x)")
#            lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] +
#                1)
#            lhei <- c(lhei[1], 0.2, lhei[2])
#        }
#        if (!missing(RowSideColors)) {
#            if (!is.character(RowSideColors) || length(RowSideColors) !=
#                nr)
#                stop("'RowSideColors' must be a character vector of length nrow(x)")
#            lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) -
#                1), 1), lmat[, 2] + 1)
#            lwid <- c(lwid[1], 0.2, lwid[2])
#        }
#        lmat[is.na(lmat)] <- 0
#    }
#    if (length(lhei) != nrow(lmat))
#        stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
#    if (length(lwid) != ncol(lmat))
#        stop("lwid must have length = ncol(lmat) =", ncol(lmat))
#    op <- par(no.readonly = TRUE)
#    on.exit(par(op))
#    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
    if (!missing(RowSideColors)) {
        par(mar = c(margins[1], 0, 0, 0.5))
        image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    }
    if (!missing(ColSideColors)) {
        par(mar = c(0.5, 0, 0, margins[2]))
        image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    }
    par(mar = c(margins[1], 0, 0, margins[2]))
    x <- t(x)
    cellnote <- t(cellnote)
    if (revC) {
        iy <- nr:1
        if (exists("ddr"))
            ddr <- rev(ddr)
        x <- x[, iy]
        cellnote <- cellnote[, iy]
    }
    else iy <- 1:nr
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 +
        c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col,
        breaks = breaks, ...)
    retval$carpet <- x
    if (exists("ddr"))
        retval$rowDendrogram <- ddr
    if (exists("ddc"))
        retval$colDendrogram <- ddc
    retval$breaks <- breaks
    retval$col <- col
    if (!invalid(na.color) & any(is.na(x))) {
        mmat <- ifelse(is.na(x), 1, NA)
        image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
            col = na.color, add = TRUE)
    }
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
        cex.axis = cexCol)
    if (!is.null(xlab))
        mtext(xlab, side = 1, line = margins[1] - 1.25)
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
        cex.axis = cexRow)
    if (!is.null(ylab))
        mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr))
        eval(substitute(add.expr))
    if (!missing(colsep))
        for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0,
            length(csep)), xright = csep + 0.5 + sepwidth[1],
            ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1,
            col = sepcolor, border = sepcolor)
    if (!missing(rowsep))
        for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) +
            1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) +
            1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1,
            col = sepcolor, border = sepcolor)
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    if (trace %in% c("both", "column")) {
        retval$vline <- vline
        vline.vals <- scale01(vline, min.scale, max.scale)
        for (i in colInd) {
            if (!is.null(vline)) {
                abline(v = i - 0.5 + vline.vals, col = linecol,
                  lty = 2)
            }
            xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
            xv <- c(xv[1], xv)
            yv <- 1:length(xv) - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (trace %in% c("both", "row")) {
        retval$hline <- hline
        hline.vals <- scale01(hline, min.scale, max.scale)
        for (i in rowInd) {
            if (!is.null(hline)) {
                abline(h = i + hline, col = linecol, lty = 2)
            }
            yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
            yv <- rev(c(yv[1], yv))
            xv <- length(yv):1 - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (!missing(cellnote))
        text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
            col = notecol, cex = notecex)
    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    }
    #else plot.new()
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    }
    else #plot.new()
    if (!is.null(main))
        title(main, cex.main = 1.5 * op[["cex.main"]])
    if (key) {
        par(mar = c(5, 4, 2, 1), cex = 0.75)
        tmpbreaks <- breaks
        if (symkey) {
            max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
            min.raw <- -max.raw
            tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
            tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
        }
        else {
            min.raw <- min(x, na.rm = TRUE)
            max.raw <- max(x, na.rm = TRUE)
        }
        z <- seq(min.raw, max.raw, length = length(col))
        image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
            xaxt = "n", yaxt = "n")
        par(usr = c(0, 1, 0, 1))
        lv <- pretty(breaks)
        xv <- scale01(as.numeric(lv), min.raw, max.raw)
        axis(1, at = xv, labels = lv)
        if (scale == "row")
            mtext(side = 1, "Row Z-Score", line = 2)
        else if (scale == "column")
            mtext(side = 1, "Column Z-Score", line = 2)
        else mtext(side = 1, "Value", line = 2)
        if (density.info == "density") {
            dens <- density(x, adjust = densadj, na.rm = TRUE)
            omit <- dens$x < min(breaks) | dens$x > max(breaks)
            dens$x <- dens$x[-omit]
            dens$y <- dens$y[-omit]
            dens$x <- scale01(dens$x, min.raw, max.raw)
            lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
                lwd = 1)
            axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
            title("Color Key\nand Density Plot")
            par(cex = 0.5)
            mtext(side = 2, "Density", line = 2)
        }
        else if (density.info == "histogram") {
            h <- hist(x, plot = FALSE, breaks = breaks)
            hx <- scale01(breaks, min.raw, max.raw)
            hy <- c(h$counts, h$counts[length(h$counts)])
            lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
                col = denscol)
            axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
            title("Color Key\nand Histogram")
            par(cex = 0.5)
            mtext(side = 2, "Count", line = 2)
        }
        else title("Color Key")
    }
    #else plot.new()
    retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
        high = retval$breaks[-1], color = retval$col)
    invisible(retval)
}
