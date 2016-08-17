cal.meth.var<-function(sample.list="sample.list",chr.list="chr.list",tissue="BM",tissue2=" ",line="H"){
read.table("chr.list")->chr.list
as.vector(chr.list[,1])->chr.list
read.table(sample.list)->s
grepl(line,s[,2]) & ( grepl(tissue,s[,3]) | grepl(tissue2,s[,3] )) ->ind
as.vector(s[,1][ind])->dir


for (j in 1:length(chr.list)){
paste(paste(dir,"all",chr.list[j],sep="/"),"cg",sep=".")->files

out<-NULL
type="sd";
read.table(files[1])->x
x<-data.frame(chr=rep(chr.list[j],length(x[,1])/2),pos=x[x[,2]==1,1],m=x[x[,2]==1,3]+x[x[,2]==4,3],c=x[x[,2]==1,4]+x[x[,2]==4,4])
out<-data.frame(chr=x[,1],pos=x[,2],meth=x[,3]/x[,4])
x[,4]<5 ->ind
out[ind,3]<-NA
#which(out[,3]<0.5) ->ind
#out[ind,3]<-0
#which(out[,3]>=0.5) ->ind
#out[ind,3]<-1

for (i in 2:length(files)) {
	read.table(files[i])->x
	data.frame(pos=x[x[,2]==1,1],m=x[x[,2]==1,3]+x[x[,2]==4,3],c=x[x[,2]==1,4]+x[x[,2]==4,4])->x
	tmp<-data.frame(x[,2]/x[,3])
	x[,3]<5 ->ind
	tmp[ind,1]<-NA
	#which(tmp[,1]<0.5) ->ind
	#tmp[ind,1]<-0
	#which(tmp[,1]>=0.5) ->ind
	#tmp[ind,1]<-1

	out<-cbind(out,tmp)
}
na.omit(out)->out

########### write ouput file for entropy calculation
tmp.line1="#1.2";
tmp.line2=paste(dim(out)[1],dim(out)[2]-2,sep="\t")
tmp.name<-paste(chr.list[j],tissue,line,"site.input2",sep=".")
write.table(paste(tmp.line1,tmp.line2,sep="\n"),file=tmp.name,col.names=F,row.names=F,quote=F)
write.table(out,file=tmp.name,col.names=T,row.names=F,quote=F,sep="\t",append=T)

########## calculate entropy
input.name=paste(chr.list[j],tissue,line,"site.input2",sep=".")
outdir<-paste(chr.list[j],tissue,line,sep=".")
dir.create(outdir)

command=paste("java -Xmx5g -jar ~/soft/QDMR/QDMR.jar infile=",input.name,",outfolder=",outdir,",SD=0.07",sep="")
try(system(command))
}
}

####################
####################
####################

plot.meth.var<-function(data,win=5000,outfilename="meth_var.pdf"){
load(data)
chr.list=unique(out.mean.sd.en[,1])
pdf(width=15,height=10,file=outfilename)
layout(matrix(c(1,2),2,1),2,1)
for (j in 1:length(chr.list)){
print(chr.list[j])
tmp<-out.mean.sd.en[out.mean.sd.en[,1]==chr.list[j],]

plot(tmp[,2],tmp$sd,xlab=chr.list[j],ylab="SD",pch=".",ylim=c(0,0.2))
points(smooth.spline(tmp[,2],tmp$sd,spar=0.25),type="l")

plot(tmp[,2],tmp$en,xlab=chr.list[j],ylab="Entropy",pch=".")
points(smooth.spline(tmp[,2],tmp$en,spar=0.25),type="l")

}
dev.off()
}


#########################
#########################
#########################

cal.sliding.meth.mean.sd.en<-function(sample.list="sample.list",chr.list="chr.list",tissue="BM",line="H",win=5000){
read.table(chr.list)->chr.list
as.vector(chr.list[,1])->chr.list
source("~/r_scripts/functions.R")
out.mean.sd.en<-NULL
for (j in 1:length(chr.list)){
print(chr.list[j])
input.name=paste(chr.list[j],tissue,line,"site.input2",sep=".")

############ sd method
read.table(input.name,skip=3)->out
apply(out[,3:dim(out)[2]],1,sd)->sd
apply(out[,3:dim(out)[2]],1,mean)->mean
data.frame(out[,1],out[,2],sd)->out
mean.sliding(out,win.size=win)->tmp.sd
data.frame(out[,1],out[,2],mean)->out
mean.sliding(out,win.size=win)->tmp.mean


############ entropy meth
path=paste(paste(chr.list[j],tissue,line,sep="."),"/Entropy.txt",sep="")
read.table(path,header=T)->en
mean.sliding(en,win.size=win)->tmp.en

tmp.out<-data.frame(chr=rep(as.character(out[1,1]),nrow(tmp.en)),pos=tmp.en$id*win-win/2,mean=tmp.mean$fst,sd=tmp.sd$fst,en=tmp.en$fst)
out.mean.sd.en<-rbind(out.mean.sd.en,tmp.out)

}
file.name<-paste("all",line,tissue,win,"mean.sd.en.rda",sep=".")
save(out.mean.sd.en,file=file.name)
}