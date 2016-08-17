extend.target.region<-function(region,probe.len=60,min.len=1){
library(GenomicRanges)
width(region) >= min.len ->ind
region[ind]->region
reduce(region)->region
width(region) < probe.len -> ind
region[ind,]->region.short
resize(region.short,probe.len,fix="center")->region.short
width(region)>= probe.len -> ind
region[ind,]->region.long
reduce(c(region.short,region.long))->new
return(new)
}

extract.cpg.cluster<-function(all.cpg,region,gap=100,min.cpg.num=1){
require(GenomicRanges)
require(bumphunter)

all.cpg %over% region ->ind
all.cpg[ind]->cpg.block.range
list(chr=as.vector(seqnames(cpg.block.range)),pos=start(cpg.block.range))->cpg.block
pns=clusterMaker(cpg.block$chr,cpg.block$pos,maxGap=gap)
y=tapply(cpg.block$pos,pns,function(x) diff(range(x))+1) # region length
z=tapply(cpg.block$pos,pns, length) #how many cpgs within one region

start=tapply(cpg.block$pos,pns,min)
end=tapply(cpg.block$pos,pns,max)
name=tapply(cpg.block$chr,pns,unique)
range=GRanges(seqnames=name,range=IRanges(start=start,end=end))
list(region.len=y,region.cpg.num=z,range=range)->region.cluster

region.cluster$region.cpg.num >= min.cpg.num ->ind
region.cluster$region.len[ind]->y
region.cluster$region.cpg.num[ind]->z
region.cluster$range[ind]->range

list(region.len=y,region.cpg.num=z,range=range)->region.cluster
return(region.cluster)
}
