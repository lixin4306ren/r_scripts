t.test.combine.bsseq<-function(i,meth,direction,correlated=TRUE){
require(RnBeads)
out<-NULL;
mapply(t.test.bsseq,tmp.list[[as.character(i)]]$subjectHits,MoreArgs=list(meth=meth,direction=direction))->t.test.pvalue
out<-combineTestPvalsMeth(t.test.pvalue, correlated = correlated)
print(i)
out
}

t.test.bsseq.no.direction<-function(i,meth,paired=TRUE){
out<-NULL;
out<-t.test(meth[i,group1],meth[i,group2],paired=paired)$p.value
out
}


t.test.bsseq<-function(i,meth,direction){
out<-NULL;
if(direction=="hypo"){
out<-t.test(meth[i,group1],meth[i,group2],paired=TRUE,alternative="less")$p.value
}
else{
out<-t.test(meth[i,group1],meth[i,group2],paired=TRUE,alternative="greater")$p.value
}
out
}
chm.combine.test.bsseq<-function(i){
library(lawstat)
out<-NULL
sapply(tmp.list[[i]]$subjectHits,chm.test.bsseq,un=a,m=b)->chi.stat
out<-combineTestPvalsMeth(pchisq(chi.stat,1,lower.tail=F), correlated = TRUE)
print(i)
out
}
