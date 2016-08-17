as.fumeric<- function(x,...) as.numeric(as.factor(x,...))
first<- function(x,k=1) x[k]
###CHECK THIS ~pmurakam/feinberg/CharmFiles/functions.R
boyorgirl <- function(A,xIndex,yIndex,plot=FALSE,id=1:ncol(A)){
  ##A is average of log int
  x=Biobase::rowMedians(t(A[xIndex,]),na.rm=TRUE)
  y=Biobase::rowMedians(t(A[yIndex,]),na.rm=TRUE)
  k=kmeans(y-x,c(min(y-x),max(y-x)))
  sex=factor(ifelse(k$cluster==which.min(k$centers),"F","M"),levels=c("M","F"))
  if(plot){
    plot(x,y,type="n")
    text(x,y,id,col=as.numeric(sex))
    legend("bottomleft",c("M","F"),col=c(1,2),pch=1)
  }
  return(sex)
}

splitit <- function(x) split(seq(along=x),x)

myplclust <- function( hclust, lab=hclust$labels, lab.col=rep(1,length(hclust$labels)), hang=0.1, ... )
{
  ## modifiction of plclust for plotting hclust objects *in colour*!
  ## Copyright Eva KF Chan 2009
  ## Arguments:
  ##    hclust:    hclust object
  ##    lab:        a character vector of labels of the leaves of the tree
  ##    lab.col:    colour for the labels; NA=default device foreground colour
  ##    hang:     as in hclust & plclust
  ## Side effect:
  ##    A display of hierarchical cluster with coloured leaf labels.
  y <- rep(hclust$height,2)
  x <- as.numeric(hclust$merge)
  y <- y[which(x<0)]
  x <- x[which(x<0)]
  x <- abs(x)
  y <- y[order(x)]
  x <- x[order(x)]
  plot( hclust, labels=FALSE,hang=hang,... )
    text( x=x, y=y[hclust$order]-(max(hclust$height)*hang), labels=lab[hclust$order],col=lab.col[hclust$order], srt=90, adj=c(1,0.5), xpd=NA, ... )
  
}

dogo <- function(names,universe,species="human"){
  if(species=="human"){
    golib="org.Hs.eg.db"
    library(golib,character.only=TRUE)
    gomap= org.Hs.egREFSEQ2EG
  } else  {
    golib="org.Mm.eg.db"
    library(golib,character.only=TRUE)
    gomap= org.Mm.egREFSEQ2EG
  }
  require(GOstats)
  x=unlist(mget(as.character(names), gomap,ifnotfound = NA))
  x=x[!is.na(x)]
  Universe=unlist(mget(as.character(universe),gomap,ifnotfound = NA))
  
  Universe=unique(c(Universe[!is.na(Universe)],unique(x)))
  
  params <- new("GOHyperGParams", geneIds = unique(x),
                universeGeneIds = Universe,
                annotation = golib,
                ontology = "BP", pvalueCutoff = 0.01, conditional = TRUE,
                testDirection="over")
  ht=hyperGTest(params)
  tab=summary(ht)
  tmp1=geneIdsByCategory(ht)
  tmp1=tmp1[tab[,1]]
  tab$IDs=sapply(tmp1,function(y) paste(names(x)[x%in%y],collapse=";"))
  return(tab)
}




affyGO <- function(names,universe=NULL,platform="hgu133plus2",species="human",pvalueCutoff=1){
  library(paste(platform,"db",sep="."),character.only=TRUE)
  gomap=get(paste(platform,"ENTREZID",sep=""))
  symmap=get(paste(platform,"SYMBOL",sep=""))
  if(species=="human"){
    golib="org.Hs.eg.db"
    library(golib,character.only=TRUE)
  } else  {
    golib="org.Mm.eg.db"
    library(golib,character.only=TRUE)
    gomap= org.Mm.egREFSEQ2EG
  }
  require(GOstats)
  x=unlist(mget(as.character(names), gomap,ifnotfound = NA))
  x=x[!is.na(x)]

  if(is.null(universe)) universe=ls(gomap)
  Universe=unlist(mget(universe,gomap,ifnotfound = NA))
  Universe=unique(c(Universe[!is.na(Universe)],unique(x)))
  
  params <- new("GOHyperGParams", geneIds = unique(x),
                universeGeneIds = Universe,
                annotation = golib,
                ontology = "BP", pvalueCutoff = pvalueCutoff,
                conditional = TRUE,
                testDirection="over")
  ht=hyperGTest(params)
  tab=summary(ht)
  tmp1=geneIdsByCategory(ht)
  tmp1=tmp1[tab[,1]]
  tab$IDs= sapply(tmp1,function(y) paste(unlist(mget(names(x)[x%in%y],symmap,ifnotfound=NA)),collapse=";"))
  return(tab)
}

myfilter2 <- function(x,filter,...){
###unlike myfilter, myfilter2 returns NAs in the edges.
  L=dim(x)[1]
  if(L>length(filter)) res=filter(x,filter,...) else{res = t(colMeans(x) %*% t(rep(1,dim(x)[1])))}

  return(res)

}


getDesc<-function(x){
 ##get gene description
  require(org.Hs.eg.db)
  genenames=sapply(mget(as.character(x),org.Hs.egREFSEQ2EG,ifnotfound = NA),function(x) x[1])
  genenames[!is.na(genenames)]= sapply(mget(genenames[!is.na(genenames)],org.Hs.egGENENAME,ifnotfound=NA),function(x) x[1])
  genenames
}


midpoints <- function(x) (x[1:(length(x)-1)]+x[2:length(x)])/2

pointsplot <- function(l,jitter=TRUE,factor=1,col=rep(1,length(l)),pch=21,...){
  if(!is.list(l)) stop("l must be a list\n")
  y=unlist(l)
  x=rep(seq(along=l),sapply(l,length))
  col=rep(col,sapply(l,length))
  if(jitter) x=jitter(x,factor)
  if(pch!=21) plot(x,y,col=col,pch=pch,...) else plot(x,y,pch=pch,bg=col,...)
}

tplot <- function(x,xlab="",ylab="",...){
  plot(as.numeric(names(x)),x,xlab=xlab,ylab=ylab,...)
}

scatterBox <- function(x,y,cuts=100,...) boxplot(split(y,cut(x,unique(quantile(x,seq(0,1,length=cuts+1))),include.lowest=TRUE)),...)


maplot <- function(x,y,...){ A=(x+y)/2;M=y-x;plot(A,M,...)}
                             
                             
mean.sliding<-function(x,chr=NULL,pos=2,col=3,win.size=10000){
start=1
end=max(x[,pos])
win.start=seq(1,end-win.size,by=win.size)
win.end=win.start+win.size-1
library(GenomicRanges)
library(data.table)

GRanges(seqnames=x[,1],ranges=IRanges(start=x[,pos],end=x[,pos]))->dat
win<-GRanges(seqnames=x[1,1],ranges=IRanges(start=win.start,end=win.end))
as.data.frame(findOverlaps(win,dat))->tmp1
which(countOverlaps(win,dat)==0)->tmp2
data.table(id=tmp1$queryHits,fst=x[tmp1$subjectHits,col])->new1
setkey(new1,id)
as.integer(new1[,mean(fst,na.rm=T),by=id]$id)->id1
data.table(id=id1,fst=new1[,mean(fst,na.rm=T),by=id]$V1)->tmp1
#tmp1<-rbind(data.frame(tmp1),data.frame(id=tmp2,fst=rep("NA",length(tmp2))))
#sort(tmp1$id,index.return = T)->ind1
#tmp1[ind1$ix,]->tmp1
return(tmp1)
}

###################
generate.sliding.window.center<-function(region,step.base=10,flanking=1000,from.start=TRUE){
library(GenomicRanges)

if(from.start==TRUE){
	flank(region,flanking,both=T)->region
}
else{
	flank(region,flanking,both=T,start=F)->region
}


out.start<-NULL
for(i in 1:((flanking*2)/step.base)){
round(start(region)+((end(region)-start(region)+1)*(step.base/(2*flanking))*(i-1)))->tmp
c(out.start,tmp)->out.start
}

matrix(out.start,nrow=length(region))->tbls.start

out.end<-NULL
for(i in 1:((flanking*2)/step.base)){
round(start(region)-1+(width(region)*(step.base/(2*flanking)*i)))->tmp
c(out.end,tmp)->out.end
}

matrix(out.end,nrow=length(region))->tbls.end

chr<-rep(seqnames(region),2*(flanking/step.base))
GRanges(seqnames=chr,range=IRanges(start=c(out.start),end=c(out.end)))->win
win
}

generate.sliding.window<-function(region,step.percent=0.01,step.base=10,up=1000,down=1000){
library(GenomicRanges)
chr<-rep(seqnames(region),1/step.percent+up/step.base+down/step.base)
start<-start(region)
end<-end(region)
generate.sliding.window.start(start,end,step.percent,step.base,up,down)->win.start
generate.sliding.window.end(start,end,step.percent,step.base,up,down)->win.end
GRanges(seqnames=chr,range=IRanges(start=c(win.start),end=c(win.end)))->win
win
}

###############generate start and end postion of region (percentage) and its flanking (base)
generate.sliding.window.start<-function(start,end,step.percent=0.01,step.base=10,up=1000,down=1000){
library(GenomicRanges)
if(up==0){
tbls.left<-NULL
}
else{
out.left<-NULL
for(i in 1:(up/step.base)){
round((start-up)+(step.base*(i-1)))->tmp
c(out.left,tmp)->out.left
}
matrix(out.left,nrow=length(start))->tbls.left
}

out.middle<-NULL

for(i in 1:(1/step.percent)){
round(start+((end-start+1)*step.percent*(i-1)))->tmp
c(out.middle,tmp)->out.middle
}
matrix(out.middle,nrow=length(start))->tbls.middle

if(down==0){
tbls.right<-NULL
}
else{
out.right<-NULL
for(i in 1:(down/step.base)){
round(end+(step.base*(i-1)))->tmp
c(out.right,tmp)->out.right
}
matrix(out.right,nrow=length(start))->tbls.right
}

cbind(tbls.left,tbls.middle,tbls.right)->tbls
tbls
}


generate.sliding.window.end<-function(start,end,step.percent=0.01,step.base=10,up=1000,down=1000){
library(GenomicRanges)

if(up==0){
tbls.left<-NULL
}
else{
out.left<-NULL
for(i in 1:(up/step.base)){
round((start-up-1)+(step.base*(i)))->tmp
c(out.left,tmp)->out.left
}
matrix(out.left,nrow=length(start))->tbls.left
}

out.middle<-NULL

for(i in 1:(1/step.percent)){
round(start-1+((end-start+1)*step.percent*(i)))->tmp
c(out.middle,tmp)->out.middle
}
matrix(out.middle,nrow=length(start))->tbls.middle

if(down==0){
tbls.right<-NULL
}
else{
out.right<-NULL
for(i in 1:(down/step.base)){
round(end-1+(step.base*(i)))->tmp
c(out.right,tmp)->out.right
}
matrix(out.right,nrow=length(start))->tbls.right
}

cbind(tbls.left,tbls.middle,tbls.right)->tbls
tbls
}


################## geneate sliding window positions for chrs
generate.sliding.window.chr<-function(region,win=1000,step=1000){
out<-NULL
	f<-function(x,a,b){
		region.start<-start(x)
		region.end<-end(x)
		seq(region.start,region.end-a+1,b)->start
		seq(region.start+a-1,region.end,b)->end
		GRanges(seqnames=seqnames(x),range=IRanges(start=start,end=end))
	}

out<-f(region[1],win,step)
for(i in 2:length(region)){
tmp.win<-f(region[i],win,step)
c(out,tmp.win)->out
}
return(out)
}


scatter.hist2<-function (x, y = NULL, smooth = TRUE, ab = FALSE, correl = TRUE,
    density = TRUE, ellipse = TRUE, digits = 2, method, cex.cor = 1,
    title = "Scatter plot + histograms", xlab = NULL, ylab = NULL,
    ...)
{
    old.par <- par(no.readonly = TRUE)
    if (missing(xlab)) {
        if (!is.null(colnames(x))) {
            xlab = colnames(x)[1]
            ylab = colnames(x)[2]
        }
        else {
            xlab = "V1"
            ylab = "V2"
        }
    }
    if (is.null(y)) {
        y <- x[, 2]
        x <- x[, 1]
    }
    else {
        if (!is.null(dim(x))) {
            x <- x[, 1, drop = TRUE]
            if (!is.null(colnames(y)))
                ylab <- colnames(y)
            if (!is.null(dim(y))) {
                y <- y[, 1, drop = TRUE]
            }
        }
    }
    xhist <- hist(x, breaks = 11, plot = FALSE)
    yhist <- hist(y, breaks = 11, plot = FALSE)
    xrange <- range(x, na.rm = TRUE)
    yrange <- range(y, na.rm = TRUE)
    nf <- layout(matrix(c(2, 4, 1, 3), 2, 2, byrow = TRUE), c(3,
        1), c(1, 3), TRUE)
    par(mar = c(5, 4, 1, 1))
    smoothScatter(x, y, xlim = xrange, ylim = yrange, xlab = xlab, ylab = ylab,
        ...)
    if (ab)
        abline(lm(y ~ x))
    if (smooth) {
        ok <- is.finite(x) & is.finite(y)
        if (any(ok))
            lines(stats::lowess(x[ok], y[ok]), col = "red")
    }
    if (ellipse) {
        ellipses(x, y, add = TRUE)
    }
    par(mar = c(0, 4, 2, 0))
    mp <- barplot(xhist$density, axes = FALSE, space = 0)
    tryd <- try(d <- density(x, na.rm = TRUE, bw = "nrd", adjust = 1.2),
        silent = TRUE)
    if (class(tryd) != "try-error") {
        d$x <- (mp[length(mp)] - mp[1] + 1) * (d$x - min(xhist$breaks))/(max(xhist$breaks) -
            min(xhist$breaks))
        if (density)
            lines(d)
    }
    title(title)
    par(mar = c(5, 0, 0, 2))
    mp <- barplot(yhist$density, axes = FALSE, space = 0, horiz = TRUE)
    tryd <- try(d <- density(y, na.rm = TRUE, bw = "nrd", adjust = 1.2),
        silent = TRUE)
    if (class(tryd) != "try-error") {
        temp <- d$y
        d$y <- (mp[length(mp)] - mp[1] + 1) * (d$x - min(yhist$breaks))/(max(yhist$breaks) -
            min(yhist$breaks))
        d$x <- temp
        if (density)
            lines(d)
    }
    par(mar = c(3, 1, 1, 1))
    if (correl) {
        plot(1, 1, type = "n", axes = FALSE)
        med.x <- median(x, na.rm = TRUE)
        med.y <- median(y, na.rm = TRUE)
        if (missing(method))
            method <- "pearson"
        r = (cor(x, y, use = "pairwise", method = method))
        txt <- format(c(r, 0.123456789), digits = digits)[1]
        if (missing(cex.cor)) {
            cex <- 0.75/strwidth(txt)
        }
        else {
            cex <- cex.cor
        }
        text(1, 1, txt, cex = cex)
    }
    par(old.par)
}


####################### countoverlap along sliding windows according to region range
calculate.number.sliding<-function(x,region,dist="base",win.base=5000,step.base=5000,up.base=200000,down.base=200000,win.percent=0.01,step.percent=0.01,up.percent=1,down.percent=1,min.num=20,min.mc.depth=0,normalize="Length",cpg.all.file="/amber3/feinbergLab/personal/xinli/no_back_up/data/R_data/hg19.cpg.all.rda",from.start=TRUE){
library(data.table)
library(GenomicRanges)

if(dist=="base"){
win<-generate.sliding.window(region,win.percent,step.base,up.base,down.base)
}

if(win.percent==1){
win<-generate.sliding.window.center(region,step.base,up.base,from.start=from.start)
}
print("starting")
countOverlaps(win,x)->count
if(normalize=="Length"){
matrix(count/width(win),nrow=length(region))->tbls
}
else if(normalize=="CpG"){
load(cpg.all.file)
countOverlaps(win,cpg.all)->count.all.cpg
count.all.cpg[count.all.cpg==0] <- NA
matrix(count/count.all.cpg,nrow=length(region))->tbls
}
else if (normalize=="None"){
matrix(count,nrow=length(region))->tbls
}

return(tbls)
}


####################### countoverlap according to region range
calculate.number.region<-function(x,region,normalize="Length",cpg.all.file="/amber3/feinbergLab/personal/xinli/no_back_up/data/R_data/hg19.cpg.all.rda"){
library(GenomicRanges)
region->win
countOverlaps(win,x)->count
if(normalize=="Length"){
matrix(count/width(win),nrow=length(region))->tbls
}
else if(normalize=="CpG"){
load("/amber3/feinbergLab/personal/xinli/no_back_up/data/R_data/hg19.cpg.all.rda")
countOverlaps(win,cpg.all)->count.all.cpg
count.all.cpg[count.all.cpg==0] <- NA
matrix(count/count.all.cpg,nrow=length(region))->tbls
}
else if (normalize=="None"){
matrix(count,nrow=length(region))->tbls
}

return(tbls)
}

####################### countoverlap overlapping length along sliding windows according to region range
calculate.length.sliding<-function(x,region,dist="base",win.base=5000,step.base=5000,up.base=200000,down.base=200000,win.percent=0.01,step.percent=0.01,up.percent=1,down.percent=1){
library(data.table)
library(GenomicRanges)

if(dist=="base"){
win<-generate.sliding.window(region,win.percent,step.base,up.base,down.base)
}

if(win.percent==1){
win<-generate.sliding.window.center(region,step.base,up.base,from.start=TRUE)
}
print("starting")
source("~/r_scripts/chip.r")
end(win)<=0 ->ind
end(win[ind])<- 1
start(win)<=0 ->ind
start(win[ind])<- 1

bedtools.coveragebed.range(bed1=x,bed2=win,mcol=FALSE)->cov

cov[,7]->count
matrix(count,nrow=length(region))->tbls

return(tbls)
}
