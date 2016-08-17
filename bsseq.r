read.mydata2<-function(file,sort=TRUE){
    library("bsseq")

    dat <- read.table(file, , row.names = NULL,
                      col.names = c("chr", "pos", "strand", "M", "Cov"),
                      colClasses = c("character", "integer", "integer",
                      "integer", "integer"))

    dat1 <- dat[dat$strand == 1,]
    if(sort==TRUE){dat1<-dat1[order(dat1$chr,dat1$pos),]}
    dat2 <- dat[dat$strand == 4,]
    if(sort==TRUE){dat2<-dat2[order(dat2$chr,dat2$pos),]}

    if (!identical (as.integer(dat1$pos),as.integer(dat2$pos-1))){stop ("unequal length for + - strand\n");}

    tmp <- dat1
    BS.forward <- BSseq(pos = tmp$pos, chr = tmp$chr, M = as.matrix(tmp$M, ncol = 1),Cov = as.matrix(tmp$Cov, ncol = 1), sampleNames = "forward")
    tmp <- dat2
    BS.reverse <- BSseq(pos = tmp$pos - 1L, chr = tmp$chr, M = as.matrix(tmp$M, ncol = 1),Cov = as.matrix(tmp$Cov, ncol = 1), sampleNames = "reverse")
    BS <- combine(BS.forward, BS.reverse)
    BS <- collapseBSseq(BS, columns = c("a", "a"))
    BS
}

read.mydata.noncg<-function(file,sort=TRUE){
    library("bsseq")

    dat <- read.table(file, , row.names = NULL,
                      col.names = c("chr", "pos", "strand", "M", "Cov"),
                      colClasses = c("character", "integer", "integer",
                      "integer", "integer"))

    tmp <- dat
    BS <- BSseq(pos = tmp$pos, chr = tmp$chr, M = as.matrix(tmp$M, ncol = 1),Cov = as.matrix(tmp$Cov, ncol = 1), sampleNames = "forward")
    BS
}

read.mydata3<-function(dat,sort=TRUE){
    library("bsseq")

    dat1 <- dat[dat$strand == 1,]
    if(sort==TRUE){dat1<-dat1[order(dat1$chr,dat1$pos),]}
    dat2 <- dat[dat$strand == 4,]
    if(sort==TRUE){dat2<-dat2[order(dat2$chr,dat2$pos),]}

    if (!identical (as.integer(dat1$pos),as.integer(dat2$pos-1))){stop ("unequal length for + - strand\n");}

    tmp <- dat1
    BS.forward <- BSseq(pos = tmp$pos, chr = tmp$chr, M = as.matrix(tmp$M, ncol = 1),
                        Cov = as.matrix(tmp$Cov, ncol = 1), sampleNames = "forward")
    tmp <- dat2
    BS.reverse <- BSseq(pos = tmp$pos - 1L, chr = tmp$chr, M = as.matrix(tmp$M, ncol = 1),
                        Cov = as.matrix(tmp$Cov, ncol = 1), sampleNames = "reverse")
    BS <- combine(BS.forward, BS.reverse)
    BS <- collapseBSseq(BS, columns = c("a", "a"))
    BS
}

merge.cg.file<-function(sample,chr.list="chr.list"){
as.vector(read.table(chr.list)$V1)->chr
out<-data.frame()
for(i in 1:length(chr)){
	file.name=paste(sample,path,chr[i],".cg",sep="")
	print(file.name)
	read.table(file.name,colClasses = c("integer", "integer", "integer"))->tmp
	tmp.column<-rep(chr[i],nrow(tmp))
	cbind(tmp.column,tmp)->tmp
	rbind(out,tmp)->out
}
colnames(out)<-c("chr", "pos", "strand", "M", "Cov")
return(out)
}

convert2bsseq<-function(sample.list="sample.list",chr.list="chr.list"){
read.table(sample.list)->s
as.vector(read.table(chr.list)$V1)->chr
print(s[1,1])
merge.cg.file(s[1,1],chr.list)->tmp.all.cg
BS<-read.mydata3(tmp.all.cg)
sampleNames(BS)<-s[1,1]
for(i in 2:nrow(s)){
print (s[i,1])
merge.cg.file(s[i,1],chr.list)->tmp.all.cg
read.mydata3(tmp.all.cg)->tmp.bs
sampleNames(tmp.bs)<-s[i,1]
BS<-combine(BS,tmp.bs)
}
save(BS, file = "BS.all.rda")
}


plotManyRegions2<-function (BSseq, regions = NULL, extend = 0, main = "", addRegions = NULL,
    annoTrack = NULL, col = NULL, lty = NULL, lwd = NULL, BSseqTstat = NULL,
    mainWithWidth = TRUE, regionCol = alpha("red", 0.1), addTicks = TRUE,
    addPoints = FALSE, pointsMinCov = 5, highlightMain = FALSE,
    verbose = TRUE)
{
    cat("preprocessing ...")
    if (!is.null(regions)) {
        if (is(regions, "data.frame"))
            gr <- data.frame2GRanges(regions, keepColumns = FALSE)
        else gr <- regions
        if (!is(gr, "GRanges"))
            stop("'regions' needs to be either a 'data.frame' (with a single row) or a 'GRanges' (with a single element)")
    }
    else {
        gr <- granges(BSseq)
    }
    gr <- resize(gr, width = 2 * extend + width(gr), fix = "center")
    BSseq <- subsetByOverlaps(BSseq, gr)
    if (!is.null(BSseqTstat))
        BSseqTstat <- subsetByOverlaps(BSseqTstat, gr)
    if (length(start(BSseq)) == 0)
        stop("No overlap between BSseq data and regions")
    if (!is.null(main) && length(main) != length(gr))
        main <- rep(main, length = length(gr))
    cat("done\n")
    for (ii in seq(along = gr)) {
        if (verbose)
            cat(sprintf("plotting region %d (out of %d)\n", ii,
                nrow(regions)))
        plotRegion2(BSseq = BSseq, region = regions[ii, ], extend = extend,
            col = col, lty = lty, lwd = lwd, main = main[ii],
            BSseqTstat = BSseqTstat, addRegions = addRegions,
            regionCol = regionCol, mainWithWidth = mainWithWidth,
            annoTrack = annoTrack, addTicks = addTicks, addPoints = addPoints,
            pointsMinCov = pointsMinCov, highlightMain = highlightMain)
    }
}

plotRegion2<-function (BSseq, region = NULL, extend = 0, main = "", addRegions = NULL,
    annoTrack = NULL, col = NULL, lty = NULL, lwd = NULL, BSseqTstat = NULL,
    mainWithWidth = TRUE, regionCol = alpha("red", 0.1), addTicks = TRUE,
    addPoints = FALSE, pointsMinCov = 5, highlightMain = FALSE)
{
    plotRects <- function(ylim) {
        if (!is.null(addRegions))
            rect(xleft = addRegions$start, xright = addRegions$end,
                ybottom = ylim[1], ytop = ylim[2], col = regionCol,
                border = NA)
    }
    restrictRegions <- function(regions, plotRange, plotChr) {
        if (is.null(regions))
            return(NULL)
        regions <- regions[regions$chr == plotChr & ((regions$start >=
            plotRange[1] & regions$start <= plotRange[2]) | (regions$end >=
            plotRange[1] & regions$end <= plotRange[2])), , drop = FALSE]
        if (nrow(regions) == 0)
            regions <- NULL
        regions
    }
    plotLines <- function(x, y, lty, col, lwd, plotRange) {
        if (sum(!is.na(y)) <= 1)
            return(NULL)
        xx <- seq(from = plotRange[1], to = plotRange[2], length.out = 2000)
        yy <- approxfun(x, y)(xx)
        lines(xx, yy, col = col, lty = lty, lwd = lwd)
    }
    plotPoints <- function(x, y, z, col) {
        points(x[z > pointsMinCov], y[z > pointsMinCov], col = col,
            pch = 16, cex = 0.5)
    }
    if (!is.null(region)) {
        if (is(region, "data.frame"))
            gr <- data.frame2GRanges(region, keepColumns = FALSE)
        else gr <- region
        if (!is(gr, "GRanges") || length(gr) != 1)
            stop("'region' needs to be either a 'data.frame' (with a single row) or a 'GRanges' (with a single element)")
    }
    else {
        gr <- GRanges(seqnames = seqnames(BSseq)[1], ranges = IRanges(start = min(start(BSseq)),
            end = max(start(BSseq))))
    }
    origWidth <- width(gr)
    gr <- resize(gr, width = 2 * extend + width(gr), fix = "center")
    plotRange <- c(start(gr), end(gr))
    plotChr <- as.character(seqnames(gr))[1]
    BSseq <- subsetByOverlaps(BSseq, gr)
    if (!is.null(BSseqTstat))
        BSseqTstat <- subsetByOverlaps(BSseqTstat, gr)
    positions <- start(BSseq)
    if (length(positions) == 0)
        stop("No overlap between BSseq data and region")
    regionCoord <- sprintf("%s: %s - %s", plotChr, format(plotRange[1],
        big.mark = ",", scientific = FALSE), format(plotRange[2],
        big.mark = ",", scientific = FALSE))
    regionWidth <- sprintf("width = %s, extended = %s", format(origWidth,
        big.mark = ",", scientific = FALSE), format(extend, big.mark = ",",
        scientific = FALSE))
    if (mainWithWidth)
        regionCoord <- sprintf("%s (%s)", regionCoord, regionWidth)
    if (is.null(main)) {
        main <- ""
    }
    else {
        if (main != "") {
            main <- sprintf("%s\n%s", main, regionCoord)
        }
        else {
            main <- regionCoord
        }
    }
    opar <- par(mar = c(0, 4.1, 0, 12), oma = c(5, 0, 4, 2), mfrow = c(1,
        1),xpd=TRUE)
    on.exit(par(opar))
    if (is.null(BSseqTstat))
        layout(matrix(1:2, ncol = 1), heights = c(2, 1))
    else layout(matrix(1:3, ncol = 1), heights = c(2, 2, 1))
    sampleNames <- sampleNames(BSseq)
    names(sampleNames) <- sampleNames
    plot(positions[1], 0.5, type = "n", xaxt = "n", yaxt = "n",
        ylim = c(0, 1), xlim = plotRange, xlab = "", ylab = "Methylation")
    axis(side = 2, at = c(0.2, 0.5, 0.8))
    if (addTicks)
        rug(positions)
    addRegions <- restrictRegions(addRegions, plotRange = plotRange,
        plotChr = plotChr)
    if (highlightMain)
        addRegions <- rbind(region[, c("chr", "start", "end")],
            addRegions[, c("chr", "start", "end")])
    if (!is.null(addRegions))
        plotRects(c(0, 1))
    smoothPs <- getMeth(BSseq, type = "smooth")
    rawPs <- getMeth(BSseq, type = "raw")
    coverage <- getCoverage(BSseq)
    if (is.null(col) & "col" %in% names(pData(BSseq)))
        col <- pData(BSseq)[["col"]]
    else col <- rep("black", ncol(BSseq))
    if (is.null(names(col)))
        names(col) <- sampleNames(BSseq)
    if (is.null(lty) & "lty" %in% names(pData(BSseq)))
        lty <- pData(BSseq)[["lty"]]
    else lty <- rep(1, ncol(BSseq))
    if (is.null(names(lty)))
        names(lty) <- sampleNames(BSseq)
    if (is.null(lwd) & "lwd" %in% names(pData(BSseq)))
        lwd <- pData(BSseq)[["lwd"]]
    else lwd <- rep(1, ncol(BSseq))
    if (is.null(names(lwd)))
        names(lwd) <- sampleNames(BSseq)
    if (addPoints) {
        sapply(sampleNames(BSseq), function(samp) {
            abline(v = positions[rawPs[, samp] > 0.1], col = "grey80",
                lty = 1)
        })
    }
    sapply(sampleNames(BSseq), function(samp) {
        plotLines(positions, smoothPs[, samp], col = col[samp],
            lty = lty[samp], lwd = lwd[samp], plotRange = plotRange)
    })
    if (addPoints) {
        sapply(sampleNames(BSseq), function(samp) {
            plotPoints(positions, rawPs[, samp], coverage[, samp],
                col = col[samp])
        })
    }
    if (!is.null(BSseqTstat)) {
        plot(positions[1], 0.5, type = "n", xaxt = "n", yaxt = "n",
            ylim = c(-8, 8), xlim = plotRange, xlab = "", ylab = "t-stat")
        axis(side = 2, at = c(-5, 0, 5))
        abline(h = 0, col = "grey60")
        plotLines(start(BSseqTstat), BSseqTstat@stats[, "tstat"],
            lty = 1, plotRange = plotRange, col = "red", lwd = 1)
        plotLines(start(BSseqTstat), BSseqTstat@stats[, "tstat.corrected"],
            lty = 2, plotRange = plotRange, col = "red", lwd = 1)
        plotLines(start(BSseqTstat), 100 * BSseqTstat@stats[,
            "tstat.sd"], lty = 2, plotRange = plotRange, col = "blue",
            lwd = 1)
    }
    if (!is.null(annoTrack))
        bsseq:::plotAnnoTrack(gr, annoTrack)
    mtext(side = 3, text = main, outer = TRUE, cex = 1)
    legend("topright",legend=sampleNames(BSseq),text.col=col,inset=c(-0.3,0))
    return(invisible(NULL))
}


meth_cluster_plot<-function(sample.list="sample.list"){
read.table(sample.list,colClasses="character",fill=T)->s
out<-NULL
name<-NULL
for(i in 1:nrow(s)){
unlist(strsplit(s[i,1],"/"))->infor
print(infor[length(infor)])
name<-append(name,infor[length(infor)])
if(file.exists(paste(infor[length(infor)],"all.bedgraph",sep="/"))){
read.table(paste(infor[length(infor)],"all.bedgraph",sep="/"))->tmp
}
else{
read.table(paste(infor[length(infor)],"all/all.bedgraph",sep="/"))->tmp
}
out<-cbind(out,tmp[,4])
}
colnames(out)<-name
na.omit(out)->out
hclust(as.dist(1-cor(out)),method="ward")->h
paste(s[,2],s[,3],sep=".")->name

source("~/r_scripts/functions.R")
pdf(file="whole_meth_cluster.pdf")
myplclust(h,lab=name,lab.col=as.integer(as.factor(name)),main="",xlab="")
dev.off()
}


dmr_pipeline<-function(sample.list="sample.list",cell.col=2,tissue.col=NULL,min.sample.num=2,mean.diff.large=0.1,mean.diff.small=0.2,m.core=6,path="/all/",cg_file_name="all.cg",name="all"){
source("~/r_scripts/bsseq.r")
library(bsseq)

############### generate BS.all.rda
if(!file.exists(paste("BS",name,".rda",sep=""))){
read.table(sample.list,colClasses="character",fill=T)->s
unlist(strsplit(s[1,1],"/"))->infor
print(infor[length(infor)])
BS<-read.mydata2(paste(infor[length(infor)],path,cg_file_name,sep=""))
sampleNames(BS)<-infor[length(infor)]

if(nrow(s)==1){
print ("saving BS.all.rda ...")
save(BS, file = paste("BS.",name,".rda",sep=""))
}
else{
for(i in 2:nrow(s)){
unlist(strsplit(s[i,1],"/"))->infor
print(infor[length(infor)])
read.mydata2(paste(infor[length(infor)],path,cg_file_name,sep=""))->tmp.bs
sampleNames(tmp.bs)<-infor[length(infor)]
BS<-combine(BS,tmp.bs)
}
print ("saving BS.all.rda ...")
save(BS, file = paste("BS.",name,".rda",sep=""))
}
}

################ generate bs.fit object
if(!file.exists("BS.fit.rda")){
load("BS.all.rda")
library(bsseq)
BS.fit.small<-BSmooth(BS,mc.cores=m.core,verbose=TRUE,n=20,h=1000,parallelBy = "chromosome")
BS.fit.large<-BSmooth(BS,mc.cores=m.core,verbose=TRUE,n=200,h=10000)
print ("saving BS.fit object ...")

save(list=c("BS.fit.small","BS.fit.large"),file="BS.fit.rda")
}

################ generate bs.stat object
if(!file.exists("BS.tstat.rda")){
load("BS.fit.rda")
read.table(sample.list)->sample
unique(sample[,tissue.col])->tissues
unique(sample[,cell.col])->cells

BS.tstat<-NULL
index<-1

if(!is.null(tissue.col)){
for(i in 1:length(tissues)){
which(sample[,2]==cells[1] & sample$V3==tissues[i])->ind.h
which(sample[,2]==cells[2] & sample$V3==tissues[i])->ind.l
as.character(sampleNames(BS.fit.small)[ind.h])->s.h
as.character(sampleNames(BS.fit.small)[ind.l])->s.l

BS.cov <- getCoverage(BS.fit.large)
keepLoci <- which(rowSums(BS.cov[, s.h] >= 3) >= min.sample.num & rowSums(BS.cov[, s.l] >= 3) >= min.sample.num & (seqnames(BS.fit.large@gr) !="chr32" & seqnames(BS.fit.large@gr) !="chrW"))
BS.fit.large.tmp <- BS.fit.large[keepLoci,]
BS.large.tstat.tmp <- BSmooth.tstat(BS.fit.large.tmp, group1 = s.h,group2 = s.l,estimate.var = "same",local.correct = TRUE,verbose =TRUE)
paste(tissues[i],"large",sep=".")->name
list(BS.large.tstat.tmp)->BS.tstat[index]
names(BS.tstat)[index]<-name
index<-index+1

BS.cov <- getCoverage(BS.fit.small)
keepLoci <- which(rowSums(BS.cov[, s.h] >= 3) >= min.sample.num & rowSums(BS.cov[, s.l] >= 3) >= min.sample.num & (seqnames(BS.fit.small@gr) !="chr32" & seqnames(BS.fit.small@gr) !="chrW"))
BS.fit.small.tmp <- BS.fit.small[keepLoci,]
BS.small.tstat.tmp <- BSmooth.tstat(BS.fit.small.tmp,group1 = s.h,group2 = s.l,estimate.var = "same",local.correct = TRUE,verbose = TRUE)
paste(tissues[i],"small",sep=".")->name
list(BS.small.tstat.tmp)->BS.tstat[index]
names(BS.tstat)[index]<-name
index<-index+1
}
}else{

for(i in 1:(length(cells)-1)){
for(j in (i+1):length(cells)){
which(sample[,2]==cells[i] )->ind.h
which(sample[,2]==cells[j] )->ind.l

as.character(sampleNames(BS.fit.small)[ind.h])->s.h
as.character(sampleNames(BS.fit.small)[ind.l])->s.l

BS.cov <- getCoverage(BS.fit.large)
keepLoci <- which(rowSums(BS.cov[, s.h] >= 3) >= min.sample.num & rowSums(BS.cov[, s.l] >= 3) >= min.sample.num & (seqnames(BS.fit.large@gr) !="chr32" & seqnames(BS.fit.large@gr) !="chrW"))
BS.fit.large.tmp <- BS.fit.large[keepLoci,]
BS.large.tstat.tmp <- BSmooth.tstat(BS.fit.large.tmp, group1 = s.h,group2 = s.l,estimate.var = "same",local.correct = TRUE,verbose =TRUE)
paste(cells[i],cells[j],"large",sep=".")->name
print(name)
list(BS.large.tstat.tmp)->BS.tstat[index]
names(BS.tstat)[index]<-name
index<-index+1

BS.cov <- getCoverage(BS.fit.small)
keepLoci <- which(rowSums(BS.cov[, s.h] >= 3) >= min.sample.num & rowSums(BS.cov[, s.l] >= 3) >= min.sample.num & (seqnames(BS.fit.small@gr) !="chr32" & seqnames(BS.fit.small@gr) !="chrW"))
BS.fit.small.tmp <- BS.fit.small[keepLoci,]
BS.small.tstat.tmp <- BSmooth.tstat(BS.fit.small.tmp,group1 = s.h,group2 = s.l,estimate.var = "same",local.correct = TRUE,verbose = TRUE)
paste(cells[i],cells[j],"small",sep=".")->name
print(name)
list(BS.small.tstat.tmp)->BS.tstat[index]
names(BS.tstat)[index]<-name
index<-index+1
}
}
}

print ("saving BS.tstat object ...")

save(BS.tstat,file="BS.tstat.rda")
}

################ generate DMRs object
if(!file.exists("BS.dmr.rda")){
load("BS.tstat.rda")
BS.dmr0<-NULL
BS.dmr<-NULL
name<-names(BS.tstat)
for(i in 1:length(name)) {
if(grepl("large",name[i])){
dmrs0.tmp <- dmrFinder(BS.tstat[[i]], cutoff = c(-2, 2),column="tstat",maxGap=10000)
dmrs.tmp <- subset(dmrs0.tmp, n >= 3 & abs(meanDiff) >= mean.diff.large & width>5000)
}
else{
dmrs0.tmp <- dmrFinder(BS.tstat[[i]], cutoff = c(-4.6, 4.6))
dmrs.tmp <- subset(dmrs0.tmp, n >= 3 & abs(meanDiff) >= mean.diff.small)
}
list(dmrs0.tmp)->BS.dmr0[i]
list(dmrs.tmp)->BS.dmr[i]
names(BS.dmr0)[i]<-name[i]
names(BS.dmr)[i]<-name[i]
}
print ("saving BS.dmr object ...")

save(BS.dmr0,file="BS.dmr0.rda")
save(BS.dmr,file="BS.dmr.rda")
}
}

sliding.count.region<-function(dat,chr=chr,start.pos=NULL,end.pos=NULL,colname=NULL,window.size,step.size,min.num=20,min.mc.depth=0){
library(data.table)
library(GenomicRanges)
library(GenomicFeatures)
#38      1       0       0
#39      4       0       0
#84      1       0       0
dat->x
x[x[,4] > 0,]->x
x[x[,3] > min.mc.depth,]->mcg

seq(start.pos,end.pos-window.size+1,step.size)->start
seq(start.pos+window.size-1,end.pos,step.size)->end
win<-GRanges(seqnames=chr,IRanges(start=start,end=end))
#transfer x, y to GRange objects
cg1<-GRanges(seqnames=chr,IRanges(start=x[,1],end=x[,1]))
mcg<-GRanges(seqnames=chr,IRanges(start=mcg[,1],end=mcg[,1]))

as.data.frame(findOverlaps(win,cg1))->tmp1
countOverlaps(win,mcg)->mcg
countOverlaps(win,cg1)->c

which(countOverlaps(win,cg1)==0)->tmp2 #regions without Cs

data.table(id=as.character(tmp1$queryHits),cov=x[,4][tmp1$subjectHits],m=x[,3][tmp1$subjectHits])->new1
setkey(new1,id)
as.integer(new1[,sum(cov),by=id]$id)->id1
data.table(id=new1[,sum(cov),by=id]$id,start=start(win)[id1],end=end(win)[id1],m=new1[,sum(m),by=id]$V1,cov=new1[,sum(cov),by=id]$V1)->tmp1
tmp1<-rbind(data.frame(tmp1),data.frame(id=as.character(tmp2),start=start(win)[tmp2],end=end(win)[tmp2],m=rep(0,length(tmp2)),cov=rep(0,length(tmp2))))
sort(as.integer(tmp1$id),index.return = T)->ind1
tmp1[ind1$ix,]->tmp1
tmp1$cov<min.num->ind
tmp1[ind,"cov"]<-0
tmp1[ind,"m"]<-0
data.frame(meth=tmp1$m/tmp1$cov)->out
colnames(out)<-colname
return(out)
}
