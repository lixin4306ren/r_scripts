block.find <- function(chip.bed.file, input.bed.file,genome.bed,index.name=NULL, bed.file.name=tmp.bed,norm.factor=1,sep=21000,cut.off.pro=0.3,cut.off.up=0.95, binsize=1000,verbose=TRUE){
		command=paste("cat ",chip.bed.file, " |genomic_scans counts -min 0 -i -v -w ", binsize, " -d ",binsize," -g ",genome.bed," |genomic_regions bed > ",index.name,".",binsize,".tmp1.density",sep="")
        try(system(command))
		command=paste("cat ",input.bed.file," |genomic_scans counts -min 0 -i -v -w ", binsize, " -d ",binsize," -g ",genome.bed," |genomic_regions bed > ",index.name,".",binsize,".tmp2.density",sep="")
        try(system(command))
		read.delim(paste(index.name,".",binsize,".tmp1.density",sep=""),header=FALSE, row.names=NULL,sep="")->chip
        read.delim(paste(index.name,".",binsize,".tmp2.density",sep=""),header=FALSE, row.names=NULL,sep="")->input
        
		chip$V4->chip.value
        input$V4->input.value
		threshold=get.cutoff2(chip.value, input.value,cut.off.pro,cut.off.up,norm.factor)
		
		ratio<-data.frame(chip[,1:3],log2(norm.factor*(chip$V4+1)/(input$V4+1)))
		names(ratio)<-c("chr","start","end","ratio")
		
		chip$V4+input$V4>10 -> index
	chr.result <- NULL
    allchr <- ratio$chr
	chrs <- sort(unique(as.vector(ratio$chr)))
    allstart <- ratio$start
	allend <- ratio$end
    for (j in 1:length(chrs)) {
        chr <- chrs[j]
        cat(chr, ". ")
        idx.chr <- which(allchr == chr)
        start.chr <- allstart[idx.chr]
		end.chr <- allend[idx.chr]
        M.chr <- ratio$ratio[idx.chr]
		
        chr.result.tmp <- NULL
		chr.result.tmp$start<-start.chr
		chr.result.tmp$end<-end.chr
		chr.result.tmp$ratio<-M.chr
        chr.result[[j]] <- chr.result.tmp
        names(chr.result)[j] <- chr
    }
	
	block<-NULL
	window.len=0
	for (i in 1:length(chr.result)) {
        if (verbose)cat(names(chr.result)[i], ".")

            win.start <- chr.result[[i]]$start
			win.end<-chr.result[[i]]$end
			y <- chr.result[[i]]$ratio
            n <- length(y)
            startidx <- NULL
            endidx <- NULL
            
			if (y[1] > threshold )startidx <- 0
            startidx <- c(startidx, which((y[1:n - 1] <= threshold ) & (y[2:n] > threshold )))
            startidx=startidx+1
			endidx <- which((y[1:n - 1] > threshold ) & (y[2:n] <=threshold ))
            if (y[length(y)] > threshold) endidx <- c(endidx, length(y))
			
			tmp.start<-win.start[startidx]
			tmp.end<-win.end[endidx]
			tmp.length<-tmp.end-tmp.start
            tmp <- cbind(tmp.start, tmp.end, tmp.length)
            block[[i]]<-tmp
			window.len=window.len+sum(tmp.length)
            names(block)[i]<-chrs[i]
    }
	
	block.merge<-NULL
	block.merge<-merge.block(block=block,sep=sep,bed.file.name=bed.file.name,chrs)
        cat ("threshold:",threshold,"\n")
        cat("befor merge LOCKs length:",window.len,"  ",window.len/(length(chip$V4[index])*binsize),"\n")
        cat("LOCKs length:",block.merge$block.length,"  ",block.merge$block.length/(length(chip$V4[index])*binsize),"\n")
        cat("LOCKs range:",range(block.merge$range))
}

get.blank.pos<-function(chr.result,threshold,chrs){
	block<-NULL
	for (i in 1:length(chr.result)) {
			names(chr.result[[i]])<-c("start","end","ratio")
            win.start <- chr.result[[i]]$start
			win.end<-chr.result[[i]]$end
			y <- chr.result[[i]]$ratio
            n <- length(y)
            startidx <- NULL
            endidx <- NULL
        	if (y[1] <= threshold ){startidx<-0}
            startidx <- c(startidx,which((y[1:n - 1] > threshold ) & (y[2:n] <= threshold )))
            startidx=startidx+1
			endidx <- which((y[1:n - 1] <= threshold ) & (y[2:n] > threshold ))
			if (y[length(y)] <= threshold) {endidx<-c(endidx,length(y))}
            if (y[1] <= threshold ){
				startidx<-startidx[-1]
				endidx<-endidx[-1]
			}
			if (y[length(y)] <= threshold) {
				startidx<-startidx[1:length(startidx)-1]
				endidx <- endidx[1:length(endidx)-1]
			}
			
            tmp <- cbind(win.start[startidx], win.end[endidx], (win.end[endidx] - win.start[startidx]))
            block[[i]]<-tmp
			#window.len=window.len+sum((win.end[endidx] - win.start[startidx]))
            names(block)[i]<-chrs[i]
    }
	return(block)
}

get.block.pos<-function(chr.result,threshold,chrs){
	block<-NULL
	for (i in 1:length(chr.result)) {
			names(chr.result[[i]])<-c("start","end","ratio")
            win.start <- chr.result[[i]]$start
			win.end<-chr.result[[i]]$end
			y <- chr.result[[i]]$ratio
            n <- length(y)
            startidx <- NULL
            endidx <- NULL
        	if (y[1] > threshold )startidx <- 0
            startidx <- c(startidx, which((y[1:n - 1] <= threshold ) & (y[2:n] > threshold )))
            startidx=startidx+1
			endidx <- which((y[1:n - 1] > threshold ) & (y[2:n] <=threshold ))
            if (y[length(y)] > threshold) endidx <- c(endidx, length(y))
			
            tmp <- cbind(win.start[startidx], win.end[endidx], (win.end[endidx] - win.start[startidx]))
            block[[i]]<-tmp
			window.len=window.len+sum((win.end[endidx] - win.start[startidx]))
            names(block)[i]<-chrs[i]
    }
	return(block)
}

merge.block <- function (block,sep = 21000,bed.file.name=tmp,chrs)
{
	if(file_test("-f",bed.file.name)){unlink(bed.file.name)}
    res<-NULL
	block.range<-NULL
	block.len=0
	for(i in 1:length(block)){
		if(length(block[[i]])==0){next}
		which(diff(block[[i]][,1])>sep)->idx.diff
		if(length(idx.diff)==0){
			start<-min(block[[i]][,1])
			end<-max(block[[i]][,2])
		}else{
			start<-c(min(block[[i]][,1]),block[[i]][idx.diff+1,1])
			end<-c(block[[i]][idx.diff,2],max(block[[i]][idx.diff,2]))
		}
		(end-start) >=20000 ->index
		res.tmp<-NULL
		res.tmp$start<-start[index]
		res.tmp$end<-end[index]
		block.len=block.len+sum(end[index]-start[index])
		res[[i]]<-res.tmp
		names(res)[i]<-chrs[i]
		#cat (i,"\n")
		if(length(res.tmp$start)>0){block.range<-rbind(block.range,range(end[index]-start[index]))}
		write.table(data.frame(rep(chrs[i],length(res[[i]]$start)),format(res[[i]]$start,scientific=FALSE,trim=TRUE),format(res[[i]]$end,scientific=FALSE,trim=TRUE)),file=bed.file.name,quote=F,col.names=F,row.names=F,append=T,sep="\t")
	}

	return(list(blocks=res,block.length=block.len,range=block.range))
}

get.cutoff2 <- function(chip.value, input.value,cut.off.down,cut.off.up,norm.factor){
    total.value <- chip.value+input.value
	ind<- total.value > 10 
	total.value<- total.value[ind]
	chip.value<-chip.value[ind]
	input.value<-input.value[ind]
	total.value<quantile(total.value,cut.off.down)->index
	chip.value[index]->chip.value2
	input.value[index]->input.value2
	(chip.value2 >0 & input.value2 > 0)->index2
	quantile(log2(norm.factor*(chip.value2[index2]+1)/(input.value2[index2]+1)),cut.off.up)->cutoff
	return(cutoff)
}
