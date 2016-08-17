NCIS.internal <- function(chip.bed.file, input.bed.file,genome.bed,index.name=NULL, cut.off.pro=0.3,min.binsize=100, max.binsize=20000, min.stop.binsize=100){

    sizevec <- rep(c(1,2,5), times = 4)*10^rep(2:5, each =3 )
    sizevec<-sort(c(rep(c(1, 2, 5), times = 4)*10^rep(2:5, each =3 ),15000))
    sizevec <- sizevec[sizevec <= max.binsize]
    sizevec <- sizevec[sizevec >= min.binsize]
    norm.est <- rep(1000, length(sizevec))

	binsize.est <- -1
    pdf(paste(index.name,".pdf",sep=""))
    for(si in 1:length(sizevec)){
        binsize <- sizevec[si]
	print(paste("binsize=",binsize,sep=" "))
	#stopifnot(FALSE)
	command=paste("cat ",chip.bed.file, " |genomic_scans counts -min 0 -i -v -op c -w ", binsize, " -d ",binsize," -g ",genome.bed," |genomic_regions bed > ",index.name,".",binsize,".tmp1.density",sep="")
	try(system(command))
	command=paste("cat ",input.bed.file," |genomic_scans counts -min 0 -i -v -op c -w ", binsize, " -d ",binsize," -g ",genome.bed," |genomic_regions bed > ",index.name,".",binsize,".tmp2.density",sep="")
	try(system(command))
	print("read density file")
		read.delim(paste(index.name,".",binsize,".tmp1.density",sep=""),header=FALSE, row.names=NULL,sep="")->chip
		read.delim(paste(index.name,".",binsize,".tmp2.density",sep=""),header=FALSE, row.names=NULL,sep="")->input

		total<-list()
		chip$V4+input$V4-> total
		ind <- total>0
		chip$V4[ind]->chip
		input$V4[ind]->input
		total[ind]->total

		nchip <- sum(chip)
		ninput <- sum(input)
		r.seq.depth <- nchip/ninput
	print("estimate norm factor");
        norm.est[si] <- est.norm.med.search(chip, input,cut.off.pro)
	print (norm.est[si])
	print ("plot figure")
	NCIS.plot(chip,input,binsize,index.name,cut.off.pro)
	title(main=paste("binsize=",binsize,"bp",sep=""))
	mtext(paste("r=",norm.est[si],sep=""))
	abline(h=log10(norm.est[si]))
	unlink(paste(index.name,".",binsize,".tmp2.density",sep=""))
	unlink(paste(index.name,".",binsize,".tmp1.density",sep=""))
	#stopping criteria
        if(si>1 & binsize.est<0){
            if(norm.est[si]>=norm.est[si-1]){
                est <- norm.est[si-1]
                binsize.est <- sizevec[si-1]
            }
            if(si==length(sizevec) & binsize.est<0){ #the end of binsize, no converge yet
                est <- norm.est[si]
                binsize.est <- sizevec[si]
            }
        } #end if(si>1 & binsize.est<0)
        if(binsize.est>0 & binsize>=min.stop.binsize){
            break
        }
    } #end for(si in 1:length(sizevec))
	dev.off()
    return(list(est=est, binsize.est=binsize.est,r.seq.depth=r.seq.depth, pi0=est/r.seq.depth))
}	
	
est.norm.med.search <- function(chip, input,cut.off){
   
    total <- chip+input
    tbl <- table(total)
    total.count <- as.integer(names(tbl))
    cum.bc <- cumsum(tbl)
    cum.prop <- cum.bc/length(total)

    if(cum.prop[1] > cut.off){
        threshold <- 1
    }else{
        #largest total before median
        threshold <- max(total.count[cum.prop < cut.off])
    }
    ind <- total <= threshold
    chip.sum.low <- sum(chip[ind])
    input.sum.low <- sum(input[ind])
    bin.count.low <- cum.bc[which(total.count==threshold)]

    chip <- chip[!ind]
    input <- input[!ind]
    od <- order(chip+input)
    chip <- chip[od]
    input <- input[od]

    #start after threshold
    cum.bc.high <- cum.bc[total.count>threshold]-bin.count.low
    
    all.cum.chip <- cumsum(chip)
    all.cum.input <- cumsum(input)
    cum.chip <- c(0, all.cum.chip[cum.bc.high])
    cum.input <- c(0, all.cum.input[cum.bc.high])
    cum.ratio <- (chip.sum.low+cum.chip)/(input.sum.low+cum.input)

    if(sum(diff(cum.ratio) >= 0) ==0) stop("No increase in cum.ratio, binding signal missing!")
    index <- min((2:length(cum.ratio))[diff(cum.ratio) >= 0])
    return(cum.ratio[index])
}

NCIS.plot<-function (chip,input,plot.file.name,index.name,cut.off){
    total <- chip+input
    tbl <- table(total)
    total.count <- as.integer(names(tbl))
    cum.bc <- cumsum(tbl)
    cum.prop <- cum.bc/length(total)
    threshold2 <- max(total.count[cum.prop < cut.off])
    threshold <- 0
 
    ind <- total <= threshold
    chip.sum.low <- sum(chip[ind])
    input.sum.low <- sum(input[ind])
    bin.count.low <- cum.bc[which(total.count==threshold)]

    chip <- chip[!ind]
    input <- input[!ind]
    od <- order(chip+input)
    chip <- chip[od]
    input <- input[od]
	total <- total[od]

    #start after threshold
    cum.bc.high <- cum.bc[total.count>threshold]
    
    all.cum.chip <- cumsum(chip)
    all.cum.input <- cumsum(input)
	all.cum.total <- cumsum(total)
    cum.chip <- c(0, all.cum.chip[cum.bc.high])
    cum.input <- c(0, all.cum.input[cum.bc.high])
	plot(log10((cum.chip+cum.input)[-1]),log10(diff(cum.chip)/diff(cum.input)),xlab="log10(total reads)",ylab="log10(chip/input)")
	abline(v=log10(all.cum.total[cum.bc[as.character(threshold2)]]))
}
