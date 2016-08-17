regionFinder3 <- function(x, chr, positions, keep, maxGap = 300, verbose = TRUE) {
    getSegments3 <- function(x, cid, verbose = TRUE){
        if(verbose) cat("[regionFinder3] segmenting\n")
        segments <- cumsum(c(1, diff(x) != 0))+cumsum(c(1, diff(cid) != 0))
        names(segments) <- NULL
        if(verbose) cat("[regionFinder3] splitting\n")
        out <- list(up = split(which(x > 0), segments[x > 0]),
                    down = split(which(x < 0), segments[x < 0]),
                    nochange = split(which(x == 0), segments[x == 0]))
        names(out[[1]]) <- NULL
        names(out[[2]]) <- NULL
        names(out[[3]]) <- NULL
        out
    }
    cid0 <- clusterMaker(chr, positions, maxGap = maxGap)
    segments <- getSegments3(x = x, cid = cid0, verbose = verbose)
    out <- vector("list", 2)
    names(out) <- c("up", "down")
    for(ii in 1:2) {
        idx <- segments[[ii]]
        if(length(idx) == 0) {
            out[[ii]] <- NULL
        } else {
            out[[ii]] <- data.frame(chr = sapply(idx, function(jj) chr[jj[1]]),
                                    start = sapply(idx, function(jj) min(positions[jj])),
                                    end = sapply(idx, function(jj) max(positions[jj])),
                                    idxStart = sapply(idx, function(jj) min(jj)),
                                    idxEnd = sapply(idx, function(jj) max(jj)),
                                    cluster = sapply(idx, function(jj) cid0[jj[1]]),
                                    n = sapply(idx, length))
        }
    }
    out
}

plotAnnoTrack <- function(gr, annoTrack) {
    ## check may need to be modified
    if(!all(sapply(annoTrack, function(xx) is(xx, "GRanges"))))
        stop("all elements in 'annoTrack' needs to be 'GRanges'")
    plot(start(gr), 1, type = "n", xaxt = "n", yaxt = "n", bty = "n",
         ylim = c(0.5, length(annoTrack) + 0.5), xlim = c(start(gr), end(gr)), xlab = "", ylab = "")
    lapply(seq(along = annoTrack), function(ii) {
        jj <- length(annoTrack) + 1- ii
        ir <- subsetByOverlaps(annoTrack[[ii]], gr)
        if(length(ir) > 0)  
            rect(start(ir)-0.5, jj - 0.15, end(ir), jj + 0.15, col = "grey60", border = NA)
        mtext(names(annoTrack)[ii], side = 2, at = jj, las = 1, line = 1)
        })

}

.bsPlotLines <- function(x, y, col, lty, lwd, plotRange) {
    if(sum(!is.na(y)) <= 1)
        return(NULL)
    xx <- seq(from = plotRange[1], to = plotRange[2], length.out = 500)
    yy <- approxfun(x, y)(xx)
    lines(xx, yy, col = col, lty = lty, lwd = lwd)
}

.bsPlotPoints <- function(x, y, z, col, pointsMinCov) {
    points(x[z>pointsMinCov], y[z>pointsMinCov], col = col, pch = 16, cex = 0.5)
}

.bsHighlightRegions <- function(regions, gr, ylim, regionCol, highlightMain) {
    if(is.data.frame(regions))
        regions <- data.frame2GRanges(regions)
    if(highlightMain)
        regions <- c(regions, gr)
    if(is.null(regions)) return(NULL)
    ## regions <- pintersect(region, rep(gr, length(regions)))
    ## regions <- regions[width(regions) == 0]
    regions <- subsetByOverlaps(regions, gr)
    regions <- pintersect(regions, rep(gr, length(regions)))
    if(length(regions) == 0)
        return(NULL)
    rect(xleft = start(regions), xright = end(regions), ybottom = ylim[1],
         ytop = ylim[2], col = regionCol, border = NA)
}

.bsGetCol <- function(object, col, lty, lwd) {
    ## Assumes that object has pData and sampleNames methods
    if(is.null(col)) {
        if("col" %in% names(pData(object)))
            col <- pData(object)[["col"]]
        else
            col <- rep("black", nrow(pData(object)))
    }
    if(length(col) != ncol(object))
        col <- rep(col, length.out = ncol(object))
    if(is.null(names(col)))
        names(col) <- sampleNames(object)
    
    if(is.null(lty)) {
        if("lty" %in% names(pData(object)))
            lty <- pData(object)[["lty"]]
        else
            lty <- rep(1, ncol(object))
    }
    if(length(lty) != ncol(object))
        lty <- rep(lty, length.out = ncol(object))
    if(is.null(names(lty)))
        names(lty) <- sampleNames(object)
    
    if(is.null(lwd)) {
        if("lwd" %in% names(pData(object)))
            lwd <- pData(object)[["lwd"]]
        else
            lwd <- rep(1, nrow(pData(object)))
    }
    if(length(lty) != ncol(object))
        lty <- rep(lty, length.out = ncol(object))
    if(is.null(names(lwd)))
        names(lwd) <- sampleNames(object)
                   
    return(list(col = col, lty = lty, lwd = lwd))
}

.bsPlotTitle <- function(gr, extend, main, mainWithWidth) {
    if(is.data.frame(gr))
        gr <- data.frame2GRanges(gr)
    if(length(gr) > 1) {
        warning("plotTitle: gr has more than one element")
        gr <- gr[1]
    }
    plotChr <- as.character(seqnames(gr))
    plotRange <- c(start(gr), end(gr))
    regionCoord <- sprintf("%s: %s - %s", plotChr, 
                           format(plotRange[1], big.mark = ",", scientific = FALSE),
                           format(plotRange[2], big.mark = ",", scientific = FALSE))
    if(mainWithWidth) {
        regionWidth <- sprintf("width = %s, extended = %s", 
                               format(width(gr), big.mark = ",", scientific = FALSE),
                               format(extend, big.mark = ",", scientific = FALSE))
        regionCoord <- sprintf("%s (%s)", regionCoord, regionWidth)
    }
    if(main != "") {
        main <- sprintf("%s\n%s", main, regionCoord)
    } else {
        main <- regionCoord
    }
    main
}

.bsGetGr <- function(object, region, extend) {
    if(is.null(region)) {
        gr <- GRanges(seqnames = seqnames(object)[1],
                      ranges = IRanges(start = min(start(object)),
                      end = max(start(object))))
    } else {
        if(is(region, "data.frame"))
            gr <- data.frame2GRanges(region, keepColumns = FALSE)
        else
            gr <- region
        if(!is(gr, "GRanges") || length(gr) != 1)
            stop("'region' needs to be either a 'data.frame' (with a single row) or a 'GRanges' (with a single element)")
        gr <- resize(gr, width = 2*extend + width(gr), fix = "center")
    }
    gr
}


plotRegion.no.smooth<-function (BSseq, region = NULL, extend = 0, main = "", addRegions = NULL,
    annoTrack = NULL, col = NULL, lty = NULL, lwd = NULL, BSseqTstat = NULL,
    stat = "tstat.corrected", stat.col = "black", stat.lwd = 1,
    stat.lty = 1, stat.ylim = c(-8, 8), mainWithWidth = TRUE,
    regionCol = alpha("red", 0.1), addTicks = TRUE, addPoints = FALSE,
    pointsMinCov = 5, highlightMain = FALSE)
{
    opar <- par(mar = c(0, 4.1, 0, 0), oma = c(5, 0, 4, 2), mfrow = c(1,
        1))
    on.exit(par(opar))
    if (is.null(BSseqTstat))
        layout(matrix(1:2, ncol = 1), heights = c(2, 1))
    else layout(matrix(1:3, ncol = 1), heights = c(2, 2, 1))
    .plotSmoothData(BSseq = BSseq, region = region, extend = extend,
        addRegions = addRegions, col = col, lty = lty, lwd = lwd,
        regionCol = regionCol, addTicks = addTicks, addPoints = addPoints,
        pointsMinCov = pointsMinCov, highlightMain = highlightMain)
    gr <- .bsGetGr(BSseq, region, extend)
    if (!is.null(BSseqTstat)) {
        BSseqTstat <- subsetByOverlaps(BSseqTstat, gr)
        plot(start(gr), 0.5, type = "n", xaxt = "n", yaxt = "n",
            ylim = stat.ylim, xlim = c(start(gr), end(gr)), xlab = "",
            ylab = "t-stat")
        axis(side = 2, at = c(-5, 0, 5))
        abline(h = 0, col = "grey60")
        mapply(function(stat, col, lty, lwd) {
            .bsPlotLines(start(BSseqTstat), BSseqTstat@stats[,
                stat], lty = lty, plotRange = c(start(gr), end(gr)),
                col = col, lwd = lwd)
        }, stat = stat, col = stat.col, lty = stat.lty, lwd = stat.lwd)
    }
    if (!is.null(annoTrack))
        plotAnnoTrack(gr, annoTrack)
    if (!is.null(main)) {
        main <- .bsPlotTitle(gr = region, extend = extend, main = main,
            mainWithWidth = mainWithWidth)
        mtext(side = 3, text = main, outer = TRUE, cex = 1)
    }
    return(invisible(NULL))
}


.plotSmoothData <- function(BSseq, region, extend, addRegions, col, lty, lwd, regionCol,
                            addTicks, addPoints, pointsMinCov, highlightMain) {
    gr <- .bsGetGr(BSseq, region, extend)
    BSseq <- subsetByOverlaps(BSseq, gr)
    
    ## Extract basic information
    sampleNames <- sampleNames(BSseq)
    names(sampleNames) <- sampleNames
    positions <- start(BSseq)
    smoothPs <- getMeth(BSseq, type = "raw")
    rawPs <- getMeth(BSseq, type = "raw")
    coverage <- getCoverage(BSseq)
        
    ## get col, lwd, lty
    colEtc <- bsseq:::.bsGetCol(object = BSseq, col = col, lty = lty, lwd = lwd)
    
    ## The actual plotting
    plot(positions[1], 0.5, type = "n", xaxt = "n", yaxt = "n",
         ylim = c(0,1), xlim = c(start(gr), end(gr)), xlab = "", ylab = "Methylation")
    axis(side = 2, at = c(0.2, 0.5, 0.8))
    if(addTicks)
        rug(positions)

    .bsHighlightRegions(regions = addRegions, gr = gr, ylim = c(0,1),
                        regionCol = regionCol, highlightMain = highlightMain)
    
    if(addPoints) {
        sapply(1:ncol(BSseq), function(sampIdx) {
            abline(v = positions[rawPs[, sampIdx] > 0.1], col = "grey80", lty = 1)
        })
    } # This adds vertical grey lines so we can see where points are plotted

    sapply(1:ncol(BSseq), function(sampIdx) {
        .bsPlotLines(positions, smoothPs[, sampIdx], col = colEtc$col[sampIdx],
                     lty = colEtc$lty[sampIdx], lwd = colEtc$lwd[sampIdx],
                     plotRange = c(start(gr), end(gr)))
    })

    if(addPoints) {
        sapply(1:ncol(BSseq), function(sampIdx) {
            .bsPlotPoints(positions, rawPs[, sampIdx], coverage[, sampIdx],
                          col = colEtc$col[sampIdx], pointsMinCov = pointsMinCov)
        })
    }
}

plotManyRegions.no.smooth <- function(BSseq, regions = NULL, extend = 0, main = "", addRegions = NULL,
                            annoTrack = NULL, col = NULL, lty = NULL, lwd = NULL, 
                            BSseqTstat = NULL, stat = "tstat.corrected", stat.col = "black",
                            stat.lwd = 1, stat.lty = 1, stat.ylim = c(-8,8),
                            mainWithWidth = TRUE, regionCol = alpha("red", 0.1), addTicks = TRUE,
                            addPoints = FALSE, pointsMinCov = 5, highlightMain = FALSE, verbose = TRUE) {
    cat("[plotManyRegions] preprocessing ...")
    if(!is.null(regions)) {
        if(is(regions, "data.frame"))
            gr <- data.frame2GRanges(regions, keepColumns = FALSE)
        else
            gr <- regions
        if(!is(gr, "GRanges"))
            stop("'regions' needs to be either a 'data.frame' (with a single row) or a 'GRanges' (with a single element)")
    } else {
        gr <- granges(BSseq)
    }
    gr <- resize(gr, width = 2*extend + width(gr), fix = "center")
    BSseq <- subsetByOverlaps(BSseq, gr)
    if(!is.null(BSseqTstat))
        BSseqTstat <- subsetByOverlaps(BSseqTstat, gr)
    
    if(length(start(BSseq)) == 0)
        stop("No overlap between BSseq data and regions")
    if(!is.null(main) && length(main) != length(gr))
        main <- rep(main, length = length(gr))
    cat("done\n")
    for(ii in seq(along = gr)) {
        if(verbose) cat(sprintf("[plotManyRegions]   plotting region %d (out of %d)\n", ii, nrow(regions)))
        plotRegion.no.smooth(BSseq = BSseq, region = regions[ii,], extend = extend,
                   col = col, lty = lty, lwd = lwd, main = main[ii], BSseqTstat = BSseqTstat,
                   stat = stat, stat.col = stat.col, stat.lwd = stat.lwd,
                   stat.lty = stat.lty, stat.ylim = stat.ylim,
                   addRegions = addRegions, regionCol = regionCol, mainWithWidth = mainWithWidth,
                   annoTrack = annoTrack, addTicks = addTicks, addPoints = addPoints,
                   pointsMinCov = pointsMinCov, highlightMain = highlightMain)
    }
}