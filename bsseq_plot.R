.bsPlotPoints <- function(x, y, z, col, pointsMinCov) {
    points(x[z>pointsMinCov], y[z>pointsMinCov], col = col, pch = 16, cex = 0.5)
}

.bsPlotLines <- function(x, y, col, lty, lwd, plotRange) {
    if(sum(!is.na(y)) <= 1)
        return(NULL)
    xx <- seq(from = plotRange[1], to = plotRange[2], length.out = 500)
    yy <- approxfun(x, y)(xx)
    lines(xx, yy, col = col, lty = lty, lwd = lwd)
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
