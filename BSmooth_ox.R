makeClusters <- function(hasGRanges, maxGap = 10^8) {
    chrOrder <- as.character(runValue(seqnames(hasGRanges)))
    if(anyDuplicated(chrOrder))
        stop("argument 'hasGRanges' is not properly order")
    grBase <- granges(hasGRanges)
    clusters <- reduce(resize(grBase, width = 2*maxGap + 1, fix = "center"))
    start(clusters) <- pmax(rep(1, length(clusters)), start(clusters))
    clusters.sp <- split(clusters, seqnames(clusters))
    stopifnot(all(sapply(clusters.sp, function(cluster.gr) {
        if(length(cluster.gr) <= 1) return(TRUE)
        all(start(cluster.gr)[-length(cluster.gr)] < end(cluster.gr)[-1])
    }))) # are the clusters ordered within the chromosome? This is probably guranteed
    clusters <- Reduce(c, clusters.sp[chrOrder])
    stopifnot(all(chrOrder == runValue(seqnames(clusters))))
    ov <- findOverlaps(grBase, clusters)
    clusterIdx <- split(as.matrix(ov)[,1], as.matrix(ov)[,2])
    names(clusterIdx) <- NULL
    clusterIdx
}
    

BSmooth.ox <- function(BSseq, oxBSseq, hmC_level, cov_cutoff=10, ns = 70, h = 1000, maxGap = 10^8, parallelBy = c("sample", "chromosome"),
                    mc.preschedule = FALSE, mc.cores = 1, keep.se = FALSE, verbose = TRUE) {
		    require(locfit)
    smooth <- function(idxes, sampleIdx) {
	#print(idxes)
	#print(sampleIdx)
        ## Assuming that idxes is a set of indexes into the BSseq object
        ## sampleIdx is a single character
        this_sample_chr <- c(sampleNames(BSseq)[sampleIdx],
                             as.character(seqnames(BSseq)[idxes[1]]))
        if(verbose >= 2)
            cat(sprintf("[BSmooth]   smoothing start: sample:%s, chr:%s, nLoci:%s\n",
                        this_sample_chr[1], this_sample_chr[2], length(idxes)))
        #Cov <- getCoverage(BSseq, type = "Cov")[idxes, sampleIdx]
	Cov <- getCoverage(BSseq, type = "Cov")[idxes, sampleIdx]+getCoverage(oxBSseq, type = "Cov")[idxes, sampleIdx]
        
	M <- hmC_level[idxes, sampleIdx]
        pos <- start(BSseq)[idxes]
        stopifnot(all(diff(pos) > 0))
        wh <- which(Cov >= cov_cutoff & !is.na(M))
        nn <- ns / length(wh)
        if(length(wh) <= ns) {
            if(keep.se)
                se.coef <- rep(NA_real_, length(Cov))
            else
                se.coef <- NULL
            return(list(coef = rep(NA_real_, length(Cov)),
                        se.coef = se.coef,
                        trans = NULL, h = h, nn = nn))
        }
        sdata <- data.frame(pos = pos[wh],
                            M = M[wh],
                            Cov = Cov[wh])
        fit <- locfit(M ~ lp(pos, nn = nn, h = h), data = sdata,weights = Cov, maxk = 10000)
        pp <- preplot(fit, where = "data", band = "local",newdata = data.frame(pos = pos))
        if(keep.se) {
            se.coef <- pp$se.fit
        } else {
            se.coef <- NULL
        }
        if(verbose >= 2)
            cat(sprintf("[BSmooth]   smoothing end: sample:%s, chr:%s, nLoci:%s, nCoveredLoci:%s\n",
                        this_sample_chr[1], this_sample_chr[2], length(idxes), nrow(sdata)))
        return(list(coef = pp$fit, se.coef = se.coef,
                    trans = pp$trans, h = h, nn = nn))
    }
    stopifnot(class(BSseq) == "BSseq")
    parallelBy <- match.arg(parallelBy)
    if(verbose) cat("[BSmooth] preprocessing ... ")
    ptime1 <- proc.time()
    clusterIdx <- makeClusters(BSseq, maxGap = maxGap)
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if(verbose) cat(sprintf("done in %.1f sec\n", stime))

    sampleNames <- sampleNames(BSseq)
    names(sampleNames) <- sampleNames

    ptime.outer1 <- proc.time()
    switch(parallelBy, "sample" = {
        if(verbose) cat(sprintf("[BSmooth] smoothing by 'sample' (mc.cores = %d, mc.preschedule = %s)\n",
                                mc.cores, mc.preschedule))
        out <- mclapply(seq(along = sampleNames), function(sIdx) {
            ptime1 <- proc.time()
            tmp <- lapply(clusterIdx, function(jj) {
                try(smooth(idxes = jj, sampleIdx = sIdx))
            })
            coef <- do.call(c, lapply(tmp, function(xx) xx$coef))
            se.coef <- do.call(c, lapply(tmp, function(xx) xx$se.coef))
            ptime2 <- proc.time()
            stime <- (ptime2 - ptime1)[3]
            if(verbose) {
                cat(sprintf("[BSmooth] sample %s (out of %d), done in %.1f sec\n",
                            sampleNames[sIdx], length(sampleNames), stime))
            }
            return(list(coef = coef, se.coef = se.coef))
        }, mc.preschedule = mc.preschedule, mc.cores = mc.cores)
        if(any(sapply(out, is, class2 = "try-error")))
            stop("BSmooth encountered smoothing errors")
        coef <- do.call(cbind, lapply(out, function(xx) xx$coef))
        se.coef <- do.call(cbind, lapply(out, function(xx) xx$se.coef))
    }, "chromosome" = {
        if(verbose) cat(sprintf("[BSmooth] smoothing by 'chromosome' (mc.cores = %d, mc.preschedule = %s)\n",
                                mc.cores, mc.preschedule))
        out <- mclapply(1:length(clusterIdx), function(ii) {
            ptime1 <- proc.time()
            tmp <- lapply(seq(along = sampleNames), function(sIdx) {
                smooth(idxes = clusterIdx[[ii]], sampleIdx = sIdx)
            })
            coef <- do.call(cbind, lapply(tmp, function(xx) xx$coef))
            se.coef <- do.call(cbind, lapply(tmp, function(xx) xx$se.coef))
            ptime2 <- proc.time()
            stime <- (ptime2 - ptime1)[3]
            if(verbose)
                cat(sprintf("[BSmooth] chr idx %d (out of %d), done in %.1f sec\n",
                            ii, length(clusterIdx), stime))
            return(list(coef = coef, se.coef = se.coef))
        }, mc.preschedule = mc.preschedule, mc.cores = mc.cores)
        if(any(sapply(out, is, class2 = "try-error")))
            stop("BSmooth encountered smoothing errors")
        coef <- do.call(rbind, lapply(out, function(xx) xx$coef))
        se.coef <- do.call(rbind, lapply(out, function(xx) xx$se.coef))
    })
    ptime.outer2 <- proc.time()
    stime.outer <- (ptime.outer2 - ptime.outer1)[3]
    if(verbose)
        cat(sprintf("[BSmooth] smoothing done in %.1f sec\n", stime.outer))
    
    rownames(coef) <- NULL
    colnames(coef) <- sampleNames(BSseq)
    if(!is.null(se.coef)) {
        rownames(se.coef) <- NULL
        colnames(se.coef) <- sampleNames(BSseq)
    }

    if(!is.null(coef))
        assay(BSseq, "coef") <- coef
    if(!is.null(se.coef))
        assay(BSseq, "se.coef") <- se.coef
    mytrans <- function(x) {
	x
    }
    environment(mytrans) <- baseenv()
    BSseq@trans <- mytrans
    parameters <- list(smoothText = sprintf("BSmooth (ns = %d, h = %d, maxGap = %d)", ns, h, maxGap),
                       ns = ns, h = h, maxGap = maxGap)
    BSseq@parameters <- parameters
    #BSseq@hmc<- hmC_level
    BSseq
}


BSmooth.ox.tstat <- function(BSseq, hmC_level,group1, group2, estimate.var = c("same", "paired", "group2"),
                          local.correct = TRUE, maxGap = NULL, qSd = 0.75, k = 101, mc.cores = 1, verbose = TRUE){
			  require(locfit)
    smoothSd <- function(Sds, k) {
        k0 <- floor(k/2)
        if(all(is.na(Sds))) return(Sds)
        thresSD <- pmax(Sds, quantile(Sds, qSd, na.rm = TRUE), na.rm = TRUE)
        addSD <- rep(median(Sds, na.rm = TRUE), k0)
        sSds <- as.vector(runmean(Rle(c(addSD, thresSD, addSD)), k = k))
        sSds
    }
    compute.correction <- function(idx, qSd = 0.75) {
        xx <- start(BSseq)[idx]
        yy <- tstat[idx]
        suppressWarnings({
            drange <- diff(range(xx, na.rm = TRUE))
        })
        if(drange <= 25000)
            return(yy)
        tstat.function <- approxfun(xx, yy)
        xx.reg <- seq(from = min(xx), to = max(xx), by = 2000)
        yy.reg <- tstat.function(xx.reg)
        fit <- locfit(yy.reg ~ lp(xx.reg, h = 25000, deg = 2, nn = 0),
                      family = "huber", maxk = 50000) 
        correction <- predict(fit, newdata = data.frame(xx.reg = xx))
        yy - correction 
    }

    estimate.var <- match.arg(estimate.var)
    stopifnot(is(BSseq, "BSseq"))
    stopifnot(hasBeenSmoothed(BSseq))
    if(is.character(group1)) {
        stopifnot(all(group1 %in% sampleNames(BSseq)))
        group1 <- match(group1, sampleNames(BSseq))
    }
    if(is.numeric(group1)) {
        stopifnot(min(group1) >= 1 & max(group1) <= ncol(BSseq))
    } else stop("problems with argument 'group1'")
    if(is.character(group2)) {
        stopifnot(all(group2 %in% sampleNames(BSseq)))
        group2 <- match(group2, sampleNames(BSseq))
    }    
    if(is.numeric(group2)) {
        stopifnot(min(group2) >= 1 & max(group2) <= ncol(BSseq))
    } else stop("problems with argument 'group2'")
    stopifnot(length(intersect(group1, group2)) == 0)
    stopifnot(length(group1) > 0)
    stopifnot(length(group2) > 0)
    stopifnot(length(group1) + length(group2) >= 3)
    if(estimate.var == "paired")
        stopifnot(length(group1) == length(group2))
    
    if(any(rowSums(getCoverage(BSseq)[, c(group1, group2)]) == 0))
        warning("Computing t-statistics at locations where there is no data; consider subsetting the 'BSseq' object first")
    
    if(is.null(maxGap))
        maxGap <- BSseq@parameters$maxGap
    if(is.null(maxGap))
        stop("need to set argument 'maxGap'")
    
    if(verbose) cat("[BSmooth.tstat] preprocessing ... ")
    ptime1 <- proc.time()
    clusterIdx <- makeClusters(BSseq, maxGap = maxGap)
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if(verbose) cat(sprintf("done in %.1f sec\n", stime))
        
    if(verbose) cat("[BSmooth.tstat] computing stats within groups ... ")
    ptime1 <- proc.time()
    #allPs <- getMeth(BSseq, type = "smooth", what = "perBase",confint = FALSE)
     allPs <- hmC_level
    group1.means <- rowMeans(allPs[, group1, drop = FALSE], na.rm = TRUE)
    group2.means <- rowMeans(allPs[, group2, drop = FALSE], na.rm = TRUE)

    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if(verbose) cat(sprintf("done in %.1f sec\n", stime))
    
    if(verbose) cat("[BSmooth.tstat] computing stats across groups ... ")
    ptime1 <- proc.time()
    switch(estimate.var,
           "group2" = {
               rawSds <- rowSds(allPs[, group2, drop = FALSE], na.rm = TRUE)
               smoothSds <- do.call(c, mclapply(clusterIdx, function(idx) {
                   smoothSd(rawSds[idx], k = k)
               }, mc.cores = mc.cores))
               scale <- sqrt(1/length(group1) + 1/length(group2))
               tstat.sd <- smoothSds * scale
           },
           "same" = {
               rawSds <- sqrt( ((length(group1) - 1) * rowVars(allPs[, group1, drop = FALSE]) +
                                (length(group2) - 1) * rowVars(allPs[, group2, drop = FALSE])) /
                              (length(group1) + length(group2) - 2))
               smoothSds <- do.call(c, mclapply(clusterIdx, function(idx) {
                   smoothSd(rawSds[idx], k = k)
               }, mc.cores = mc.cores))
               scale <- sqrt(1/length(group1) + 1/length(group2))
               tstat.sd <- smoothSds * scale
           },
           "paired" = {
               rawSds <- rowSds(allPs[, group1, drop = FALSE] - allPs[, group2, drop = FALSE])
               smoothSds <- do.call(c, mclapply(clusterIdx, function(idx) {
                   smoothSd(rawSds[idx], k = k)
               }, mc.cores = mc.cores))
               scale <- sqrt(1/length(group1))
               tstat.sd <- smoothSds * scale
           })
    tstat <- (group1.means - group2.means) / tstat.sd
    is.na(tstat)[tstat.sd == 0] <- TRUE
    if(local.correct) {
        tstat.corrected <- do.call(c, mclapply(clusterIdx,
                                               compute.correction, qSd = qSd,
                                               mc.cores = mc.cores))
    }
    ptime2 <- proc.time()
    stime <- (ptime2 - ptime1)[3]
    if(verbose) cat(sprintf("done in %.1f sec\n", stime))
    
    if(local.correct) {
        stats <- cbind(rawSds, tstat.sd, group2.means, group1.means,
                       tstat, tstat.corrected)
        colnames(stats) <- c("rawSds", "tstat.sd",
                             "group2.means", "group1.means", "tstat",
                             "tstat.corrected")
 
    } else {
        stats <- cbind(rawSds, tstat.sd, group2.means, group1.means,
                       tstat)
        colnames(stats) <- c("rawSds", "tstat.sd",
                             "group2.means", "group1.means", "tstat")
    }
    
    parameters <- c(BSseq@parameters,
                    list(tstatText = sprintf("BSmooth.tstat (local.correct = %s, maxGap = %d)",
                         local.correct, maxGap),
                         group1 = group1, group2 = group2, k = k, qSd = qSd,
                         local.correct = local.correct, maxGap = maxGap))
    out <- BSseqTstat(gr = granges(BSseq), stats = stats, parameters = parameters)
    out
}
