###Remove sequence effect:
remove_gc_effect <- function(M, seq) {
    require(Biostrings)
    stopifnot(nrow(M)==length(seq))
    #dset = DNAStringSet(seq)
    freq = alphabetFrequency(seq)
    gccontent = rowSums(freq[,c("C", "G")]) / seq@ranges@width
    gc = round(gccontent*100); gc[gc<=20]<-20; gc[gc>=80]<-80
    
    sapply(1:ncol(M),function(i) {
      cat(".")
      meds=tapply(M[,i],gc,median)
      mads=tapply(M[,i],gc,mad)
      Index=match(gc,as.numeric(names(meds)))
      bias=meds[Index]
      sds=mads[Index]
      return((M[,i]-bias)/sds*mad(M[,i]))
    })
}
