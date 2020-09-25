show.sig.genes = function(seg, pmax, genes) {
    sigseg = seg[seg$p<pmax & !is.na(seg$p),]
    nsig = nrow(sigseg)
    for (i in 1:nsig) {
        inside = genes$seqid==sigseg$chr[i] & genes$start<=sigseg$pos[i] & genes$end>=sigseg$pos[i]
        print(paste(sigseg$chr[i], sigseg$id[i], sigseg$pos[i], sigseg$p[i], genes$name[inside], genes$seqid[inside], genes$start[inside], genes$end[inside]), quote=F)
    }
}
