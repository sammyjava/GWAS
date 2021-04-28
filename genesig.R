## count the significant loci per gene
##                 seqid  start    end strand length          name
## ENSG00000186092     1  69091  70008      +    917         OR4F5
##
## colnames(seg) = c("chr","pos","id","ref","alts","caseString","controlString","caseNoCalls","controlNoCalls","statistic","p")

genesig = function(seg) {
    ## scan through genes DF
    for (i in 1:nrow(genes)) {
        loci = seg[seg$chr==genes$seqid[i] & seg$pos>=genes$start[i] & seg$pos<=genes$end[i] & !is.nan(seg$p),]
        nsig = nrow(loci[loci$p<5e-8,])
        nhigh = nrow(loci[loci$p<1e-2,])
        if (nhigh>10) {
            print(paste(genes$name[i],genes$seqid[i],":",genes$start[i],"-",genes$end[i],nrow(loci),nsig,nhigh))
        }
    }
}

    
    
