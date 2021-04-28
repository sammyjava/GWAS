source("snpData.R")
##
## plot log Fisher p value of each seg call on a single chromosome in the given start,end range
## we focus on loci with positive odds ratio since we're interested in ALTs that lead to the condition
##
## http://myvariant.info/v1/query?q=rs727503873
##
## colnames(seg) = c("chr","pos","id","ref","alts","caseString","controlString","caseNoCalls","controlNoCalls","statistic","p")
##
plot.p.region = function(seg=seg, chr="0",
                         start=0, end=0, idStart=NULL, idEnd=NULL, geneStart=NULL, geneEnd=NULL,
                         gene=NULL, title=NULL,
                         label=FALSE, labelSNP=FALSE, showGenes=FALSE, showDensity=FALSE,
                         ymin=0, ymax=0, pSig=5e-8, minMAF=0.0,
                         posList=c()) {

    sigColor = "red"
    listColor = "blue"

    if (!is.null(gene)) {
        geneRecord = genes[genes$name==gene,]
        chr = geneRecord$seqid
        start = geneRecord$start
        end = geneRecord$end
    }

    if (chr=="0") {
        chr = seg$chr[1]
        if (start==0) {
            start = min(seg$pos)
        }
        if (end==0) {
            end = max(seg$pos)
        }
    } else if (!is.null(idStart) && !is.null(idEnd)) {
        start = seg$pos[seg$id==idStart]
        end = seg$pos[seg$id==idEnd]
    } else if (!is.null(geneStart) && !is.null(geneEnd)) {
        start = genes$start[genes$name==geneStart]
        end = genes$end[genes$name==geneEnd]
    } else {
        if (start==0) {
            start = min(seg$pos[seg$chr==chr])
        } else {
            start = min(seg$pos[seg$chr==chr & seg$pos>=start])
        }
        if (end==0) {
            end = max(seg$pos[seg$chr==chr])
        } else if (is.null(gene)) {
            end = max(seg$pos[seg$chr==chr & seg$pos<=end])
        }
    }
    
    ## collect the points of interest
    pts = seg$pos>0
    if (chr!="0") {
        pts = pts & seg$chr==chr & seg$pos>=start & seg$pos<=end
    }
    if (minMAF>0.0) {
        pts = pts & seg$MAF>minMAF
    }

    numLoci = length(seg$pos[pts])

    ## significant points
    ptsSig = pts & seg$p<pSig & !is.na(seg$p)
    hasSig = length(seg$p[ptsSig]) > 0
    ptsSigNum = length(seg$pos[ptsSig])
    ptsHigh = pts & seg$p<1e-2 & !is.na(seg$p)
    ptsHighNum = length(seg$pos[ptsHigh])

    ## infinite points
    ptsInf = is.infinite(seg$mlog10p)
    hasInf = length(seg$p[ptsInf]) > 0

    ## limits
    xlim = c(start, end)
    if (ymax==0) ymax = max(seg$mlog10p[pts & is.finite(seg$mlog10p)])
    if (hasInf) {
        ymax = ymax + 1
        seg$mlog10p[ptsInf] = ymax
    }
    ylim = c(ymin, ymax)

    ## plot title
    if (is.null(gene) && is.null(title)) {
        title = paste(chr,":",start,"-",end, numLoci,"loci","MAF >",minMAF,"N(p<5e-8) =",ptsSigNum,"N(p<0.01) =",ptsHighNum)
    } else if (is.null(title)) {
        title = paste(gene, numLoci,"loci","MAF >",minMAF,"N(p<5e-8) =",ptsSigNum,"N(p<0.01) =",ptsHighNum)
    }

    ## plot
    plot(seg$pos[pts], seg$mlog10p[pts],
         xlab=paste("Chr",chr,"position"),
         ylab="-log10(p)",
         xlim=xlim,
         ylim=ylim,
         pch=1, cex=0.8, col="black",
         ## main=paste(deparse(substitute(seg)),chr,":",start,"-",end)
         main = title
         )
    ## use different symbol for infinite points
    points(seg$pos[ptsInf], seg$mlog10p[ptsInf], pch=17, cex=1.5, col=sigColor)

    ## ## vertical chromosome lines if plotting full genome
    ## if (chr=="0") {
    ##     segpts = seg[pts,]
    ##     currentChr = "0"
    ##     for (i in 1:nrow(segpts)) {
    ##         if (segpts$chr[i]!=currentChr) {
    ##             currentChr = segpts$chr[i]
    ##             lines(c(i,i), ylim, lwd=1, col="gray")
    ##             text(i, ylim[2], currentChr, cex=0.5)
    ##         }
    ##     }
    ## }

    ## highlight somewhat significant p values
    points(seg$pos[ptsHigh], seg$mlog10p[ptsHigh], pch=1, cex=0.8, col=sigColor)

    ## highlight highly significant p values
    points(seg$pos[ptsSig], seg$mlog10p[ptsSig], pch=19, cex=0.9, col=sigColor)

    ## highlight loci in the posList
    ptsPosList = seg$pos %in% posList
    points(seg$pos[ptsPosList], seg$mlog10p[ptsPosList], pch=19, cex=0.9, col=listColor)
    text(seg$pos[ptsPosList], seg$mlog10p[ptsPosList], seg$id[ptsPosList], col=listColor, pos=4, cex=0.8, offset=0.3)

    ## line at pSig
    ## lines(xlim, -log10(c(pSig,pSig)), col="gray", lty=2)

    ## label significant points if requested
    if (hasSig && label) {
        ## snpInfo = c()
        ## for (rsId in seg$id[ptsSig]) {
        ##     clinvar.clinsig = ""
        ##     if (labelSNP && startsWith(rsId,"rs")) {
        ##         ## query the SNP API
        ##         json = snpData(rsId)
        ##         for (hit in json$hits) {
        ##             if ("clinvar" %in% colnames(hit)) {
        ##                 clinvar = hit$clinvar
        ##                 for (clinsig in clinvar$clinsig) {
        ##                     if (!is.na(clinsig)) clinvar.clinsig = paste(clinvar.clinsig, clinsig)
        ##                 }
        ##             }
        ##         }
        ##     }
        ##     snpInfo = c(snpInfo, clinvar.clinsig)
        ## }

        ## control first since genotypes has REF first
        ## textStr = paste(seg$pos[ptsSig], seg$id[ptsSig], seg$genotypes[ptsSig], seg$controlString[ptsSig], seg$caseString[ptsSig])
        textStr = seg$id[ptsSig]
        ## ## add odds ratio if just two genotypes
        ## for (i in 1:length(textStr)) {
        ##     if (seg$ngenotypes[ptsSig][i]==2) {
        ##         caseCounts = as.numeric(strsplit(schiz$caseString[ptsSig][i], "|", TRUE)[[1]])
        ##         controlCounts = as.numeric(strsplit(schiz$controlString[ptsSig][i], "|", TRUE)[[1]])
        ##         caseRatio = caseCounts[2]/caseCounts[1]
        ##         controlRatio = controlCounts[2]/controlCounts[1]
        ##         oddsRatio = caseRatio/controlRatio
        ##         textStr[i] = paste(textStr[i], round(oddsRatio,2))
        ##     }
        ## }

        ## if (length(snpInfo)>0) textStr = paste(textStr, snpInfo)
        
        text(seg$pos[ptsSig], seg$mlog10p[ptsSig], textStr, col=sigColor, pos=4, cex=0.8, offset=0.3)
    }

    ## show gene or genes if requested
    ## REQUIRES load-genes!!
    ypos = par("yaxp")[2]
    bar = c(ypos-(ymax-ymin)*0.01, ypos+(ymax-ymin)*0.01)
    ## if (!is.null(gene)) {
    ##     lines(xlim, c(ypos,ypos))
    ##     x = (start+end)/2
    ##     text(x, ypos, gene, pos=1, cex=0.7, offset=0.4)
    ##     ## end bars
    ##     lines(rep(start,2), bar)
    ##     lines(rep(end,2), bar)
    ##     if (geneRecord$strand=="-") {
    ##         text((start+end)/2, ypos, "<")
    ##     } else {
    ##         text((start+end)/2, ypos, ">")
    ##     }
    ## }
    if (showGenes) {
        within = genes$seqid==chr & genes$end>=start & genes$start<=end
        genesWithin = genes[within,]
        for (i in 1:nrow(genesWithin)) {
            lines(c(genesWithin$start[i],genesWithin$end[i]), c(ypos,ypos))
            x = (genesWithin$start[i]+genesWithin$end[i])/2
            text(x, ypos, genesWithin$name[i], pos=1, cex=0.6)
            lines(rep(genesWithin$start[i],2), bar)
            lines(rep(genesWithin$end[i],2), bar)
            if (genesWithin$strand[i]=="-") {
                text((genesWithin$start[i]+genesWithin$end[i])/2, ypos, "<")
            } else {
                text((genesWithin$start[i]+genesWithin$end[i])/2, ypos, ">")
            }
        }
    }

    if (showDensity) {
        width = 1000000
        for (i in seq(start,end,by=width)) {
            j = i + width
            n = length(seg$pos[seg$chr==chr & seg$pos>=i & seg$pos<j])
            if (n>0) {
                lines(c(i,j), log10(c(n,n)), col="blue", lwd=3)
                text((i+j)/2, log10(n), n, col="blue", pos=3)
            }
            
        }
    }
}

