##
## read the segregation data from a txt or txt.gz file
## chr    pos      id          genotypes         MAF               noCalls caseString    controlString statistic           p
## 6      25726329 rs112943240 C/C|C/CTT|CTT/CTT 0.307875605815832 0       2330|2126|513 3012|2676|557 -1.4112760160634399 0.1581632435465825
## chr  pos     id          genotypes               caseString    controlstring  noCalls  statistic          p   
## 1	866511	rs60722469  CCCCT/CCCCT|C/CCCCT|C/C 2450|2115|401 3059|2670|513  6        0.3458728500395499 0.7294382733518878
read.seg = function(file="seg.txt.gz") {
    seg = read.table(file=file, header=F, sep="\t")
    ## new
    ## colnames(seg) = c("chr","pos","id","genotypes","MAF","noCalls","caseString","controlString","statistic","p")
    ## old
    colnames(seg) = c("chr","pos","id","genotypes","caseString","controlString","noCalls","statistic","p")
    ## store the number of genotypes, 2=SNP without ALT HOM so we can compute odds ratio
    seg$ngenotypes = 0
    seg$mcc = 0
    foo = strsplit(seg$genotypes, "|", TRUE)
    for (i in 1:nrow(seg)) {
        seg$ngenotypes[i] = length(foo[[i]])
        ## calculate a very faux MCC based on the first two or three genotypes 
        ## caseString    controlString
        ## 2450|2115|401 3059|2670|513
        if (seg$ngenotypes[i]==2) {
            caseBar = strsplit(seg$caseString[i], "|", TRUE)
            controlBar = strsplit(seg$controlString[i], "|", TRUE)
            caseREF = as.numeric(caseBar[[1]][1])
            controlREF = as.numeric(controlBar[[1]][1])
            caseHET = as.numeric(caseBar[[1]][2])
            controlHET = as.numeric(controlBar[[1]][2])
            TP = caseHET
            FP = controlHET
            TN = controlREF
            FN = caseREF
            seg$mcc[i] = (TP*TN-FP*FN) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
        } else if (seg$ngenotypes[i]==3) {
            caseBar = strsplit(seg$caseString[i], "|", TRUE)
            caseREF = as.numeric(caseBar[[1]][1])
            caseHET = as.numeric(caseBar[[1]][2])
            caseHOM = as.numeric(caseBar[[1]][3])
            controlBar = strsplit(seg$controlString[i], "|", TRUE)
            controlREF = as.numeric(controlBar[[1]][1])
            controlHET = as.numeric(controlBar[[1]][2])
            controlHOM = as.numeric(controlBar[[1]][3])
            TP = caseHET + caseHOM
            FP = controlHET + controlHOM
            TN = controlREF
            FN = caseREF
            seg$mcc[i] = (TP*TN-FP*FN) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
        }
    }
    ## store -log10(p) to save time later
    seg$mlog10p = -log10(seg$p)
    ## get the chromosomes into a simple list
    chrs = unique(seg$chr)
    return(seg)
}


