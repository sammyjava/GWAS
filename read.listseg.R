##
## read the segregation data from a listseg.txt file
##
## Contig  ID                  HET  CaseREF  ControlREF  CaseHET  ControlHET  CaseHOM  ControlHOM  CaseNC  ControlNC  StdStat             p                    OR
## 6       AA_A_9_30018537_FS  AP   246      219         147      152         7        29          0       0          1.9349060283303026  0.05300182831300198  0.682520325203252

read.listseg = function(file="listseg.txt.gz") {
    seg = read.table(file=file, header=T, sep="\t")
    seg$Total = seg$CaseREF + seg$ControlREF + seg$CaseHET + seg$ControlHET + seg$CaseHOM + seg$ControlHOM
    seg$mlog10p = -log10(seg$p)
    return(seg)
}


