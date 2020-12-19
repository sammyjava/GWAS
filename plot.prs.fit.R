library(fitdistrplus)

################################################################################################
## plot a histogram and normal distribution fits of the polygenic risk scores of the individuals
##
## input: dataframe of "sample", "label", "score"
##        label ("case" or "ctrl")
################################################################################################

plot.prs.fit = function(prs, label) {
    ## just in case
    colnames(prs) = c("sample", "label", "score")
    ## use method `plot.fitdist`
    fit = fitdist(prs$score[prs$label=="ctrl"], "norm")
    plot(fit)
}
