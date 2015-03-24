source("setup.R")
pdf("unique_comparisons.pdf")
layout(matrix(1:3, ncol = 1))
for (i in 1:3) {
    plotVenn(uniquely_K4_membership[, c(i, 4)], plainNumeric, F)
    plotVenn(uniquely_K4_membership[, c(i, 5)], plainNumeric, F)
    plotVenn(uniquely_K4_membership[, c(i, 6)], plainNumeric, F)
}
for (i in 4:6) {
    plotVenn(uniquely_K4_membership[, c(1, i)], plainNumeric, F)
    plotVenn(uniquely_K4_membership[, c(2, i)], plainNumeric, F)
    plotVenn(uniquely_K4_membership[, c(3, i)], plainNumeric, F)
}

dev.off() 
