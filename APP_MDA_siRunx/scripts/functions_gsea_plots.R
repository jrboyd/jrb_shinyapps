

gsea.plotVennDiagrams = function() {
    # plot venn diagrams of +/- pairs of marks
    
    
    
    keep = apply(cbind(k27me3_tc, k4me3_tc), 1, function(x) {
        return(max(x) > log2(1.1))
    })
    
    k27me3_tc_bool = k27me3_tc[keep, ] > 0.3
    k4me3_tc_bool = k4me3_tc[keep, ] > 1
    k4ac_tc_bool = k4ac_tc[keep, ] > 1
    
    keep = apply(k27me3_tc_bool, 1, any) | apply(k4me3_tc_bool, 1, any)
    
    
    # layout(matrix(1:3, ncol = 3))
    pdf("bivalent_venn.pdf", width = 7)
    vennDiagram(cbind(k27me3_tc_bool[, 1], k4me3_tc_bool[, 1])[keep, ], names = c("H3K27me3", "H3K4me3"))
    title("MCF10A")
    vennDiagram(cbind(k27me3_tc_bool[, 2], k4me3_tc_bool[, 2])[keep, ], names = c("H3K27me3", "H3K4me3"))
    title("MCF7")
    vennDiagram(cbind(k27me3_tc_bool[, 3], k4me3_tc_bool[, 3])[keep, ], names = c("H3K27me3", "H3K4me3"))
    title("MDA231")
    dev.off()
}

# par(bg = rgb(0,0,0,0)) evaluates subsets of genes that are marked and not expressed. vice versa
gsea.falsePosNeg = function() {
    res = list()
    pdf("false_posneg.pdf")
    for (l in lines) {
        
        
        exprThresh = 1
        markThresh = 1
        falsePos4me3 = getExpr(l) < exprThresh & get4me3(l) > markThresh  #no expr with 4me3+
        falseNeg4me3 = getExpr(l) > exprThresh & get4me3(l) < markThresh  #+ expr with 4me3-
        falsePos4ac = getExpr(l) < exprThresh & get4ac(l) > markThresh  #no expr with 4me3+
        falseNeg4ac = getExpr(l) > exprThresh & get4ac(l) < markThresh  #+ expr with 4me3-
        bothFalsePos = falsePos4me3 & falsePos4ac
        bothFalseNeg = falseNeg4me3 & falseNeg4ac
        
        
        
        a = cbind(falsePos4me3, falseNeg4me3, falsePos4ac, falseNeg4ac)
        b = cbind(a, bothFalsePos, bothFalseNeg)
        
        res[[length(res) + 1]] = colSums(b)
        
        layout(mat = matrix(1:3, nrow = 3, ncol = 1), heights = c(0.2, 0.4, 0.4))
        # textplot(l,cex = 10)
        plot(1, 1)
        vennDiagram(a[, c(1, 3)])
        mtext("marked but not expressed", side = 3, adj = 0.5, padj = 0.9)
        vennDiagram(a[, c(2, 4)])
        mtext("expressed but not marked", side = 3, adj = 0.5, padj = 0.9)
        layout(1)
        vennDiagram(a)
        minx = par("usr")[1]
        maxx = par("usr")[2]
        text(l, x = (maxx - minx)/2 + minx, y = (maxx - minx) * 7/8 + minx)
    }
    dev.off()
    
    resMatrix = matrix(unlist(res), ncol = 3)
    rownames(resMatrix) = names(res[[1]])
    colnames(resMatrix) = lines
    write.csv(resMatrix, "false_posneg.csv")
}



gsea.plotVolcano = function(ensg_inPath, pdfName = "volcano.pdf", startCol = 1, endCol = 3, minFE = 0, minFC = 0, 
    filterSubset = T, doPlots = T, rectangles = NULL) {
    print(length(ensg_inPath))
    bothMarks = getStandardTable()
    
    bothMarks_inPath = getIntersect(bothMarks, ensg_inPath)
    magnitudeMetric = function(x) {
        # out = apply(x, 1, max) #row max of x
        out = apply(x, 1, min)  #row min of x
        # out = log2(rowMeans(2^bothMarks_inPath[,c(1,2)])) #row mean of x
        return(out)
    }
    
    get_x = function(dat, a, b) {
        return(dat[, a] - dat[, b])
    }
    
    get_y = function(dat, a, b) {
        return(magnitudeMetric(dat[, c(a, b)]))
    }
    
    bg_plot = function(dat, a, b) {
        # draw background with alpha value
        plot(get_x(dat, a, b), get_y(dat, a, b), pch = 16, col = rgb(0, 0, 0, 0.1), xlab = paste("log2 FC", eName, 
            "from", sName), ylab = paste("Minimum log2 FE of", sName, "and", eName, "over input"),
            xlim = c(-4,4), ylim = c(0,6))
    }
    
    fg_plot = function(dat, a, b) {
        # draw foreground points over background
        points(get_x(dat, a, b), get_y(dat, a, b), col = "red")
    }
    
    
    
    start = startCol
    end = endCol
    
    sName = colnames(bothMarks)[start]
    eName = colnames(bothMarks)[end]
    
    bg_dat = bothMarks
    fg_dat = bothMarks_inPath
    
    
    if (filterSubset) {
        keep = apply(bg_dat[, c(start, end)], 1, max) > minFE
        bg_dat = bg_dat[keep, ]
        
        keep = apply(fg_dat[, c(start, end)], 1, max) > minFE
        fg_dat = fg_dat[keep, ]
    }
    
    if (doPlots) {
#         bg_plot(bg_dat, end, start)
#         
#         bg_plot(bg_dat, end, start)
#         fg_plot(fg_dat, end, start)
        
        bg_plot(bg_dat, end, start)
        fg_plot(fg_dat, end, start)
    }
    
    
    # return(list(res_bg, res_fg))
}

gsea.plotVolcanoRandom = function(pdfName = "volcanoRandom.pdf", startCol = 1, endCol = 3, minFE = 0.5, minFC = 0.7, 
    filterSubset = T) {
    bothMarks = getStandardTable()
    if (filterSubset) {
        start = startCol
        end = endCol
        keep = apply(bothMarks[, c(start, end)], 1, max) > minFE
        bothMarks = bothMarks[keep, ]
    }
    numTrials = nrow(getIntersect(bothMarks, fullList))
    # randi = runif(nrow(bothMarks)) o = order(randi) randList = rownames(bothMarks)[o[1:numTrials]]
    randList = getRandSubsetList(numTrials, bothMarks)
    
    return(gsea.plotVolcano(randList, pdfName, startCol, endCol, minFE, minFC, filterSubset))
}

gsea.plotVSbg = function(ensg_inPath, pdfName = "triplot.pdf") {
    pdf(pdfName, width = 12, height = 6)
    # png('triplot.png', width = 1600, height = 1200)
    layout(matrix(c(0, 4:6, 0, 7, 1:3, 0, rep(8, 5)), ncol = 5, nrow = 3, byrow = T), widths = c(0.3, 1, 1, 1, 
        0.1), heights = c(0.3, 1, 0.1))
    par(mai = rep(0, 4))
    cex.axis = 3.2
    bothMarks = cbind(k4ac_tc, k4me3_tc)
    keep = apply(bothMarks, 1, max) > 1
    for (l in lines) {
        d = getLine(l)
        
        d_inPath = d[ensg_inPath, ]
        xi = ncol(d)
        yi = ncol(d) - 1
        
        x_name = colnames(d)[xi]
        y_name = colnames(d)[yi]
        
        # keep = apply(d[,c(xi,yi)], 1, function(x){return(max(x) > 1)})
        d = d[keep, ]
        
        x_dat = d[, xi]
        names(x_dat) = rownames(d)
        y_dat = d[, yi]
        names(y_dat) = rownames(d)
        
        plot(x_dat, y_dat, pch = 16, col = rgb(0, 0, 0, 0.2), xlim = c(0, 6), ylim = c(0, 6), axes = T, xaxt = "n", 
            yaxt = "n", xlab = "", ylab = "")
        axis(side = 1, cex.axis = cex.axis, padj = 0.5)
        if (l == "MCF10A") 
            axis(side = 2, cex.axis = cex.axis, padj = -0.3)
        points(x_dat[ensg_inPath], y_dat[ensg_inPath], pch = 16, col = rgb(1, 0, 0, 1))
    }
    # for(l in lines){ d = getLine(l) d_inPath = d[ensg_inPath,] xi = ncol(d) - 2 yi = ncol(d) - 1
    # plot(d[,c(xi,yi)], pch = 16, col = rgb(0,0,0,.2), xlim = c(0,2.8), ylim = c(0,6), axes = T, xaxt = 'n',
    # yaxt = 'n', xlab = '', ylab = '') axis(side = 1, cex.axis = cex.axis, padj = .5) if(l == 'MCF10A')
    # axis(side = 2, cex.axis = cex.axis, padj = -.3) points(d_inPath[,c(xi,yi)], pch = 16, col = rgb(1,0,0,1)) }
    # textplot(object = 'MCF10A') textplot(object = 'MCF7') textplot(object = 'MDA231') textplot(object = 'K4ac
    # vs K4me3') textplot(object = 'H3K4me3') textplot(object = 'H3K4ac')
    
    dev.off()
}

gsea.plotBoxplots = function(ensg_inPath, pdfName = "boxplots.pdf") {
    r = 0.95
    
    
    
    bothMarks = getStandardTable()
    pdf(pdfName, width = 12, height = 12)
    par(bg = rgb(0, 0, 0, 0))
    par(mai = c(0, 0, 0, 0))
    layout(matrix(1:4, nrow = 4))
    for (l in lines) {
        d = getLine(l)
        
        d_inPath = d[ensg_inPath, ]
        
        
        xi = ncol(d)
        yi = ncol(d) - 1
        
        x_name = colnames(d)[xi]
        y_name = colnames(d)[yi]
        
        # keep = apply(d[,c(xi,yi)], 1, function(x){return(max(x) > 1)})
        d = bothMarks
        ensg_rand = getRandSubsetList(nrow(d), d)
        
        x_dat = d[, xi]
        names(x_dat) = rownames(d)
        y_dat = d[, yi]
        names(y_dat) = rownames(d)
        
        
        boxplot(list(x_dat[ensg_rand], x_dat[ensg_inPath], y_dat[ensg_rand], y_dat[ensg_inPath]), at = 1:4, names = rep("", 
            4), range = r, boxwex = 0.5)
        if (l == "MDA231") 
            axis(side = 1, labels = c(x_name, paste(x_name, "subset"), y_name, paste(y_name, "subset")), at = 1:4, 
                las = 2)
        # textplot(object = 'H3K4me3') textplot(object = 'H3K4ac')
        
    }
    
    
    
    bothMarks = getStandardTable()
    bothMarks_inPath = bothMarks[intersect(ensg_inPath, rownames(bothMarks)), ]
    
    boxplot.cellLines = function(index, drawLabels = F) {
        a = bothMarks
        a[, 1:3] = a[, 1:3] - a[, index]
        a[, 4:6] = a[, 4:6] - a[, index + 3]
        a_inPath = a[intersect(ensg_inPath, rownames(a)), ]
        a = a[intersect(ensg_rand, rownames(a)), ]
        boxplot(a, at = 1:6, names = rep("", 6), range = r)
        if (drawLabels) 
            axis(side = 1, labels = colnames(a), at = 1:6, las = 2)
        boxplot(a_inPath, at = 1:6, names = rep("", 6), range = r)
        if (drawLabels) 
            axis(side = 1, labels = paste(colnames(a), "subset"), at = 1:6, las = 2)
    }
    
    layout(matrix(c(1:6, 7, 7), nrow = 4, ncol = 2, byrow = T))
    par(mai = c(0, 0, 0, 0))
    for (i in 1:3) {
        boxplot.cellLines(i, drawLabels = i == 3)
    }
    
    # layout(matrix(1:2, nrow = 2)) a= bothMarks a[,1:3] = a[,1:3] - a[,1] a[,4:6] = a[,4:6] - a[,4] a_inPath =
    # a[intersect(ensg_inPath,rownames(a)),] boxplot(a) boxplot(a_inPath) a= bothMarks a[,1:3] = a[,1:3] - a[,2]
    # a[,4:6] = a[,4:6] - a[,5] a_inPath = a[intersect(ensg_inPath,rownames(a)),] boxplot(a) boxplot(a_inPath) a=
    # bothMarks a[,1:3] = a[,1:3] - a[,3] a[,4:6] = a[,4:6] - a[,6] a_inPath =
    # a[intersect(ensg_inPath,rownames(a)),] boxplot(a) boxplot(a_inPath)
    boxplot.marks = function(index, drawLabels = F) {
        a = bothMarks
        index = index - 1
        a[, c(1, 4)] = a[, c(1, 4)] - a[, 1 + 3 * index]
        a[, c(2, 5)] = a[, c(2, 5)] - a[, 2 + 3 * index]
        a[, c(3, 6)] = a[, c(3, 6)] - a[, 3 + 3 * index]
        a_inPath = a[intersect(ensg_inPath, rownames(a)), ]
        a = a[intersect(ensg_rand, rownames(a)), ]
        boxplot(a, at = 1:6, names = rep("", 6), range = r)
        if (drawLabels) 
            axis(side = 1, labels = colnames(a), at = 1:6, las = 2)
        boxplot(a_inPath, at = 1:6, names = rep("", 6), range = r)
        if (drawLabels) 
            axis(side = 1, labels = paste(colnames(a), "subset"), at = 1:6, las = 2)
    }
    
    
    layout(matrix(c(1:4, 5, 5), nrow = 3, ncol = 2, byrow = T))
    par(mai = c(0, 0, 0, 0))
    for (i in 1:2) {
        boxplot.marks(i, drawLabels = i == 2)
    }
    
    
    
    # a= bothMarks a[,c(1,4)] = a[,c(1,4)] - a[,1] a[,c(2,5)] = a[,c(2,5)] - a[,2] a[,c(3,6)] = a[,c(3,6)] -
    # a[,3] a_inPath = a[intersect(ensg_inPath,rownames(a)),] boxplot(a) boxplot(a_inPath) a= bothMarks
    # a[,c(1,4)] = a[,c(1,4)] - a[,4] a[,c(2,5)] = a[,c(2,5)] - a[,5] a[,c(3,6)] = a[,c(3,6)] - a[,6] a_inPath =
    # a[intersect(ensg_inPath,rownames(a)),] boxplot(a) boxplot(a_inPath)
    dev.off()
}

# Function to plot color bar
color.bar <- function(lut, min, max = -min, nticks = 11, ticks = seq(min, max, len = nticks), title = "") {
    scale = (length(lut) - 1)/(max - min)
    
    # dev.new(width=1.75, height=5)
    plot(c(0, 10), c(min, max), type = "n", bty = "n", xaxt = "n", xlab = "", yaxt = "n", ylab = "", main = title)
    axis(2, las = 1)
    for (i in 1:(length(lut) - 1)) {
        y = (i - 1)/scale + min
        rect(0, y, 10, y + 1/scale, col = lut[i], border = NA)
    }
} 

gsea.plotListMembership = function(gseaLists, shared_threshold = 0.15, doVennMatrix = F) {
  load("gsea_membershipTable.save")
  passed_membership = membership[, gseaLists]  #calculate subset that passed high p value threshold
  keep = rowSums(passed_membership) > 0
  passed_membership = passed_membership[keep, ]
  
  do = ncol(passed_membership)
  
  a = matrix(1, do, do)  #construct matrices of overlaps between gene lists
  
  rownames(a) = colnames(passed_membership)
  colnames(a) = colnames(passed_membership)
  
  if (doVennMatrix) {
    overlapL = a
    overlapR = a
    
    for (i in 1:do^2) {
      rowNum = floor((i - 1)/do) + 1
      colNum = (i + rowNum - 1)%%(do + 1)
      print(paste(rowNum, colNum))
      if (rowNum == colNum) {
        next
      } else {
        vc = vennCounts(passed_membership[, c(rowNum, colNum)])
        lVal = vc[4, 3]/(vc[2, 3] + vc[4, 3])
        rVal = vc[4, 3]/(vc[3, 3] + vc[4, 3])
        overlapL[rowNum, colNum] = lVal
        overlapR[rowNum, colNum] = rVal
      }
    }
    
    
    
    layout(matrix(1:do^2, ncol = do, byrow = T))  #plot matrix of overlaps
    par(mai = rep(0, 4), xpd = T)
    
    a = overlapL + overlapR
    o = order(rowSums(a), decreasing = T)
    overlapL = overlapL[o, o]
    overlapR = overlapR[o, o]
    
    for (r in 1:nrow(overlapL)) {
      for (c in 1:ncol(overlapL)) {
        if (r == c) {
          plot(c(0, 1), c(0, 1), axes = F, type = "n")
          name = colnames(overlapL)[c]
          w = 8
          # if(w > 2){
          suff = substr(name, nchar(name) - w, nchar(name))
          mid = substr(name, w + 1, min(w + w, nchar(name)))
          mid2 = substr(name, 2 * w + 1, min(3 * w, nchar(name)))
          begin = strtrim(name, w)
          
          # }
          
          text(0.5, 0.85, begin, cex = 0.5)  #, cex = 1/w/2)
          if (nchar(name) > w) 
            text(0.5, 0.6, mid, cex = 0.5)  #, cex = 1/w/2)
          if (nchar(name) > 2 * w) 
            text(0.5, 0.35, mid2, cex = 0.5)  #, cex = 1/w/2)
          if (nchar(name) > 3 * w) 
            text(0.5, 0.1, suff, cex = 0.5)
        } else if (r > c) 
          plot(c(0, 1), c(0, 1), axes = F, type = "n") else {
            lVal = overlapL[r, c]
            lColor = (c(1, 1 - lVal, 1 - lVal))
            lColor = rgb(matrix(lColor, ncol = 3))
            rVal = overlapR[r, c]
            rColor = (c(1, 1 - rVal, 1 - rVal))
            rColor = rgb(matrix(rColor, ncol = 3))
            barplot(c(lVal, rVal), ylim = c(0, 1), axes = F, col = c(lColor, rColor), width = 0.4, xlim = c(-1, 
                                                                                                            2))
          }
      }
    }
  }
  
  
  perc_passed = rowSums(passed_membership)/ncol(passed_membership)
  o = order(perc_passed, decreasing = T)
  perc_passed = perc_passed[o]
  keep = perc_passed > shared_threshold
  perc_passed = perc_passed[keep]
  heatmap(passed_membership[names(perc_passed), ], Rowv = F, Colv = F, scale = "n", col = c("white", "red"), 
          cexCol = 0.7, margins = c(18, 3), cexRow = 17/length(perc_passed))
} 