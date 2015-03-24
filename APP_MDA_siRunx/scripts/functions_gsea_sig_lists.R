if(F){#debugging parameters
  pairIndexes = 1:2;
  membership = membership;
  pdfName = NULL
  threshold = 5 
  doPercentagePlot = F
  doPvalPlot = T
}

plot_listMembership = function(gseaLists, membership, shared_threshold = 0.01, doVennMatrix = F, listDend = F, ...) {
  #plots list membership as a binary heatmap
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
  passed_membership = t(passed_membership[names(perc_passed),])
  o = order(colSums(passed_membership), decreasing = T)
  passed_membership = passed_membership[,o]
  if(is.matrix(passed_membership) & min(dim(passed_membership) >= 2)){
    heatmap(passed_membership, Rowv = listDend , Colv = NA, scale = "n", col = c("white", "black"), 
            cexCol = min(28/ncol(passed_membership),1), margins = c(8, 16), cexRow = min(28/nrow(passed_membership),1),  ...)
  }
  else{
    plot(c(0,1),c(0,1))
    text(.5,.5, 'Not enough overlap between lists')
  }
} 

plot_lists = function(listNames, pairIndex, membership, pdfName = NULL, doPercentagePlot = F, doPvalPlot = T, sharedCutoffs = c(.25,.1)){
  if (!is.null(pdfName)) 
    pdf(pdfName, width = 12, height = 10)
  print('a')
  lgp_dat = t(data[, pairIndex, listNames, 1])
  obs_dat = t(data[, pairIndex, listNames, 2])
  exp_dat = t(data[, pairIndex, listNames, 3])
  print('b')
  if(doPercentagePlot){
    toPlot = (obs_dat - exp_dat)/rowSums(obs_dat)
    cr = colorRamp(c("blue", "white", "red"))
    MAX = max(abs(toPlot))
    colors = cr((0:100/100))
    colors = rgb(colors/255)
    
    a = heatmap.2(toPlot, key.xlab = "Fraction Shifted", trace = "n", Colv = F, dendrogram = "r", col = colors, 
                  main = paste("            ", pairNames[pairIndex]), margins = c(20, 26), cexRow = min(22/nrow(toPlot),                   1), cexCol = 12/ncol(toPlot))
  }   
  if(doPvalPlot){
    toPlot = ifelse((obs_dat - exp_dat) > 0, lgp_dat, -lgp_dat)
    #     toPlot = ifelse(toPlot > threshold + 2, threshold + 2, toPlot)
    #     toPlot = ifelse(toPlot < -threshold - 2, -threshold - 2, toPlot)
    cr = colorRamp(c("blue", "white", "red"))
    MAX = max(abs(toPlot))
    colors = cr((0:100/100))
    colors = rgb(colors/255)
    plotTitle = pairNames[pairIndex]
    plotTitle = sub('to', '\nto', plotTitle)
    xlab = colnames(toPlot)
    start = strsplit(pairNames[pairIndex], ' to ')[[1]][1]
    start = sub('From ', '', start)
    end = strsplit(pairNames[pairIndex], ' to ')[[1]][2]
    xlab = sub('start', start, xlab)
    xlab = sub('end', end, xlab)
    colnames(toPlot) = xlab
    a = heatmap.2(toPlot, key.xlab = "Directional -log10 pval", trace = "n", Colv = F,labCol = F, key.title = "",
                  dendrogram = "r", col = colors, main = plotTitle , margins = c(2, 24), 
                  cexRow = min(22/nrow(toPlot), 1), cexCol = min(12/ncol(toPlot), 1))
  }
  #for(sc in sharedCutoffs)
  #  plot_listMembership(gseaLists = rownames(toPlot), membership = membership, listDend = a$rowDendrogram, shared_threshold = sc, main = paste('At least', sc, 'shared'))
  if (!is.null(pdfName)) 
    dev.off()
}
get_sigLists = function(pair, membership, threshold = 5) {
  # plot heatmap tables of most sig shifted lists for each mark pairing.  the lists are rows, regions as columns. 
  # color indicates direction (enriched = red, depleted = blue) with -log p labels
  #pairIndexes are the 2nd dimension index identifying histone mark pair being tested
  #membership is a binary table identifying what genes are in which list
  #if pdfName is set, output goes to a new pdf instead of current device
  #threshold is the -log10 pvalue required for a list to be kept
  #doPercentagePlot will use percentage of all genes in list that shift as heatmap values
  #doPvalPlot will use -log10 pvalues as heatmap value
  #recommend only using one of these settings
  lgp_dat = t(data[, pair, , 1])
  obs_dat = t(data[, pair, , 2])
  exp_dat = t(data[, pair, , 3])
  keep = apply(lgp_dat, 1, max) > threshold
  if(sum(keep) == 0){
    layout(1)
    plot(c(0,1),c(0,1))
    text(.5,.5,(paste('no lists met pvalue cutoff for pair:\n', dimnames(data)$pair[pair])))
    next
  }
  list_names_toPlot = rownames(lgp_dat)[keep]
  return(list_names_toPlot)
}