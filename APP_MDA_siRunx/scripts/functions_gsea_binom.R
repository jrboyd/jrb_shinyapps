getRegionParams = function(grpNum, pairNum, gseaListNum) {
    # output is list of paramter tables by region in grp, pair, and list
    pvalTab = all_groups[[grpNum]][, pairNum, gseaListNum, 1]
    obsTab = all_groups[[grpNum]][, pairNum, gseaListNum, 2]
    expTab = all_groups[[grpNum]][, pairNum, gseaListNum, 3]
    return(list(lgpval = pvalTab, observed = obsTab, expected = expTab))
}

getListParams = function(grpNum, pairNum, regionNum) {
    # output is list of paramter tables by gene list in grp, pair, and region
    pvalTab = all_groups[[grpNum]][regionNum, pairNum, , 1]
    obsTab = all_groups[[grpNum]][regionNum, pairNum, , 2]
    expTab = all_groups[[grpNum]][regionNum, pairNum, , 3]
    return(list(lgpval = pvalTab, observed = obsTab, expected = expTab))
}

getParamsByList = function(grpNum) {
    # output is list of gsea lists in input group number each entry is 3 parameter tables where rows are regions,
    # columns are pairs
    nLists = dim(all_groups[[grpNum]])[3]
    output = list()
    for (i in 1:nLists) {
        pvalTab = all_groups[[grpNum]][, , i, 1]
        obsTab = all_groups[[grpNum]][, , i, 2]
        expTab = all_groups[[grpNum]][, , i, 3]
        output[[length(output) + 1]] = list(lgpval = pvalTab, observed = obsTab, expected = expTab)
    }
    names(output) = all_listNames[[grpNum]]
    return(output)
}

gsea.binomTest = function(bg_dat, fg_dat, rectangles, sName = "start", eName = "end", doPlots = F) {
    all_rects = rectangles
    
    magnitudeMetric = function(x) {
        # out = apply(x, 1, max) #row max of x
        out = apply(x, 1, min)  #row min of x
        # out = log2(rowMeans(2^bothMarks_inPath[,c(1,2)])) #row mean of x
        return(out)
    }
    
    get_x = function(dat, a = 2, b = 1) {
        return(dat[, a] - dat[, b])
    }
    
    get_y = function(dat, a = 2, b = 1) {
        return(magnitudeMetric(dat[, c(a, b)]))
    }
    if (doPlots) {
        plot(get_x(bg_dat), get_y(bg_dat), pch = 16, col = rgb(0, 0, 0, 0.1))
        points(get_x(fg_dat), get_y(fg_dat), pch = 16, col = rgb(1, 0, 0, 0.7))
    }
    
    
    res_bg = list()
    res_fg = list()
    i = 0
    for (r in all_rects) {
        i = i + 1
        if (doPlots) 
            rects.draw(r)  #, col = colors[i])
        res_bg[[length(res_bg) + 1]] = rects.countIn(get_x(bg_dat), get_y(bg_dat), r)
        res_fg[[length(res_fg) + 1]] = rects.countIn(get_x(fg_dat), get_y(fg_dat), r)
    }
    res = list()
    for (i in 1:length(res_bg)) {
        bg = res_bg[[i]]
        fg = res_fg[[i]]
        exp_p = sum(bg)/length(bg)
        obs_n = length(fg)
        obs_suc = sum(fg)
        a = binom.test(x = obs_suc, n = obs_n, p = exp_p)
        res[[length(res) + 1]] = c(a$p.value, a$statistic, a$null.value * a$parameter)
    }
    names(res) = names(all_rects)
    return(res)
}

gsea.trackEnrichment = function(bg_dat, fg_dat, rectangles, sName = "start", eName = "end", doPlots = F) {
    all_rects = rectangles
    
    
    xRegNames = c("Down---", "Down--", "Down-", "No Change", "Up+", "Up++", "Up+++")
    yRegNames = c("both detected", "only")
    names(all_rects)[1:7] = paste(eName, xRegNames, "from", sName, ":", yRegNames[1])
    names(all_rects)[8:10] = paste(eName, xRegNames[1:3], "from", sName, ":", yRegNames[2], sName, "detected")
    names(all_rects)[11] = paste("Neither Detected")
    names(all_rects)[12:14] = paste(eName, xRegNames[5:7], "from", sName, ":", yRegNames[2], eName, "detected")
    
    # #divide area into rectangles rect_noSignal = list(xmin = -minFC, xmax = minFC, ymin = MIN, ymax = minFE)
    # rect_noChange = list(xmin = -minFC, xmax = minFC, ymin = minFE, ymax = MAX) rect_up1 = list(xmin = minFC,
    # xmax = 2*minFC, ymin = minFE, ymax = MAX) rect_up2 = list(xmin = 2*minFC, xmax = 3*minFC, ymin = minFE,
    # ymax = MAX) rect_up3 = list(xmin = 3*minFC, xmax = MAX, ymin = minFE, ymax = MAX) rect_down1 = list(xmin =
    # -2*minFC, xmax = -minFC, ymin = minFE, ymax = MAX) rect_down2 = list(xmin = -3*minFC, xmax = -2*minFC, ymin
    # = minFE, ymax = MAX) rect_down3 = list(xmin = MIN, xmax = -3*minFC, ymin = minFE, ymax = MAX)
    # rect_upSingleton = list(xmin = minFC, xmax = MAX, ymin = MIN, ymax = minFE) rect_downSingleton = list(xmin
    # = MIN, xmax = -minFC, ymin = MIN, ymax = minFE) #setup iterable list of rectangles all_rects = list(
    # rect_noSignal, rect_noChange, rect_up1, rect_up2, rect_up3, rect_down1, rect_down2, rect_down3,
    # rect_upSingleton, rect_downSingleton) names(all_rects) = c('No Signal', 'No Change', 'Up1', 'Up2', 'Up3',
    # 'Down1', 'Down2', 'Down3', 'Only in final', 'Only in start')
    
    
    
    
    # library(RColorBrewer) colors = RColorBrewer::brewer.pal(n = length(all_rects), name = 'Set3') colors =
    # paste(colors, '80', sep = '')
    
    
    
    magnitudeMetric = function(x) {
        # out = apply(x, 1, max) #row max of x
        out = apply(x, 1, min)  #row min of x
        # out = log2(rowMeans(2^bothMarks_inPath[,c(1,2)])) #row mean of x
        return(out)
    }
    
    get_x = function(dat, a = 2, b = 1) {
        return(dat[, a] - dat[, b])
    }
    
    get_y = function(dat, a = 2, b = 1) {
        return(magnitudeMetric(dat[, c(a, b)]))
    }
    if (doPlots) {
        plot(get_x(bg_dat), get_y(bg_dat), pch = 16, col = rgb(0, 0, 0, 0.1))
        points(get_x(fg_dat), get_y(fg_dat), pch = 16, col = rgb(1, 0, 0, 0.7))
    }
    
    
    res_bg = list()
    res_fg = list()
    i = 0
    for (r in all_rects) {
        i = i + 1
        if (doPlots) 
            rects.draw(r)  #, col = colors[i])
        res_bg[[length(res_bg) + 1]] = rects.countIn(get_x(bg_dat), get_y(bg_dat), r)
        res_fg[[length(res_fg) + 1]] = rects.countIn(get_x(fg_dat), get_y(fg_dat), r)
    }
    res = list()
    for (i in 1:length(res_bg)) {
        bg = res_bg[[i]]
        fg = res_fg[[i]]
        exp_p = sum(bg)/length(bg)
        obs_n = length(fg)
        obs_suc = sum(fg)
        obs_p = obs_suc/obs_n
        if (obs_p > obs_n) 
            a = binom.test(x = obs_suc, n = obs_n, p = exp_p)
        res[[length(res) + 1]] = c(a$p.value, a$statistic, a$null.value * a$parameter)
    }
    names(res) = names(all_rects)
    return(res)
} 
