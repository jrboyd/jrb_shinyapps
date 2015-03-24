


rects = rects.get(vertLines = c(-1,1), horLines = 1, XMIN = -4, XMAX = 4, YMAX = 6)

plot_stainedGlassOverlay = function(gsea_name, pairName){
  
  query = pairName
  #   for(group_i in 1:length(grps)){
  #firstGrp = getParamsByList(group_i)
  grp_data = data[,,gsea_name,,drop = F]
  for(list_i in 1:length(gsea_name)){
    
    obs_dat = grp_data[,,list_i,2]
    exp_dat = grp_data[,,list_i,3]
    
    paramTabs = grp_data[,,list_i,]
    ratios = (obs_dat - exp_dat)
    csums = colSums(obs_dat)
    ratios = ratios / csums
    ratios = ratios[,query, drop = F]
    pvals = grp_data[,,list_i,1]
    pvals = pvals[,query, drop = F]
    
    obs_dat = obs_dat[,query, drop = F]
    exp_dat = exp_dat[,query, drop = F]
    
    MIN = min(ratios)
    MAX = max(ratios)
    
    
    #layout(matrix(1:(2*length(query)), ncol = 2*length(query)),widths = rep(c(5,1), length(query)))
    for(q in 1:length(query)){
      
      #plot(c(-4,4), c(0,6), type = 'n', xlab = 'log2 Fold Change (FC)', ylab = 'log2 minimum Fold Enrichment (FE)')
      #title(paste(gsea_name, query, sep = '\n'),sub = 'Labels are -log10 pvalues')
      rat = ratios[,q]
      pv = pvals[,q]
      obs = obs_dat[,q]
      exp = round(exp_dat[,q],digits = 1)
      cfactor = 1/max(abs(ratios))
      cr = function(x){
        x = x * cfactor
        x = x / 2 + .5
        x = colorRamp(colors = c('#2166ac', "#FFFFFF", '#b2182b'))(x)
        x = rgb(x/255)
        return(x)
      }
      colors = sapply(rat, function(x){x = x * cfactor; return(ifelse(x > 0, yes = rgb(1-x,1,1-x), no = rgb(1,1+x,1+x)))} )
      colors = cr(rat)
      colors = paste(colors, '4D', sep = '')
      for(i in 1:length(rects)){
        r = rects[[i]]
        center = rects.center(r)
        rects.draw(r, col = colors[i])
        
        #text(center[1], center[2]-.2, round(rat[i],1))
      }
      for(i in 1:length(rects)){
        r = rects[[i]]
        center = rects.center(r)
        if(pv[i] >= 2){
          expFactor = 1
          draw.circle(center[1],center[2], r = .3*expFactor, col = rgb(1,1,1))
          text(center[1], center[2], paste(obs[i],'/', exp[i]), cex = expFactor)
        }
        else
          text(center[1], center[2], paste(obs[i],'/', exp[i]))
      }
      
      #color.bar(colorRampPalette(c('#2166ac', "white", '#b2182b'))(100),-(max(abs(ratios))))
      #title(sub = 'Color Key', cex.main = .8,ylab = 'Fraction of List Shifted')
    }
    #     }
    
  }
}

plot_stainedGlass = function(gsea_name, pairName){

  query = pairName
  #   for(group_i in 1:length(grps)){
  #firstGrp = getParamsByList(group_i)
  grp_data = data[,,gsea_name,,drop = F]
  for(list_i in 1:length(gsea_name)){
    
    obs_dat = grp_data[,,list_i,2]
    exp_dat = grp_data[,,list_i,3]
    
    paramTabs = grp_data[,,list_i,]
    ratios = (obs_dat - exp_dat)
    csums = colSums(obs_dat)
    ratios = ratios / csums
    ratios = ratios[,query, drop = F]
    pvals = grp_data[,,list_i,1]
    pvals = pvals[,query, drop = F]
    
    
    MIN = min(ratios)
    MAX = max(ratios)
    
    
    layout(matrix(1:(2*length(query)), ncol = 2*length(query)),widths = rep(c(5,1), length(query)))
    for(q in 1:length(query)){
      
      plot(c(-4,4), c(0,6), type = 'n', xlab = 'log2 Fold Change (FC)', ylab = 'log2 minimum Fold Enrichment (FE)')
      title(paste(gsea_name, query, sep = '\n'),sub = 'Labels are -log10 pvalues')
      rat = ratios[,q]
      pv = pvals[,q]
      
      cfactor = 1/max(abs(ratios))
      cr = function(x){
        x = x * cfactor
        x = x / 2 + .5
        x = colorRamp(colors = c('#2166ac', "white", '#b2182b'))(x)
        x = rgb(x/255)
        return(x)
      }
      colors = sapply(rat, function(x){x = x * cfactor; return(ifelse(x > 0, yes = rgb(1-x,1,1-x), no = rgb(1,1+x,1+x)))} )
      colors = cr(rat)
      for(i in 1:length(rects)){
        r = rects[[i]]
        center = rects.center(r)
        rects.draw(r, col = colors[i])
        
        #text(center[1], center[2]-.2, round(rat[i],1))
      }
      for(i in 1:length(rects)){
        r = rects[[i]]
        center = rects.center(r)
        if(pv[i] >= 2){
          expFactor = 1
          draw.circle(center[1],center[2], r = .3*expFactor, col = rgb(1,1,1))
          text(center[1], center[2], round(pv[i],1), cex = expFactor)
        }
        else
          text(center[1], center[2], round(pv[i],1))
      }
      
      color.bar(colorRampPalette(c('#2166ac', "white", '#b2182b'))(100),-(max(abs(ratios))))
      title(sub = 'Color Key', cex.main = .8,ylab = 'Fraction of List Shifted')
    }
    #     }
    
  }
}
