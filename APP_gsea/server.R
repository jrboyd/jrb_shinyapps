#loads binom test results for 6 region analysis
#identify most significantly changed lists
#group lists by biology
#for each group: for each mark pair: plot overall heatmap and stained glass plot for individual list
#use dimnames(data) to check matrix contents
#membershipTable must be loaded
#T = list contains gene, F = list does not contain gene

#plot_stainedGlass(sigLists, pairIndexes = pi)

rect_key = 1:6


shinyServer(function(input, output, session) {
  
  sigLists = reactive(get_sigLists(input$pair_selector, membership, threshold = input$pval_threshold))
  sigRegions = reactive({
    p_name = input$pair_selector
    threshold = input$pval_threshold
    l_name = input$list_selector
    if(is.null(l_name) || is.null(p_name)){
      return(NULL)
    }
    dat = data[,p_name,l_name,, drop = F]
    dat = dat[1:6,,,1:3]
    key = 1:nrow(dat)
    keep = dat[,1] > threshold #must pass pvalue threshold
    dat = dat[keep,, drop = F]
    key = key[keep]
    keep = dat[,2] > dat[,3] #must be enriched, not depleted
    dat = dat[keep,, drop = F]
    key = key[keep]
    if(nrow(dat) < 1){
      return(character())
    }
    start = pair2starts[p_name]
    end = pair2ends[p_name]
    rownames(dat) = gsub('start', start,  rownames(dat))
    rownames(dat) = gsub('end', end,  rownames(dat))
    rect_key = key
    
    return(rownames(dat))
    })
  
  sigRects = reactive({
    p_name = input$pair_selector
    threshold = input$pval_threshold
    l_name = input$list_selector
    if(is.null(l_name) || is.null(p_name)){
      return(NULL)
    }

    dat = data[,p_name,l_name,, drop = F]
    dat = dat[1:6,,,1:3]
    key = 1:nrow(dat)
    keep = dat[,1] > threshold #must pass pvalue threshold
    dat = dat[keep,, drop = F]
    key = key[keep]
    keep = dat[,2] > dat[,3] #must be enriched, not depleted
    dat = dat[keep,, drop = F]
    key = key[keep]
    if(nrow(dat) < 1){
      return(NA)
    }
    start = pair2starts[p_name]
    end = pair2ends[p_name]
    rownames(dat) = gsub('start', start,  rownames(dat))
    rownames(dat) = gsub('end', end,  rownames(dat))
    names(key) = rownames(dat)
    return(key)
  })
  
  sigGenes = reactive({
    p_name = input$pair_selector
    l_name = input$list_selector
    r_name = input$region_selector
    if(is.null(l_name) || is.null(p_name)){
      return(NULL)
    }
    rect = rects.get()[[sigRects()[r_name]]]
    
    start = pair2starts[p_name]
    end = pair2ends[p_name]
    keep = rects.countIn(x = get_x(markData_4me3_4ac, end, start),
                  y = get_y(markData_4me3_4ac, end, start),
                  rect)
    inRegion = rownames(markData_4me3_4ac)[keep]
    keep = membership[,l_name]
    inList = rownames(membership)[keep]
    final = intersect(inRegion, inList)
    
    return(final)
    
  })
  
  output$list_selector = renderUI(selectInput('list_selector', label = 'Plotted List', choices = sigLists()[order(sigLists())], width = '100%'))
  
  volcanoPlot = function(){
    plotType = input$plot_type
    
    if(plotType == 'Hotspot'){
      plot_stainedGlass(input$list_selector, input$pair_selector)
    }
    else if(plotType == 'Volcano'){
      keep = membership[,input$list_selector]
      ensg = rownames(membership)[keep]
      start = name2col[pair2starts[input$pair_selector]]
      end = name2col[pair2ends[input$pair_selector]]
      gsea.plotVolcano(ensg_inPath = ensg, startCol = start, endCol = end)
      if(input$volcano_labels){
        ensg = intersect(rownames(markData_4me3_4ac), ensg)
        xs = get_x(markData_4me3_4ac[ensg,], end, start)
        ys = get_y(markData_4me3_4ac[ensg,], end, start)
        for(i in 1:length(xs)){
          x = xs[i]
          y = ys[i]
          j = .18
          rect(x - j*2.8, y - j, x + j*2.8, y + j, col = 'white', border = NA)
        }
        text(xs, ys, ensg2sym[ensg], col = 'dark red')
      }
    }
    else{#plotBoth
      keep = membership[,input$list_selector]
      ensg = rownames(membership)[keep]
      start = name2col[pair2starts[input$pair_selector]]
      end = name2col[pair2ends[input$pair_selector]]
      gsea.plotVolcano(ensg_inPath = ensg, startCol = start, endCol = end)
      plot_stainedGlassOverlay(input$list_selector, input$pair_selector)
    }
  }
  
  conditionalPlot = reactive({
    plotType = input$plot_type
    
    if(plotType == 'Hotspot'){
      plot_stainedGlass(input$list_selector, input$pair_selector)
    }
    else if(plotType == 'Volcano'){
      keep = membership[,input$list_selector]
      ensg = rownames(membership)[keep]
      start = name2col[pair2starts[input$pair_selector]]
      end = name2col[pair2ends[input$pair_selector]]
      gsea.plotVolcano(ensg_inPath = ensg, startCol = start, endCol = end)
      if(input$volcano_labels){
        ensg = intersect(rownames(markData_4me3_4ac), ensg)
        xs = get_x(markData_4me3_4ac[ensg,], end, start)
        ys = get_y(markData_4me3_4ac[ensg,], end, start)
        for(i in 1:length(xs)){
          x = xs[i]
          y = ys[i]
          j = .18
          rect(x - j*2.8, y - j, x + j*2.8, y + j, col = 'white', border = NA)
        }
        text(xs, ys, ensg2sym[ensg], col = 'dark red')
      }
    }
    else{#plotBoth
      keep = membership[,input$list_selector]
      ensg = rownames(membership)[keep]
      start = name2col[pair2starts[input$pair_selector]]
      end = name2col[pair2ends[input$pair_selector]]
      gsea.plotVolcano(ensg_inPath = ensg, startCol = start, endCol = end)
      plot_stainedGlassOverlay(input$list_selector, input$pair_selector)
    }
  })
  
  output$plots = renderPlot(conditionalPlot())
  
  output$sigLists = renderPlot(plot_lists(sigLists(), input$pair_selector, membership))
  
  output$list_selector = renderUI({
    
    selectInput('list_selector', label = 'Plotted List', choices = sigLists()[order(sigLists())], width = '100%')
    
  })
  
  output$regions = renderUI({
    
    selectInput('region_selector', label = 'Enriched Region', choices = sigRegions(), width = '100%')
  })
  
  output$geneList_action = renderUI({
    type = input$geneList_type 
    gl = sigGenes()
    str = paste(ensg2sym[gl], collapse = ',')#csv symbols is default
    if(type == gl_types[2]){
      str = paste(gl, collapse = ',')#raw ensgs
    }
    else if(type == gl_types[3]){#cut ensgs
     str =  paste(ensg2cut[gl], collapse = ',')
    }
    else if(type == gl_types[4]){#david query
      str =  'TODO'
    }
    textInput(inputId = 'sigGenes_text', 'Copyable Gene List:', str)
  })
  
  pretty_name = reactive({
    fname = paste(
      'Adjusted p-value <', 10^(0-input$padj_threshold), 
      '\nFold Change', 2^input$fc_threshold, input$updown,
      '\nDetection requirement ', 2^(input$maxes_threshold),
      sep = ' '
    )
    fname = sub('either', 'up or down', fname)
    if(!is.null(input$updown) && input$updown == 'no change')
      fname = sub('no change', 'less than', fname)
    fname = 'tmp'
    return(fname)
  })
  
  dl_name = reactive({
    fname = input$pair_selector
    fname = paste(
      fname, 
      'pval',
      input$pval_threshold,
      sep = '_'
    )
    return(fname)
  })
  
  dl_listname = reactive({
    fname = dl_name()
    fname = paste(
      fname, 
      input$list_selector,
      sep = '_')
    return(fname)
  })
  
  dl_regname = reactive({
    fname = dl_listname()
    fname = paste(
      fname, 
      input$region_selector,
      sep = '_')
    return(fname)
  })
  
  output$dl_list = downloadHandler(
    filename = reactive({
      fname = dl_regname()
      print(fname)
      fname = paste('list_',fname, '.txt', sep = '')
    }),
    content = function(file){
      ensg = sigGenes()
      write.table(ensg2sym[ensg], file = file, row.names = F, col.names = F, quote = F)
    }
  )
  
  output$dl_table = downloadHandler(
    filename = reactive({
      fname = dl_regname()
      fname = paste('table_',fname, '.csv', sep = '')
    }),
    content = function(file){
      ensg = sigGenes()
      out = cbind(ensg, ensg2sym[ensg])
      colnames(out) = c('ensg', 'gene_symbol')
      out = cbind(out, markData_4me3_4ac[ensg,])
      write.table(out, file = file, row.names = F, col.names = T, quote = F, sep =',')
    }
  )
  
  output$dl_volcano = downloadHandler(
    filename = reactive({
      fname = dl_listname()
      fname = paste('volcano_',fname, '.pdf', sep = '')
    }),
    content = function(file){
      pdf(file, width = 16, height = 16)
      volcanoPlot()
      title(paste('Shift of', input$list_selector, input$pair_selector, sep ='\n'))
      dev.off()
    }
  )
  
  output$dl_sigTable_hmap = downloadHandler(
    filename = reactive({
      fname = dl_name()
      fname = paste('list_heatmap_',fname, '.pdf', sep = '')
    }),
    content = function(file){
      pdf(file, width = 12, height = 8)
      plot_lists(sigLists(), input$pair_selector, membership)
      title(paste('Signficantly shifted GSEA lists\n', dl_name()))
      dev.off()
    }
  )
  output$dl_sigTable_csv = downloadHandler(
    filename = reactive({
      fname = dl_name()
      fname = paste('list_table_',fname, '.csv', sep = '')
    }),
    content = function(file){
#       pdf(file, width = 12, height = 8)
#       plot_lists(sigLists(), input$pair_selector, membership)
#       dev.off()
      out = t(data[,input$pair_selector, sigLists(),1])
      write.table(out, file = file, row.names = T, col.names = T, quote = F, sep =',')
    }
  )
  
})