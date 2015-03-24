
shinyServer(function(input, output, session) {
  
  active_columns = reactive({
    a = paste('MCF7 vs', input$mcf7_updown)
    b = paste('MDA231 vs', input$mda231_updown)
    c = paste('MDA231-siNS vs', input$mda231siNS_updown)
    selection = character()
    if(length(input$mcf7_updown) > 0)
      selection = append(selection, a)
    if(length(input$mda231_updown) > 0)
      selection = append(selection, b)
    if(length(input$mda231siNS_updown) > 0)
      selection = append(selection, c)
    return(selection)
  })
  
  active_rows = reactive({
    type = input$type
    keep = rep(T, nrow(padj))
    if(type == 'Protein Coding'){
      keep = type2groups[ensg2type[rownames(padj)]] == "coding"
    }
    else if(type == 'lncoding'){
      keep = type2groups[ensg2type[rownames(padj)]] == "lnoncoding"
    }
    return(keep)
  })
  
  test_results = reactive({
    padj_t = input$padj_threshold
    fc_t = input$fc_threshold
    maxes_t = input$maxes_threshold
    
    sel = active_columns()
    keep = active_rows()
    
    sel_padj = padj[keep,sel, drop = F]
    sel_fc = fc[keep,sel, drop = F]
    sel_maxes = maxes[keep,sel, drop = F]
    
    pass_padj = sel_padj > padj_t
    fc_dir = input$updown
    fc_test = function(){
      return(abs(sel_fc) > fc_t)
    }
    if(fc_dir == 'up'){
      fc_test = function(){
        return(sel_fc >= fc_t)
      }
    }
    else if(fc_dir == 'down'){
      fc_test = function(){
        return(sel_fc <= -fc_t)
      }
    }
    else if(fc_dir == 'either'){
      fc_test = function(){
        return(abs(sel_fc) >= fc_t)
      }
    }
    else if(fc_dir == 'no change'){
      fc_test = function(){
        return(abs(sel_fc) < fc_t)
      }
    }
    else{
      print('unrecognized fc direction!')
    }
    pass_fc = fc_test()
    
    pass_maxes = sel_maxes > maxes_t
    
    pass_all = pass_padj & pass_fc & pass_maxes
    return(pass_all)
  })
  
  get_vennCounts = reactive({
    return(vennCounts(test_results()))
  })
  
  output$venn = renderPlot({
    vc = get_vennCounts()
    vc[1,ncol(vc)] = NA
    vennDiagram(vc, cex = c(1,1,.5))
    
  })  
  
  output$volcano = renderPlot({
    ydat = rowMeans(maxes[active_rows(), active_columns(), drop = F])
    xdat = rowMeans(fc[active_rows(), active_columns(), drop = F])
    ensgs = selected_list()
    plot(xdat, ydat, col = rgb(0,0,0,.2), pch = 16, xlab = 'log2 FC', ylab = 'log2 max expression')
    points(xdat[ensgs], ydat[ensgs], col = 'red', pch = 16, cex = 1.2)
    
  })
  
  output$select_gene_list = renderUI({
    vc = isolate(get_vennCounts())
    active_columns()
    keys = vc[2:nrow(vc), 1:(ncol(vc)-1), drop = F]
    
    res = test_results()
    all_choices = character()
    for(i in 1:nrow(keys)){
      k = keys[i,] > 0
      choice = paste(colnames(res)[k], collapse = ' & ')
      all_choices = append(all_choices, choice)
    }
    return(selectInput(width = '50%',
      inputId = 'selected_gene_list', 
      label = 'Select group for export', 
      choices = all_choices,
      selected = all_choices[1]
      ))    
  })
  
  dl_name = reactive({
    fname = input$selected_gene_list
    fname = paste(
                  fname, 
                  'padj', input$padj_threshold, 
                  'fc', input$fc_threshold, input$updown,
                  'atLeast', input$maxes_threshold,
                  sep = '_'
                  )
    return(fname)
  })
  
  pretty_name = reactive({
    fname = paste(
      'Adjusted p-value <', 10^-input$padj_threshold, 
      '\nFold Change', 2^input$fc_threshold, input$updown,
      '\nDetection requirement ', 2^input$maxes_threshold,
      sep = ' '
    )
    fname = sub('either', 'up or down', fname)
    if(input$updown == 'no change')
      fname = sub('no change', 'less than', fname)
    return(fname)
  })
  
  dl_listname = reactive({
    fname = dl_name()
    fname = paste('list_',fname, '.txt', sep = '')
  })
  
  dl_tablename = reactive({
    fname = dl_name()
    fname = paste('table_',fname, '.csv', sep = '')
  })
  
  dl_vennname = reactive({
    fname = dl_name()
    fname = paste('venn_',fname, '.pdf', sep = '')
  })
  
  selected_list = reactive({
    
    sel = input$selected_gene_list
    if(is.null(sel)) return()
    tmp = strsplit(sel, split = ' & ')[[1]]
    tr = test_results()
    key = rep(F, ncol(tr))
    names(key) = colnames(tr)
    key[tmp] = T
    keep = t(apply(tr, 1, function(x)return(all(x == key))))
    return(rownames(tr)[keep])
  })
  
  content_list = function(file){
    ensg = selected_list()
    write.table(ensg2sym[ensg], file = file, row.names = F, col.names = F, quote = F)
  }
  
  content_table = function(file){
    ensg = selected_list()
    ac = active_columns()
    out = cbind(ensg, ensg2sym[ensg])
    colnames(out) = c('ensg', 'gene_symbol')
    for(i in 1:length(ac)){
      out = cbind(out, padj[ensg,ac[i], drop = F], fc[ensg,ac[i], drop = F], maxes[ensg,ac[i], drop = F])
      colnames(out)[(3+3*(i - 1)):(5+3*(i - 1))] = paste(ac[i], c('log10 padj', 'log2 fc', 'log2 max'))
    }
    write.table(out, file = file, row.names = F, col.names = T, quote = F, sep =',')
  }
  
  content_venn = function(file){
    pdf(file)
    vc = get_vennCounts()
    vc[1,ncol(vc)] = NA
    vennDiagram(vc, cex = c(.6,.5,.2))
    title(pretty_name())
    dev.off()
  }
  
  output$dl_list = downloadHandler(
    filename = dl_listname,
    content = content_list
  )
  
  output$dl_table = downloadHandler(
    filename = dl_tablename,
    content = content_table
  )
  
  output$dl_venn = downloadHandler(
    filename = dl_vennname,
    content = content_venn
  )
  
  dl_volcname = reactive({
    fname = dl_name()
    fname = paste('volcano_',fname, '.pdf', sep = '')
  })
  
  content_volc = function(file){
    pdf(file)
    ydat = rowMeans(maxes[active_rows(), active_columns(), drop = F])
    xdat = rowMeans(fc[active_rows(), active_columns(), drop = F])
    ensgs = selected_list()
    plot(xdat, ydat, col = rgb(0,0,0,.2), pch = 16, xlab = 'log2 FC', ylab = 'log2 max expression')
    points(xdat[ensgs], ydat[ensgs], col = 'red', pch = 16, cex = 1.2)
    title(pretty_name())
    dev.off()
  }
  
  output$dl_volcano = downloadHandler(
    filename = dl_volcname,
    content = content_volc
  )
  
  
  
  
})