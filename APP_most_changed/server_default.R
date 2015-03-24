source('heatmap.3-split.R')
# markData_4me3_4ac = readRDS('data/markData_4me3_4ac.rds')
# assign('markData_4me3_4ac', markData_4me3_4ac, envir = .GlobalEnv)
ensg2sym = readRDS('data/ensg2sym.rds')
assign('ensg2sym', ensg2sym, envir = .GlobalEnv)

shinyServer(function(input, output, session) {
  fetchClass=function(n, res){
    # data: original data
    # n: the class number from line plot(1 at top, last at bottom)
    # classDesc: [[1]] returned by makeFigure
    # classAssign: [[2]] returned by makeFigure
    # return list of values for all class members, values relative to column 1, and mean of each column's values relative to column 1
    classSizes = res[[1]]
    asPlotted = res[[3]]
    start=1
    if(n > 1)
      start = sum(classSizes[1:(n-1)])+1
    end=sum(classSizes[1:n])
    nth_class = asPlotted[start:end,]
    return(nth_class)
  }
  # Combine the selected variables into a new data frame
  get_markData = reactive({
    markData_4me3_4ac = readRDS('data/markData_4me3_4ac.rds')
    return(markData_4me3_4ac)
  })
  
  getDat <- reactive({
    get_markData()[, input$p2] - get_markData()[, input$p1]
  })
  
  getNumClasses <- reactive({
    input$clusters
  })
  
  getKeep = reactive({
    dat = getDat()
    threshold = getThreshold()
    dat > log2(threshold) | dat < -log2(threshold)
  })
  
  getThreshold = reactive({
    input$threshold
  })
  
  
  #   output$plot1 <- renderPlot({
  #     par(mar = c(5.1, 4.1, 0, 1))
  #     plot(selectedData(),
  #          col = clusters()$cluster,
  #          pch = 20, cex = 3)
  #     points(clusters()$centers, pch = 4, cex = 4, lwd = 4)
  #   })
  getTitle = reactive({
    paste(
      sum(getKeep()),'promoters\n',
      'meeting', getThreshold(), 'fold change\nfrom', input$p2, 'to', input$p1)
  })
  ia = 0
  getRes = function(){
    return(get('res', .GlobalEnv))
  }
  
  do_copyENSGs = observeEvent(input$copyENSGs,{
    if(is.null(input$classSelector))
      return(NULL)
    selected = as.numeric(isolate(input$classSelector))
    res = getRes()
    i = 1
    n_dat = fetchClass(selected[i], res)
    while(i < length(selected)){
      i = i +1
      n_dat = rbind(n_dat, fetchClass(selected[i], res))
    }
    output = cbind(rownames(n_dat), ensg2sym[rownames(n_dat)])
    ensgs = output[,1]
    ensgs = matrix(unlist(strsplit(ensgs, '\\.')), nrow = 2)[1,]
    writeClipboard(ensgs)
  })
  
  do_copyNames = observeEvent(input$copyNames,{
    if(is.null(input$classSelector))
      return(NULL)
    selected = as.numeric(isolate(input$classSelector))
    res = getRes()
    i = 1
    n_dat = fetchClass(selected[i], res)
    while(i < length(selected)){
      i = i +1
      n_dat = rbind(n_dat, fetchClass(selected[i], res))
    }
    output = cbind(rownames(n_dat), ensg2sym[rownames(n_dat)])
    writeClipboard(output[,2])
  })
  
  do_gotoDAVID = observeEvent(input$gotoDAVID,{
    if(is.null(input$classSelector))
      return(NULL)
    selected = as.numeric(isolate(input$classSelector))
    res = getRes()
    i = 1
    n_dat = fetchClass(selected[i], res)
    while(i < length(selected)){
      i = i +1
      n_dat = rbind(n_dat, fetchClass(selected[i], res))
    }
    output = cbind(rownames(n_dat), ensg2sym[rownames(n_dat)])
    ensgs = output[,1]
    ensgs = matrix(unlist(strsplit(ensgs, '\\.')), nrow = 2)[1,]#trim after .
    type = 'ENSEMBL_GENE_ID'
    annot = 'GOTERM_BP_FAT,KEGG_PATHWAY,PANTHER_PATHWAY,'
    tool = 'term2term'
    ids = paste(ensgs,collapse = ',')
    address = paste('http://david.abcc.ncifcrf.gov/api.jsp?',
                    'type=',type,
                    '&ids=',ids,
                    '&tool=',tool,
                    '&annot=',annot, sep = '')
    writeClipboard(address)
    browseURL(address)
  })
  
  do_writePdf = observeEvent(input$writePdf,{
    fname = input$pdfName
    sub('.pdf','',fname)
    fname = paste(fname, '.pdf', sep = '')
    pdf(fname)
    if(input$p1 == input$p2){
      plot(c(0,1), c(0,1), type = 'n')
      text(.5,.5,'selections must be different!')
      return()
    }
    nclasses = getNumClasses()
    dat = getDat()
    threshold = getThreshold()
    keep = getKeep()
    if(sum(keep) < nclasses){
      plot(c(0,1), c(0,1), type = 'n')
      text(.5,.5,'threshold too stringent!')
      return()
    }
    myclust = function(x){
      hclust(x/rowSums(x))
    }
    options(warn = -1)
    sink('NUL')
    res_tmp = heatmap.3(get_markData()[getKeep(), ], hclustfun = myclust, trace = "n", classCount = getNumClasses(), main = getTitle(), margins = c(.05,2), cexCol = 3, key = F, colsep = 3, key.xlab = 'log2 FE')
    tmp = res_tmp[[3]]
    rownames(tmp) = ensg2sym[rownames(tmp)]
    res_tmp[[3]] = tmp
    plot.HeatmapLists(res_tmp)
    sink()
    options(warn = 0)
    dev.off()
  })
  
  
  output$heatmap = renderPlot({
    if(input$p1 == input$p2){
      plot(c(0,1), c(0,1), type = 'n')
      text(.5,.5,'selections must be different!')
      return()
    }
    nclasses = getNumClasses()
    dat = getDat()
    threshold = getThreshold()
    keep = getKeep()
    if(sum(keep) < nclasses){
      plot(c(0,1), c(0,1), type = 'n')
      text(.5,.5,'threshold too stringent!')
      return()
    }
    options(warn = -1)
    
    
    sink('NUL')
    res = heatmap.3(get_markData()[getKeep(), ], trace = "n", classCount = getNumClasses(), main = getTitle(), margins = c(.05,2), cexCol = 3, key = F, colsep = 3)
    assign('res', res, envir = .GlobalEnv)
    sink()
    options(warn = 0)
  })
  
  
  
  output$subsetList = renderTable({
    getDat()
    input$p1#attach line/mark selectors
    input$p2
    input$threshold
    input$clusters
    if(is.null(input$classSelector))
      return(NULL)
    selected = as.numeric(input$classSelector)
    res = getRes()
    i = 1
    n_dat = fetchClass(selected[i], res)
    while(i < length(selected)){
      i = i +1
      n_dat = rbind(n_dat, fetchClass(selected[i], res))
    }
    output = cbind(rownames(n_dat), ensg2sym[rownames(n_dat)])
    rownames(output) = NULL
    return(output)
    
  })
  
  #   output$subsetList = renderPlot({
  #     justList = T
  #     selected = as.numeric(input$classSelector)
  #     res = getRes()
  #     i = 1
  #     print(selected[i])
  #     n_dat = fetchClass(selected[i], res)
  #     while(i < length(selected)){
  #       i = i +1
  #       print(selected[i])
  #       n_dat = rbind(n_dat, fetchClass(selected[i], res))
  #     }
  #     plotTitle = paste(nrow(n_dat), 'selected promoters')
  #     if(justList){
  #       res_sub = invisible(heatmap.3(n_dat, trace = "n", classCount = 4, main = plotTitle, margins = c(.05,2), cexCol = 3, key = F, colsep = 3))
  #     plot.HeatmapLists(res_sub)
  #     }else{
  #       heatmap.3(n_dat, trace = "n", classCount = 4, main = plotTitle, margins = c(.05,2), cexCol = 3, key = F, colsep = 3)
  #     }
  #   })
  
  output$classSelector = renderUI({
    getDat()
    input$p1#attach line/mark selectors
    input$p2
    input$threshold
    input$clusters
    getNumClasses()
    getThreshold()
    res = getRes()
    plottedColors = res[[4]]
    selectInput(inputId = 'classSelector', label = 'Select clusters for table and copying', choices = 1:length(res[[1]]), multiple = T, selected = 1)
  })
  
  output$numVehiclesTable <- renderUI({
    getDat()
    input$p1#attach line/mark selectors
    input$p2
    input$threshold
    input$clusters
    #print(getNumClasses())
    res = getRes()
    #print(is.null(res))
    if(is.null(res)){
      return(NULL)
    }
    plottedClasses = res[[1]]
    plottedColors = res[[4]]
    asPlotted = res[[3]]
    
    # Create a Bootstrap-styled table
    tbody = ""
    for(i in 1:length(plottedClasses)){
      entry = tags$tr(
        tags$td(i),
        tags$td(span(style = sprintf(
          "width:1.1em; height:1.1em; background-color:%s; display:inline-block;",
          plottedColors[i]
        ))),
        tags$td(plottedClasses[i]))
      tbody = paste(tbody, entry, sep = '\n')
    }
    
    tags$table(class = "table",
               tags$thead(tags$tr(
                 tags$th("Cluster ID"),
                 tags$th("Color"),
                 tags$th("Size")
               )),
               tags$tbody(
                 HTML(tbody)
               )
               #,
               #                tags$tr(class = "active",
               #                        tags$td(),
               #                        tags$td("Total"),
               #                        tags$td(length(res[[1]]))
               #                )
    )
    
  })
})