source('setup.R')


shinyUI(fluidPage(
  headerPanel('GSEA Plots'),
  fluidRow(
    
    
    #   sidebarLayout(
    #     sidebarPanel(),
    #     mainPanel(plotOutput('comparisonCloud'), plotOutput('commonCloud'))
    #   ),
           column(width = 7,
                  selectInput('pair_selector', label = 'Select Mark Pair', choices = pairNames,width = '100%')
           ),
           column(width = 5,sliderInput('pval_threshold', label = 'Pvalue Threshold', min = 0, max = 10, value = 6))
           
    ),
  fluidRow(
    column(width = 8, plotOutput('sigLists')),
    column(width = 4, 
           downloadButton('dl_sigTable_hmap', 'Download List Heatmap'),
           downloadButton('dl_sigTable_csv', 'Download List Table'))
  ),
  fluidRow(
    column(width = 12, 
           uiOutput('list_selector'),
           radioButtons('plot_type',label = 'Select Plot Type', choices = c('Volcano', 'Hotspot', 'Both'),inline = T))
    ),
  
  fluidRow(
    column(width = 8, plotOutput('plots')), 
    column(width = 4, 
           downloadButton('dl_volcano', 'Download Volcano Plot'),
           checkboxInput('volcano_labels', label = 'Label Selected', value = F))
  ),
  fluidRow(
    uiOutput('regions')
  ),
  radioButtons('geneList_type',label = 'Copyable gene list format:', choices = gl_types, inline = T),
  uiOutput('geneList_action'),
  fluidRow(
    downloadButton('dl_list', 'Download Selected List'),
    downloadButton('dl_table', 'Download Selected Table')
    
  )
)

)