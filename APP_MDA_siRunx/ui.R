source('setup.R')

shinyUI(fluidPage(
  headerPanel('MDA Runx1 Knockdown'),
  fluidRow(
    column(width = 6, radioButtons(inputId = 'type', label = 'Select Gene Type', choices = c('Protein Coding', 'lncoding', 'Both'), selected = 'Protein Coding'))
  ),
  fluidRow(
    column(width = 3,radioButtons(inputId = 'updown', label = 'Fold change direction', choices = c('up', 'down', 'either', 'no change'), selected = 'either')),
    column(width = 3,checkboxGroupInput('mcf7_updown', label = 'MCF7 to ...', choices = comp_to[1:3])),
    column(width = 3,checkboxGroupInput('mda231_updown', label = 'MDA231 to ...', choices = comp_to[4:5], selected = 'MDA231-siRunx1')),
    column(width = 3,checkboxGroupInput('mda231siNS_updown', label = 'MDA231-siNS to ...', choices = comp_to[6], selected = 'MDA231-siRunx1'))
  ),
  fluidRow(
    column(width = 4, sliderInput('padj_threshold', label = '-log10 p-value threshold', min = 0, max = 9, value = 2)),
    column(width = 4, sliderInput('fc_threshold', label = 'log2 fold-change threshold', min = 0, max = 10, value = 1)),
    column(width = 4, sliderInput('maxes_threshold', label = 'log2 max threshold', min = 0, max = 16, value = 2))
  ),
  fluidRow(
    column(plotOutput('venn'), width = 6),
    column(plotOutput('volcano'), width = 6)
  ),
  uiOutput('select_gene_list'),
  fluidRow(
    downloadButton('dl_list', 'Download Selected List'),
    downloadButton('dl_table', 'Download Selected Table'),
    downloadButton('dl_venn', 'Download Venn Diagram'),
    downloadButton('dl_volcano', 'Download Volcano Plot')
  )
)

)