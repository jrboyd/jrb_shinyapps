#ui.R
library(shiny)
markData_4me3_4ac = readRDS('data/markData_4me3_4ac.rds')
ensg2sym = readRDS('data/ensg2sym.rds')

shinyUI(fluidPage(
  headerPanel('Most Changed Heatmaps'),
  fluidRow(
    sidebarPanel(
    
    
      selectInput('p1', 'Numerator Mark', colnames(markData_4me3_4ac), selected = colnames(markData_4me3_4ac)[2]),
      selectInput('p2', 'Denominator Mark', colnames(markData_4me3_4ac),
                  selected=names(iris)[[2]]),
      sliderInput('clusters', 'Cluster count', 4,
                  min = 2, max = 8, step = 1),
      sliderInput('threshold', 'Set Threshold Fold Change', 4,
                  min = 1.5, max = 10, step = .5),
    
    uiOutput("numVehiclesTable"),
    uiOutput("classSelector"),
  #  actionButton("copyENSGs", 'Copy ENSGs'),
  #  actionButton("copyNames", 'Copy Names'),
  #  actionButton("gotoDAVID", 'Go to DAVID Query'),
  #  textInput("pdfName", 'PDF name', 'heatmap.pdf'),
  #  actionButton("writePdf", 'Write PDF')),
    uiOutput('davidURL'),
    downloadButton('downloadData', 'Download')),
    mainPanel(
      column(width = 12,
      plotOutput('heatmap'),
      #plotOutput('subsetList'),
      tableOutput('subsetList')
    ))
  )
))