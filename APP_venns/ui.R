source("setup.R")
source("globals.R")

shinyUI(fluidPage(headerPanel("Venn Diagrams"), radioButtons("dataSrc", "Select Data Source", 
    choices = c(src_k4marks, src_perc75, src_uniq), selected = src_perc75, inline = T), 
    sidebarLayout(sidebarPanel(width = 6, tags$h3("Plot Controls"), radioButtons("plotType", 
        "Select Plot Type", choices = c(plainNumeric, proportional), selected = plainNumeric, 
        inline = T), uiOutput("groupCheckboxes"), checkboxGroupInput("typeFilter", label = "Filter by gene type", 
        choices = unique(type2groups), selected = "coding", inline = T), sliderInput("threshold", 
        label = "Adjust Threshold", min = 0.5, max = 5, value = 1, step = 0.5)), mainPanel(width = 6, 
        plotOutput("venn"), fluidRow(downloadButton("dl_venn", "Download Plotted Venn"), 
            checkboxInput("labels_on_venn", label = "Label Venn", value = F)))), fluidRow(h1("List Export"), 
        column(width = 4, h3("Filter"), uiOutput("list_conditions")), column(width = 5, 
            h3("Copy"), textOutput("list_size"), radioButtons("geneList_type", label = "Copyable gene list format:", 
                choices = gl_types, inline = T), uiOutput("list_copyable")), column(width = 3, 
            h3("Download"), downloadButton("dl_list", "Download Selected List"), downloadButton("dl_table", 
                "Download Selected Table"))))) 
