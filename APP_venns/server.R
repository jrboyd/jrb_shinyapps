plainNumeric = "Plain Numeric"
proportional = "Proportional"

shinyServer(function(input, output, session) {
    
    getMembership = reactive({
        sel = input$dataSrc
        out = NULL
        if (sel == src_k4marks) {
            out = markData
        } else if (sel == src_perc75) {
            pass75
            pass98
            percData = cbind(pass75, pass98)
            colnames(percData) = c(paste("75th-Percentile", colnames(pass75)), paste("98th-Percentile", 
                colnames(pass98)))
            out = percData
        } else if (sel == src_uniq) {
            out = uniquely_K4_membership
        }
        return(out)
    })
    
    split_groups = function(cnames) {
        # splits input names (likely colnames) by '_' and ' ' input names must have same
        # number of elements when split first element defines group remaining elements should
        # be unique within group returns list of character arrays names of list items are
        # group names and array contents are individual names
        tmp = strsplit(cnames, "[_ ]")
        tmp = matrix(unlist(tmp), nrow = length(tmp[[1]]))
        grpNames = tmp[1, , drop = F]
        uniq = unique(grpNames[1, ])
        content = apply(tmp[2:nrow(tmp), , drop = F], 2, function(x) return(paste(x, 
            collapse = " ")))
        out = list()
        for (un in uniq) {
            keep = grpNames == un
            out[[length(out) + 1]] = content[keep]
            names(out)[length(out)] = un
        }
        return(out)
    }
    
    output$groupCheckboxes = renderUI({
        memb = getMembership()
        grps = split_groups(colnames(memb))
        tbody = ""
        for (i in 1:length(grps)) {
            grp = grps[[i]]
            nam = names(grps)[i]
            entry = (checkboxGroupInput(inputId = paste("checkboxes_", nam, sep = ""), 
                label = nam, choices = grp, selected = grp[1], inline = T))
            tbody = paste(tbody, entry)
        }
        return(column(width = 12, HTML(tbody)))
    })
    
    
    testType = function(ensgs, groups_toKeep) {
        # ensgs : character vector of ensg.x ids groups_toKeep : character vector of gene
        # types to keep returns boolean vector of rows to keep that match any gene type in
        # groups_toKeep
        groups = type2groups[ensg2type[ensgs]]
        keep = sapply(groups, function(x) {
            return(any(x == groups_toKeep))
        })
        names(keep) = NULL
        keep[is.na(keep)] = F
        return(keep)
    }
    
    getSelectedMarks = reactive({
        # returns currently selected items from checkbox groups
        memb = getMembership()
        grps = split_groups(colnames(memb))
        all_selected = character()
        for (i in 1:length(grps)) {
            nam = names(grps)[i]
            key = paste("checkboxes_", nam, sep = "")
            sel = input[[key]]
            sel = paste(nam, sel)
            print(sel)
            all_selected = append(all_selected, sel)
        }
        toPlot = all_selected
        keep = !substr(toPlot, nchar(toPlot), nchar(toPlot)) == " "
        toPlot = toPlot[keep]
    })
    
    getVennTable = reactive({
        toPlot = getSelectedMarks()
        # if(length(toPlot) > 3){ plot(c(0,1),c(0,1)) text(.5,.5,'only 3 groups allowed!')
        # return() }
        if (length(toPlot) < 1) {
            plot(c(0, 1), c(0, 1))
            text(0.5, 0.5, "at least 2 group required!")
            return()
        }
        memb = getMembership()
        if (class(memb[1, 1]) == "numeric") 
            memb = memb > input$threshold
        keep = testType(rownames(memb), input$typeFilter)
        final = memb[keep, toPlot]
        return(final)
    })
    
    plotVenn = function(vennTable, plotType, labels_on_venn) {
        final = vennTable
        if (plotType == plainNumeric) {
            if (labels_on_venn) {
                vennDiagram(final)
            } else {
                colors = RColorBrewer::brewer.pal(ncol(final), name = "Dark2")
                vennDiagram(final, names = "", circle.col = colors)
                legend(x = "topleft", legend = colnames(final), fill = colors)
            }
        } else {
            ve = venneuler(final)
            colors = hsv(ve$colors, s = 0.2)
            if (labels_on_venn) {
                plot(ve)
            } else {
                ve$labels = rep("", length(ve$labels))
                plot(ve)
                legend(x = "topleft", legend = colnames(final), fill = colors)
            }
        }
    }
    
    output$venn = renderPlot({
        final = getVennTable()
        plotVenn(final, input$plotType, input$labels_on_venn)
        # 
    })
    
    output$dl_venn = downloadHandler(filename = {
        fname = dl_name()
        fname = paste("venn_", fname, ".pdf", sep = "")
    }, content = function(file) {
        pdf(file)
        final = getVennTable()
        plotVenn(final, input$plotType, input$labels_on_venn)
        dev.off()
    })
    
    output$list_size = renderText({
        gl = geneList()
        return(paste(length(gl), "genes in result"))
    })
    
    getConditions = reactive({
        toPlot = getSelectedMarks()
        if (length(toPlot) < 1) {
            return()
        }
        conditions = toPlot
        names(conditions) = toPlot
        for (tp in toPlot) {
            key = paste("radio_", tp, sep = "")
            conditions[tp] = input[[key]]
        }
        return(conditions)
    })
    
    geneList = reactive({
        conditions = getConditions()
        dat = getVennTable()
        keep = conditions != "either"
        if (sum(keep) == 0) {
            # all are either, no filtering done
            return(rownames(dat))
        }
        dat = dat[, keep, drop = F]
        conditions = conditions[keep]
        keys = conditions == "pos"
        keep = apply(dat, 1, function(x) return(all(x == keys)))
        dat = dat[keep, , drop = F]
        return(rownames(dat))
    })
    
    output$list_conditions = renderUI(expr = {
        toPlot = getSelectedMarks()
        if (length(toPlot) == 0) {
            return()
        }
        tbody = ""
        for (i in 1:length(toPlot)) {
            tp = toPlot[i]
            entry = (radioButtons(inputId = paste("radio_", tp, sep = ""), label = tp, 
                choices = c("pos", "neg", "either"), selected = "pos", inline = T))
            tbody = paste(tbody, entry)
        }
        return(column(width = 12, HTML(tbody)))
    })
    
    dl_name = reactive({
        
        
    })
    
    output$dl_list = downloadHandler(filename = {
        fname = dl_name()
        fname = "gene_list.txt"
    }, content = function(file) {
        gl = geneList()
        out = ensg2sym[gl]
        write.table(out, file = file, row.names = F, col.names = F, quote = F, sep = "\t")
    })
    
    output$dl_table = downloadHandler(filename = {
        fname = dl_name()
        fname = "membership_table.csv"
    }, content = function(file) {
        dat = getVennTable()
        gl = geneList()
        out = dat[gl, ]
        out = cbind(ensg2sym[gl], gl, out)
        colnames(out)[1:2] = c("Gene Symbol", "ENSG id")
        write.table(out, file = file, row.names = F, col.names = T, quote = F, sep = ",")
    })
    
    
    output$list_copyable = renderUI({
        type = input$geneList_type
        gl = geneList()
        str = paste(ensg2sym[gl], collapse = ",")  #csv symbols is default
        if (type == gl_types[2]) {
            str = paste(gl, collapse = ",")  #raw ensgs
        } else if (type == gl_types[3]) {
            # cut ensgs
            str = paste(ensg2cut[gl], collapse = ",")
        } else if (type == gl_types[4]) {
            # david query
            str = "TODO"
        }
        textInput(inputId = "sigGenes_text", "Copyable Gene List:", str)
    })
}) 
