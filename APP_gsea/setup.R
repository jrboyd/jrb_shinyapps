load('data//mycounts_data.save')
load('data//dictionaries.save')
cut2ensg = names(ensg2cut)
names(cut2ensg) = ensg2cut
getStandardTable = function() return(markData_4me3_4ac)
getIntersect = function(table, ensgs){
  keep = intersect(rownames(table), ensgs)
  return(table[keep,])
}
load('data//data binomTable Feb 27.save')


library(shiny)
library(tm)
library(wordcloud)
library(gplots)
library(plotrix)
library(RColorBrewer)
tmp = sapply(dir('scripts', full.names = T), source)
pairNames = dimnames(data)[[2]]
name2col = 1:ncol(markData_4me3_4ac)

names(name2col) = colnames(markData_4me3_4ac)
tmp = sapply(strsplit(pairNames, ' to '), function(x)return(x[1]))
tmp = sub('From ', '', tmp)
pair2starts = tmp
names(pair2starts) = pairNames

pair2ends = sapply(strsplit(pairNames, ' to '), function(x)return(x[2]))
names(pair2ends) = pairNames

load('data//gsea_membershipTable_ensgs.save')

gl_types = c('Symbols', 'ENSG.#', 'ENSG', 'DAVID')
