load('data//dictionaries.save')
ensg_list = list()
foldchange_list = list()
padj_list = list()
max_list = list()

grow = function(l, x){#append x to list l
  l[[length(l) + 1]] = x
  return(l)
}

for(f in dir('data//DESeq', full.names = T)){
  tmp = read.table(f, sep = '\t', stringsAsFactors = F, comment.char = "", header = T)  
  ensgs = tmp[,1]
  padj = tmp[,10]
  bmean = tmp[,4]
  a = tmp[,5]
  b = tmp[,6]
  bmax = apply(cbind(a, b), 1, max)
  fc = log2(b + 1) - log2(a + 1)
  
  ensg_list = grow(ensg_list, ensgs)
  foldchange_list = grow(foldchange_list, fc)
  padj_list = grow(padj_list, padj)
  max_list = grow(max_list, bmax)
}

tmp = 'ENSG00000168490.11'
for(e in ensg_list){
  tmp = union(tmp, e)
}
all_ensg = tmp

cnames = dir('data//DESeq')[1:6*2]
tmp = strsplit(cnames, '_')
tmp = matrix(unlist(tmp), nrow = 6)
cnames = paste(tmp[2,], tmp[3,], tmp[4,])
comp_from = tmp[2,]
comp_to = tmp[4,]

assemble_matrix = function(data_list, def = 0){
  all_data = matrix(0, nrow = length(all_ensg), ncol = 0)
  rownames(all_data) = all_ensg
  for(i in 1:(length(ensg_list)/2)){
    j = (i*2)-1
    e_p = ensg_list[[j]]
    fc_p = data_list[[j]]
    e_ln = ensg_list[[j+1]]
    fc_ln = data_list[[j+1]]
    all_data = cbind(all_data, rep(def, nrow(all_data)))
    all_data[e_p,i] = fc_p
    all_data[e_ln,i] = fc_ln
  }
  colnames(all_data) = cnames
  return(all_data)
}

fc = assemble_matrix(foldchange_list)
padj = assemble_matrix(padj_list, def = 1)
maxes = assemble_matrix(max_list)

padj = -log10(padj+10^-9)
maxes = log2(maxes + 1)

save(fc, padj, maxes, comp_from, comp_to, file = 'deseq_data.save')
