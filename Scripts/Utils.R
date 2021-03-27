# this method computes enrichments of a given set of scRNA-seq marker genes
# among the genes from a specific module
module_enrichment = function(markers_df, mod_genes, mod_name){
  cells = unique(as.vector(markers_df$cluster))
  identities = vector()
  num_tests = length(cells)
  cells_col = vector()
  cell_cnts_col = vector()
  pvals_col = vector()
  genes_col = vector()
  counts_col = vector()
  for(i in 1:length(cells)){
    cell = cells[i]
    cell_genes = as.vector(markers_df[markers_df$cluster==cell,]$gene)
    cell_identity = unique(as.vector(markers_df[markers_df$cluster==cell,]$Cell_Identity))
    cells_col = c(cells_col,cell)
    identities = c(identities, cell_identity)
    cell_cnts_col = c(cell_cnts_col, length(cell_genes))
    overlap = intersect(cell_genes,mod_genes)
    genes_col = c(genes_col,paste(overlap,collapse=","))
    a = length(overlap)
    counts_col = c(counts_col,a)
    b = length(mod_genes) - a
    c = length(cell_genes) - a
    d = 20000 - (a + b + c)
    contin_table = matrix(c(a,b,c,d), ncol = 2, byrow = T)
    p = fisher.test(contin_table)$p.value
    pvals_col = c(pvals_col,p)
  }
  pTable = as.matrix(p.adjust(pvals_col, method = 'BH'))
  pTable = -log10(pTable)
  CountTbl = as.matrix(counts_col)
  dimnames(CountTbl) = list(cells, mod_name)
  dimnames(pTable) = list(cells, mod_name)
  enrichment_results = data.frame(cbind(cells_col, identities, cell_cnts_col, counts_col,
                                        pvals_col, p.adjust(pvals_col, method = 'BH'),genes_col))
  colnames(enrichment_results) = c("Cell Type","Cell_Identity", "Number of Markers",
                                   "Overlap_Count", "Pval", "qval","Genes")
  enrichment_results$logp = -log10(as.numeric(as.vector(enrichment_results$qval)))
  enrichment_results$DEG = rep(mod_name,dim(enrichment_results)[1])
  enrichment_results$DEG_Size = rep(length(mod_genes),dim(enrichment_results)[1])
  return(list("enrichments" = enrichment_results, "CountTable" = CountTbl, "pTable" = pTable))
}

# phenotype enrichments within a given geneset
trait_enrichment = function(trait_df, mod_genes, mod_name, bg_cnt = 20000){
  traits = unique(as.vector(trait_df$Trait))
  num_tests = length(traits)
  traits_col = vector()
  trait_cnts_col = vector()
  pvals_col = vector()
  genes_col = vector()
  counts_col = vector()
  for(i in 1:length(traits)){
    trait = traits[i]
    trait_genes = unique(as.vector(trait_df[trait_df$Trait==trait,]$Gene))
    traits_col = c(traits_col,trait)
    trait_cnts_col = c(trait_cnts_col, length(trait_genes))
    overlap = intersect(trait_genes,mod_genes)
    genes_col = c(genes_col,paste(overlap,collapse=","))
    a = length(overlap)
    counts_col = c(counts_col,a)
    b = length(mod_genes) - a
    c = length(trait_genes) - a
    d = bg_cnt - (a + b + c)
    contin_table = matrix(c(a,b,c,d), ncol = 2, byrow = T)
    p = fisher.test(contin_table)$p.value
    pvals_col = c(pvals_col,p)
  }
  pTable = as.matrix(p.adjust(pvals_col, method = 'BH'))
  pTable = -log10(pTable)
  CountTbl = as.matrix(counts_col)
  dimnames(CountTbl) = list(traits, mod_name)
  dimnames(pTable) = list(traits, mod_name)
  enrichment_results = data.frame(cbind(traits_col, trait_cnts_col, counts_col,
                                        pvals_col, p.adjust(pvals_col, method = 'BH'),genes_col))
  colnames(enrichment_results) = c("Trait", "Number of Genes",
                                   "Overlap_Count", "Pval", "qval","Genes")
  enrichment_results$logp = -log10(as.numeric(as.vector(enrichment_results$qval)))
  enrichment_results$DEG = rep(mod_name,dim(enrichment_results)[1])
  enrichment_results$DEG_Size = rep(length(mod_genes),dim(enrichment_results)[1])
  return(list("enrichments" = enrichment_results, "CountTable" = CountTbl, "pTable" = pTable))
}

# function to plot enrichment (labeled) heatmaps using WGCNA package
# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/
plot_enrichments = function(pTable, countTable, ModTotals1,
                            ModTotals2, ylabels, xlabels, mar,
                            xangle = 45, cex_text = 1.0, cex_lab = 1.0){
  pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)]);
  pTable[pTable>50 ] = 50
  
  # Actual plotting
  sizeGrWindow(12,6 )
  par(mfrow=c(1,1))
  par(mar=mar, bg=NA)#c(bottom,left,top,right)
  # Use function labeledHeatmap to produce the color-coded table with all the trimmings
  labeledHeatmap(Matrix = pTable,
                 yLabels = paste(ylabels,":",ModTotals1),
                 xLabelsAdj = 1.0,
                 xLabels = xlabels,
                 colorLabels = FALSE,
                 ySymbols = paste(ylabels,": ",ModTotals1),
                 xSymbols = xlabels,
                 textMatrix = countTable,
                 colors = blueWhiteRed(100)[50:100],
                 main = "",
                 xLabelsAngle = xangle,
                 cex.text = cex_text, cex.lab = cex_lab, setStdMargins = FALSE)
}

# to compute overlaps between two separate gene (DEG) lists
deg_overlaps = function(deg_list1, deg_list2, deg_names1, deg_names2, 
                        bg_cnt = 20000, alt='two.sided'){
  degs1 = vector()
  degs2 = vector()
  cnts1 = vector()
  cnts2 = vector()
  pvals_col = vector()
  overlaps_col = vector()
  counts_col = vector()
  #pTable = matrix(0, nrow = length(deg_list1), ncol = length(deg_list2))
  #CountTbl = matrix(0, nrow = length(deg_list1), ncol = length(deg_list2))
  for(i in 1:length(deg_list1)){
    deg1 = deg_names1[i]
    genes1 = deg_list1[[i]]
    for(j in 1:length(deg_list2)){
      deg2 = deg_names2[j]
      genes2 = deg_list2[[j]]
      degs1 = c(degs1, deg1)
      degs2 = c(degs2, deg2)
      cnts1 = c(cnts1, length(genes1))
      cnts2 = c(cnts2, length(genes2))
      overlap = intersect(genes1, genes2)
      overlaps_col = c(overlaps_col,paste(overlap,collapse=","))
      a = length(overlap)
      counts_col = c(counts_col,a)
      b = length(genes1) - a
      c = length(genes2) - a
      d = bg_cnt - (a + b + c)
      contin_table = matrix(c(a,b,c,d), ncol = 2, byrow = T)
      p = fisher.test(contin_table, alternative = alt)$p.value
      pvals_col = c(pvals_col,p)
    }
  }
  pTable = matrix(p.adjust(pvals_col, method = 'BH'), nrow = length(deg_list2), ncol = length(deg_list1))
  pTable = -log10(pTable)
  CountTbl = matrix(counts_col, nrow = length(deg_list2), ncol = length(deg_list1))
  dimnames(CountTbl) = list(deg_names2, deg_names1)
  dimnames(pTable) = list(deg_names2, deg_names1)
  enrichment_results = data.frame(cbind(degs1, degs2, cnts1, cnts2, counts_col, pvals_col, 
                                        p.adjust(pvals_col, method = 'BH'), overlaps_col))
  colnames(enrichment_results) = c("DEG1","DEG2","Size_DEG1", "Size_DEG2", "Overlap_Count",
                                   "Pval", "qval", "Genes")
  enrichment_results$logp = -log10(as.numeric(as.vector(enrichment_results$qval)))
  return(list("enrichments" = enrichment_results, "CountTable" = CountTbl, "pTable" = pTable))
}

# fetch broad cell-categories for a given list of cell types
map_identities = function(markers_df, cells){
  identities = c()
  for(cell in cells){
    identity = unique(as.vector(markers_df[markers_df$cluster==cell,]$Cell_Identity))
    identities = c(identities, identity)
  }
  return(identities)
}

# generates a log-pvalue matrix based on overlaps between one or more gene lists
# and a given set of scRNA-seq marker genes
generate_logps_degs = function(markers_df, degs, deg_names){
  enrichments = data.frame()
  logp = list()
  print(length(degs))
  print(deg_names)
  for(i in 1:length(degs)){
    results = module_enrichment(markers_df, degs[[i]], deg_names[i])
    enrichments = rbind(enrichments, results$enrichments)
    logp[[i]] = results$pTable
  }
  logp = do.call(cbind, logp)
  logp = cbind(logp, map_identities(markers_df, row.names(logp)))
  sig_cells = unique(as.vector(enrichments[as.numeric(as.vector(enrichments$qval)) <= 0.05,'Cell Type']))
  
  if(length(sig_cells)==1){
    logp = t(as.matrix(logp[sig_cells,], nrow = 1))
    row.names(logp) = sig_cells
  }
  else{
    logp = logp[sig_cells,]
  }
  return(logp)
  logp = logp[sig_cells,]
  return(logp)
}

#### Reading Markers Info ####
str_split_fixed = function(clusters, delim){
  results = c()
  for(cluster in clusters){
    results = c(results,strsplit(cluster, "_")[[1]][1])
  }
  return(results)
}
map_marker_type = function(row, fc_col){
  log_fc = as.numeric(row[[fc_col]][1])
  cluster = row[['cluster']][1]
  if(log_fc > 0)
    cluster = paste(cluster, "-","UP",sep = '')
  else
    cluster = paste(cluster, "-", "DN", sep = '')
  return(cluster)
}
read_markers = function(){
  # reading scRNA-seq markers from Travaglini et.al
  # https://www.nature.com/articles/s41586-020-2922-4
  all_lung_markers = read.table('input_data/Lung_Markers/lung_markers_Krasnow.txt', 
                                sep = '\t', header = T)
  all_lung_markers = all_lung_markers[,-c(7,8,9,10,11,12,13,16)]
  all_lung_markers = all_lung_markers[all_lung_markers$p_val_adj<=0.05,]
  all_lung_markers = all_lung_markers[all_lung_markers$avg_logFC>=0.5,]
  all_lung_markers$cluster = as.vector(all_lung_markers$cluster)
  
  # retaining only 10X markers 
  lung_markers_krasnow = all_lung_markers[!grepl("SS2", all_lung_markers$cluster),]
  lung_markers_krasnow = lung_markers_krasnow[order(lung_markers_krasnow$Cell_Identity),]
  
  # reading markers published by Adams et.al
  # https://advances.sciencemag.org/content/6/28/eaba1983
  lung_markers_kaminski = read.table('input_data/Lung_Markers/lung_markers_Kaminski.txt', 
                                     sep = '\t', header = T)
  lung_markers_kaminski$Cell_Identity = lung_markers_kaminski$Cell.Category
  lung_markers_kaminski$Cell.Category = NULL
  lung_markers_kaminski$cluster = lung_markers_kaminski$cellType
  lung_markers_kaminski$cellType = NULL
  
  lung_markers_kaminski = lung_markers_kaminski[lung_markers_kaminski$logFC >= 0.5,]
  lung_markers_kaminski = lung_markers_kaminski[order(lung_markers_kaminski$Cell_Identity),]
  
  # reading lung markers published by Habermann et.al
  # https://advances.sciencemag.org/content/6/28/eaba1972
  lung_markers_banovich = read.table('input_data/Lung_Markers/lung_markers_Banovich.txt', 
                                     sep = '\t', header = T)
  lung_markers_banovich = lung_markers_banovich[lung_markers_banovich$avg_logFC >= 0.5,]
  lung_markers_banovich = lung_markers_banovich[lung_markers_banovich$p_val <= 0.05,]
  lung_markers_banovich = lung_markers_banovich[order(lung_markers_banovich$Cell_Identity),]

  return(list(lung_markers_krasnow, lung_markers_kaminski, lung_markers_banovich))
}

