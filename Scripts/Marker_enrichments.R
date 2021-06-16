library(optparse)
#workingDir = '/path/to/git/repository/' # set the working directory
#setwd(workingDir)

marker_enrichment = function(markers_df, mod_genes, mod_name){
  cells = unique(as.vector(markers_df$cell))
  num_tests = length(cells)
  cells_col = vector()
  cell_cnts_col = vector()
  pvals_col = vector()
  genes_col = vector()
  counts_col = vector()
  for(i in 1:length(cells)){
    cell = cells[i]
    cell_genes = as.vector(markers_df[markers_df$cell==cell,]$gene)
    cells_col = c(cells_col,cell)
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
  enrichment_results = data.frame(cbind(cells_col, cell_cnts_col, counts_col,
                                        pvals_col, p.adjust(pvals_col, method = 'BH'),genes_col))
  colnames(enrichment_results) = c("Cell Type","Number of Markers",
                                   "Overlap_Count", "Pval", "qval","Genes")
  enrichment_results$logp = -log10(as.numeric(as.vector(enrichment_results$qval)))
  enrichment_results$DEG = rep(mod_name,dim(enrichment_results)[1])
  enrichment_results$DEG_Size = rep(length(mod_genes),dim(enrichment_results)[1])
  return(list("enrichments" = enrichment_results, "CountTable" = CountTbl, "pTable" = pTable))
}

option_list = list(
  make_option(c("-a", "--marker_file"), type="character", default=NULL, 
              help="File containing scRNA-seq marker information 
              (cell type, gene symbol/ID, logFC, p-value)",
              metavar="character"),
  make_option(c("-p", "--p_value"), type="numeric", default=0.05, 
              help="P-value threshold for filtering the marker associations)",
              metavar="character"),
  make_option(c("-l", "--logFC"), type="numeric", default=0.5, 
              help="P-value threshold for filtering the marker associations)",
              metavar="character"),
  make_option(c("-c", "--cluster_file"), type="character", default=NULL, 
              help="Two-column file (Gene-ClusterID) containing module memberships", 
              metavar="character"),
  make_option(c("-o", "--outpath"), type="character", default='/', 
              help="Outpath", metavar="character")
)

start_time = Sys.time()
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
print(opt)

########## Marker enrichment in gene modules from SARS-CoV-2 DEG interactome - Figure 4 ########

# reading gene clusters from conserver DEG-PPI interactome in SARS-CoV-2.
# Steps followed to obtain the gene modules:
  # 1. conserved DEGs are combined with SARS-CoV2-Human PPI genes (336 genes)
  # 2. The combined gene set was uploaded to STRING (https://string-db.org/) to identify SARS-CoV-2 interactome
  # 3. The combined interactome is analysed using the Cytoscape (https://cytoscape.org/) tool.
  # 4. clusterMaker2 (http://www.rbvi.ucsf.edu/cytoscape/clusterMaker2/) plugin was used to implement MCL clustering

node_stats = read.csv(opt$cluster_file, sep = "\t")
# # removing unclustered geness
node_stats = node_stats[!node_stats$cluster == '',]
# # candidate modules => at least 5 genes
selected_clusters = names(table(node_stats$cluster)[table(node_stats$cluster)>4])

logfc = as.numeric(opt$logFC)
pval = as.numeric(opt$p_value)

markers_df = read.table(opt$marker_file, sep = '\t', header = T)
markers_df = markers_df[markers_df$logFC >= logfc,]
markers_df = markers_df[markers_df$pval_adj <= pval,]

cluster_genes = list()
for(i in selected_clusters){
  genes = as.vector(node_stats[node_stats$cluster==i,]$gene)
  cluster_genes[[i]] = genes
}

enc_results = data.frame()

for(i in 1:length(cluster_genes)){
  enc_results = rbind(enc_results, marker_enrichment(markers_df, cluster_genes[[i]], 
                                                selected_clusters[i])$enrichments)
}
enc_results = enc_results[as.numeric(as.vector(enc_results$Overlap_Count))>0,]
outfile = paste(opt$outpath,"/", "module_marker_enrichments.txt", sep="")
print(outfile)
write.table(enc_results_combined, outfile, sep = "\t", row.names = F, quote = F)

end_time = Sys.time()
print(end_time - start_time)