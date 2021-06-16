library(igraph)
library(MCL)
library(mclust)
library(optparse)
library(stringi)
library(dplyr)

option_list = list(
  make_option(c("-d", "--deg_file"), type="character", default=NULL, 
              help="File containing consensus transcriptome (and optional virus-host proteome) one symbol for each line ",
              metavar="character"),
  make_option(c("-p", "--PPI_file"), type="character", default=NULL, 
              help="Optional file containing human PPI network", metavar="character"),
  make_option(c("-f", "--filter"), type="character", default="combined_score >= 900 or experimental >= 700", 
              help="Filtering condition (optional) ", metavar="character"),
  make_option(c("-i", "--inflation_value"), type="character", default="2.5", 
              help="The inflation parameter for the MCL clustering (default = 2.5) ", metavar="character"),
  make_option(c("-m", "--max_iter"), type="character", default="100", 
              help="The inflation parameter for the MCL clustering ", metavar="character"),
  make_option(c("-o", "--outpath"), type="character", default='/', 
              help="Outpath", metavar="character")
)

start_time = Sys.time()
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
print(opt)
if(length(opt$deg_file)==0){
  print(paste("File containing conserved DEGs is missing. Please try again"))
  quit()
}
if(length(opt$PPI_file)==0){
  print(paste("PPI file is not provided. Using the default PPI file from STRING(v11)"))
  opt$PPI_file = "../input_data/other data/filtered_PPI.txt"
}

# # optional script to implement MCL clustering
# 
# # reading conserved signatures + virus-host proteome (optional) -> differentially expressed in at least 2 out 3 studies
degs = read.table(opt$deg_file, sep="\t", header = F)
degs = as.vector(degs$V1)
print(paste("Number of DEGs:", length(degs), sep=" "))

# # PPI links
PPI_data = read.table(opt$PPI_file, sep="\t", header = T)
print(paste("*****", "PPI Links", "*****"))
print(head(PPI_data))
print("Dimension (before filtering)")
print(dim(PPI_data))
condition = opt$filter
condition = stri_replace_all_fixed(condition, c(" OR ", " or ", " AND ", " and "), 
                                   c(" | ", " | ", " & ", " & "), vectorize_all = F)
PPI_data = PPI_data %>% filter(eval(parse(text=condition)))
print("Dimension (after filtering)")
print(dim(PPI_data))

# # filtering PPI among DEGs
PPI_data_fil = PPI_data[(PPI_data$protein1 %in% degs) & (PPI_data$protein2 %in% degs),]
PPI_data_fil = PPI_data_fil[!duplicated(t(apply(PPI_data_fil,1,sort))),]
PPI_data_fil = PPI_data_fil[!PPI_data_fil$protein1==PPI_data_fil$protein2,]

# # constructing a graph based on the interactome
PPI_data_fil$combined_score = NULL
PPI_data_fil$experimental = NULL
graph_fil = graph_from_data_frame(PPI_data_fil, directed = F)
# # nodes in the interactome
nodes_fil = names(V(graph_fil))
print(paste("Number of nodes in the final interactome:", length(nodes_fil), 
      "Weighted:", is.weighted(graph_fil)))

# # MCL clustering step
adj_fil = as_adjacency_matrix(graph_fil)
inflation = as.numeric(opt$inflation_value)
max_iter = as.numeric(opt$max_iter)
print(paste("*****", "MCL Clustering", "*****"))
print(paste("Parameters:", "Inflation:", inflation, "Maxiters:", max_iter))
mcl_results = mcl(adj_fil, inflation = inflation, addLoops = T, max.iter = max_iter)
mod_score = modularity(graph_fil, as.vector(mcl_results$Cluster+1))
print(paste("Number of clusters:", length(unique(mcl_results$Cluster))))

# # final clustering results
clustering_results = as.data.frame(cbind(nodes_fil, mcl_results$Cluster))
names(clustering_results) = c("Gene", "MCL_Clustering")
outfile = paste(opt$outpath,"/", "MCL_Clusters.txt", sep="")
write.table(clustering_results, outfile, sep = "\t", row.names = F, quote = F)
end_time = Sys.time()
print(end_time - start_time)