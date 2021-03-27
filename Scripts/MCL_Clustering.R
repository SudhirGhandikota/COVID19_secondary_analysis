library(igraph)
library(MCL)
library(mclust)

# optional script to implement MCL clustering
workingDir = '/path/to/git/repository/' # set the working directory
setwd(workingDir)

# reading conserved signatures -> differentially expressed in at least 2 out 3 studies
# Supplemental Table S2
cons_up = as.vector(read.table("input_data/SARS-CoV2_DEGS/Conserved_UP.txt", header = F)$V1)
cons_dn = as.vector(read.table("input_data/SARS-CoV2_DEGS/Conserved_DN.txt", header = F)$V1)
PPI_genes = read.table('input_data/SARS-CoV2_DEGs/336 PPI.txt')$V1

# PPI links from STRING(v11) filtered based on 
# total score >= 900 or experimental score >= 700
PPI_data = read.table("input_data/other data/filtered_PPI.txt",sep="\t", header = T)

# filtering SARS-CoV2-interactome
PPI_data_sars2 = PPI_data[(PPI_data$Gene1 %in% c(cons_up,cons_dn,PPI_genes))
                          & (PPI_data$Gene2 %in% c(cons_up,cons_dn,PPI_genes)),]
PPI_data_sars2 = PPI_data_sars2[!duplicated(t(apply(PPI_data_sars2,1,sort))),]
PPI_data_sars2 = PPI_data_sars2[!PPI_data_sars2$Gene1==PPI_data_sars2$Gene2,]
PPI_data_sars2$Exp_Score = NULL
# constructing a graph based on the interactome
graph_sars2 = graph_from_data_frame(PPI_data_sars2, directed = F)
# nodes in the interactome
nodes_sars2 = names(V(graph_sars2))
print(is.weighted(graph_sars2))

# MCL clustering step
adj_sars2 = as_adjacency_matrix(graph_sars2)
mcl_results = mcl(adj_sars2, inflation = 2.5, addLoops = T)
mod_score = modularity(graph_sars2, as.vector(mcl_results$Cluster+1))
length(mcl_results$Cluster)
length(nodes_sars2)

# final clustering results
clustering_results = as.data.frame(cbind(nodes_sars2, mcl_results$Cluster))
names(clustering_results) = c("Gene", "MCL_Clustering")
#write.table(clustering_results, "SARS-CoV-2-Cons_MCL_Clusters", sep = "\t", row.names = F, quote = F)