library(optparse)
#workingDir = dirname(getwd()) # set the working directory to the git repository
#setwd(workingDir)
source('Utils.R')
getwd()

option_list = list(
  make_option(c("-a", "--assoc_file"), type="character", default=NULL, 
              help="File containing consensus phenotype-genotype associations from PheGenI 
              (https://www.ncbi.nlm.nih.gov/gap/phegeni) ",
              metavar="character"),
  make_option(c("-p", "--p_value"), type="numeric", default=0.00001, 
              help="P-value threshold for filtering the associations)",
              metavar="character"),
  make_option(c("-g", "--min_genes"), type="numeric", default=5, 
              help="Minimum number of genes in a candidate cluster",
              metavar="character"),
  make_option(c("-c", "--cluster_file"), type="character", default=NULL, 
              help="Two-column file (Gene-ClusterID) containing module memberships", 
              metavar="character"),
  make_option(c("-r", "--remove_intergenic"), action="store_true", type="logical", default=NULL, 
              help="Flag to indicate the removal of intergenic associations", 
              metavar="character"),
  make_option(c("-o", "--outpath"), type="character", default='/', 
              help="Outpath", metavar="character")
)
start_time = Sys.time()
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
print(opt)

phegeni_file = opt$assoc_file
phegen_info = read.csv(phegeni_file, sep = "\t", header = T)
if("remove_intergenic" %in% opt){
  print("Removing intergenic associations...")
  phegen_info = phegen_info[!phegen_info$Context=='intergenic',]
}

threshold = opt$p_value
phegen_info = phegen_info[phegen_info$P.Value<opt$p_value,]

node_stats = read.csv(opt$cluster_file, sep = "\t")
# # removing unclustered genes
node_stats = node_stats[!node_stats$cluster == '',]
# # candidate modules => at least 5 genes
selected_clusters = names(table(node_stats$cluster)[table(node_stats$cluster)>=opts$min_genes])

bg_cnt = length(union(unique(phegen_info$Gene), unique(phegen_info$Gene.2)))
cluster_trait_enrichments = data.frame()
for(cid in selected_clusters){
  print(cid)
  genes = as.vector(node_stats[node_stats$cluster==cid,]$name)
  enrichments = trait_enrichment(phegen_info, genes, cid, bg_cnt = bg_cnt)$enrichments
  enrichments = enrichments[as.numeric(as.vector(enrichments$Overlap_Count))>0,]
  cluster_trait_enrichments = rbind(cluster_trait_enrichments, enrichments)
}
# filtering enrichments
cluster_trait_enrichments = cluster_trait_enrichments[as.numeric(as.vector(cluster_trait_enrichments$Overlap_Count))>0,]
outfile = paste(opt$outpath,"/", "module_PheGenI_enrichments.txt", sep="")
write.table(cluster_trait_enrichments, outfile, sep = "\t")

end_time = Sys.time()
print(end_time - start_time)