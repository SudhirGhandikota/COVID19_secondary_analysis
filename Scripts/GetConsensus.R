if (!require("biomaRt")) BiocManager::install("biomaRt")
library(optparse)
library(biomaRt)
option_list = list(
  make_option(c("-f", "--files"), type="character", default=NULL, 
              help="result RUVSeq files(comma-separated)", metavar="character"),
  make_option(c("-t", "--outpath"), type="character", default='/', 
              help="Outpath", metavar="character"),
  make_option(c("-o", "--org_assemblies"), type="character", default="GRCh38.p13", 
              help="Ensembl assemblies associated with each result file(comma-separated).
              https://uswest.ensembl.org/info/about/species.html", 
              metavar="character"),
  make_option(c("-l", "--logFC"), type="double", default="0.5", 
              help="logFC value for threshold", metavar="character"),
  make_option(c("-p", "--pvalue"), type="double", default="0.05", 
              help="logFC value for threshold", metavar="character"),
  make_option(c("-k", "--k"), type="numeric", default="2", 
              help="Minimum number of studies the gene has to be differentially expressed (atleast 2 or more)",
              metavar="character")
)
get_human_orthologs = function(assembly, genes){
  #assembly = 'ChlSab1.1'
  #genes = genes
  ensembl = useMart("ensembl")
  all_datasets = listDatasets(ensembl) 
  dataset = all_datasets[all_datasets$version==assembly,]$dataset
  if(length(dataset)==0){
    print(paste("Unable to recognize the provided Ensembl assembly:", assembly, 
                "Using the original symbols"))
    return(genes)
  }
  ensembl = useDataset(dataset,mart = ensembl)
  homologs = getBM(attributes=c('hsapiens_homolog_associated_gene_name'), 
                   filters = 'external_gene_name', values = genes, mart = ensembl)
  homologs = unique(as.vector(homologs[,'hsapiens_homolog_associated_gene_name']))
  return(homologs)
}

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
print(opt)
files = opt$files
files = unlist(strsplit(files, ","))

org_ids = opt$org_assemblies
# assuming all studies are based on human samples if ensembl assemblies aren't provided
if(org_ids == 'GRCh38.p13'){
  print(paste("No species assemblies have been provided.",
  "Assuming all symbols belong to the same species"))
  org_ids = rep('GRCh38.p13', length(files))
}else{
  org_ids = unlist(strsplit(org_ids, ","))
  if(length(org_ids) != length(files)){
    print(paste("Number of species assemblies do not match the number of input files provided.",
                "Please try again", sep=""))
    quit()
  }
  
}

logfc = as.numeric(opt$logFC)
pval = as.numeric(opt$pvalue)

up_degs = list()
down_degs = list()
system.time(for(i in 1:length(files)){
  infile = files[i]
  org_id = org_ids[i]
  deg_result = read.table(infile, sep="\t", header = T, row.names = 1, check.names = F)
  degs_up = row.names(deg_result[deg_result$logFC>=logfc & deg_result$FDR<pval,])
  degs_dn = row.names(deg_result[deg_result$logFC<=-logfc & deg_result$FDR<pval,])
  if(org_id != 'GRCh38.p13'){
    degs_up = get_human_orthologs(org_id, degs_up)
    degs_dn = get_human_orthologs(org_id, degs_dn)
    #degs_up = as.vector(homologene(degs_up, inTax = org_id, outTax = 9606)[,'9606'])
  }
  up_degs[[i]] = degs_up
  down_degs[[i]] = degs_up
  write.table(degs_up, paste(dirname(infile), "/", "Upregulated.txt", sep=""),
              sep="\t", row.names = F, quote = F, col.names = F)
  write.table(degs_dn, paste(dirname(infile), "/", "Downregulated.txt", sep=""),
              sep="\t", row.names = F, quote = F, col.names = F)
})

# consensus signature
k = opt$k
if(k==length(files)){
  cons_up = Reduce(intersect, up_degs)
  cons_dn = Reduce(intersect, down_degs)
}else{
  combs = combn(1:length(files), 2)
  cons_up = unique(unlist(apply(combs, 2, function(x) intersect(up_degs[[x[1]]], up_degs[[x[2]]]))))
  cons_dn = unique(unlist(apply(combs, 2, function(x) intersect(down_degs[[x[1]]], down_degs[[x[2]]]))))
}

outpath = opt$outpath
write.table(cons_up, paste(outpath,"/", "Consensus_up.txt", sep=""),
            row.names = F, col.names = F, quote = F)
write.table(cons_dn, paste(outpath,"/", "Consensus_down.txt", sep=""),
            row.names = F, col.names = F, quote = F)


