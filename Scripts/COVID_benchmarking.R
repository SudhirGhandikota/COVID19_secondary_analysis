#BiocManager::install('RUVSeq')
#BiocManager::install('devtools')
library(RUVSeq)
library(limma)
library(edgeR)
library(combinat)
library(ggplot2)
library(httr)
library(jsonlite)
library(stringr)
library(stringi)
library(doParallel)

# returhs differentially expressed genes
get_degs = function(count_data, NumberofControls, NumberofTreatments, 
                    Counts, NumberofSamples){
  filter <- apply(count_data, 1, function(x) length(x[x>Counts])>=NumberofSamples)
  filtered <- count_data[filter,]
  genes <- rownames(filtered)
  aa <- length(genes)

  Controlrep <- rep("Ctl",NumberofControls)
  Treatmentrep <- rep("Trt",NumberofTreatments)
  x <- as.factor(c(Controlrep,Treatmentrep))
  set <- newSeqExpressionSet(as.matrix(filtered), phenoData = data.frame(x, row.names=colnames(filtered)))
  
  ##upper quartile normalization
  set <- betweenLaneNormalization(set, which="upper")
  ##Differential expression analysis using Empirical RUV
  design <- model.matrix(~x, data=pData(set))
  y <- DGEList(counts=counts(set), group=x)
  y <- calcNormFactors(y, method="upperquartile")
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  fit <- glmFit(y, design)
  lrt <- glmLRT(fit, coef=2)

  topTags(lrt)
  result1 <- topTags(lrt, n=dim(y)[1]+1, adjust.method="BH", sort.by="logFC")
  result1 = result1$table
  degs_up = result1[(result1$logFC>=0.6) & (result1$FDR<0.05),]
  degs_dn = result1[(result1$logFC<=-0.6) & (result1$FDR<0.05),]
  return(list("Up" = degs_up, "Dn" = degs_dn))
}

# generated randomized counts by randomly permuting samples (columns)
randomize_counts = function(count_data){
  samples = colnames(count_data)
  labels = get_pheno_labels(samples)
  comp_labels = ifelse(labels==0, 1, 0)
  while(TRUE){
    random_samples = sample(samples)
    if(all(samples==random_samples))
      next
    ran_data = count_data[,random_samples]
    ran_labels = get_pheno_labels(random_samples)
    if(all(labels==ran_labels) | (all(comp_labels==ran_labels)))
      next
    return(ran_data)
  }
  
}
get_pheno_labels = function(samples){
  pheno_labels = c()
  for(name in samples){
    if(grepl("Control", name, fixed = T))
      pheno_labels = c(pheno_labels, 0)
    else if(grepl("Mock", name, fixed = T))
      pheno_labels = c(pheno_labels, 0)
    else
      pheno_labels = c(pheno_labels, 1)
  }
  return(pheno_labels)
}
get_overlap_pval = function(genes1, genes2){
  overlap = intersect(genes1, genes2)
  if(length(overlap)==0)
    return(c(0, 0, 1.0))
  a = length(overlap)
  b = length(genes1) - a
  c = length(genes2) - a
  d = 20000 - (a + b + c)
  contin_table = matrix(c(a,b,c,d), ncol = 2, byrow = T)
  p = fisher.test(contin_table, alternative = 'greater')$p.value
  jaccard_score = a/length(union(genes1, genes2))
  return(c(length(overlap), jaccard_score, p))
}

# generates a line plot of DEG counts
count_plots = function(results_list, num_perms, up_intercept, 
                       down_intercept, ylim=100, scale_val = 500, title='', subtitle='',
                       legend_bool = TRUE, filename = NA, width=800, height=700){
  
  up_cnt = unlist(lapply(results_list, '[[', 'Up_cnt'))
  dn_cnt = unlist(lapply(results_list, '[[', 'Dn_cnt'))
  df = data.frame("Trial" = seq(1:num_perms), "DEG_Count" = c(up_cnt, dn_cnt),
                  "Type" = c(rep("Up", num_perms), rep("Down", num_perms)))
  if(!is.na(filename))
    tiff(filename, res = 144, width = width, height = height, units = 'px')
  ylim = max(max(up_intercept, down_intercept), ylim)
  p = ggplot(df, aes(x=Trial, y=DEG_Count)) + 
    geom_line(aes(color = Type), size=0.3, linetype="dashed") + 
    scale_color_manual(values = c("blue", "darkred")) +
    labs(title=title, subtitle = subtitle) + lims(y = c(0, ylim)) + 
    scale_y_continuous(breaks=seq(0, ylim, scale_val)) + theme_bw() +
    theme(plot.title = element_text(size=16, face="bold", hjust = 0.5),
          plot.subtitle = element_text(size=12, hjust = 0.5),
          axis.title.x = element_text(size=14, face="bold"),
          axis.title.y = element_text(size=14, face="bold"),
          axis.text.x = element_text(size=13),
          axis.text.y = element_text(size=13),
          plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
  if(up_intercept>0)
    p = p + geom_hline(yintercept=up_intercept, linetype="solid", color = "red", size=0.5)
  if(down_intercept>0)
    p = p + geom_hline(yintercept=down_intercept, linetype="solid", color = "blue", size=0.5)
  if(legend_bool){
    p = p + theme(legend.title = element_text(size=12, face="bold"))
    p = p + theme(legend.text = element_text(size=12))
    p = p + theme(legend.background = element_rect(fill="lightgrey", size=0.5, linetype="solid"))
    p = p + theme(legend.position='bottom')
    p = p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  }
  print(p)
  if(!is.na(filename))
    dev.off()
}
other_plot = function(results_list, colname, num_perms, intercept, xlabel,
                      ylabel, title, subtitle='', log_transform = FALSE, filename=NA){
  results = c()
  for(res in results_list)
    results = c(results, as.vector(unlist(res[colname]))[1])
  results = as.numeric(results)
  if(log_transform)
    results = -log10(results)
  df = as.data.frame(cbind(seq(1:num_perms),results))
  colnames(df) = c(xlabel, ylabel)
  if(!is.na(filename))
    tiff(filename, res = 144, width = 800, height = 700, units = 'px')
  p = ggplot(df, aes_string(x=xlabel, y=ylabel)) + 
    geom_line(size=0.2, linetype="dashed") + 
    geom_hline(yintercept=intercept, linetype="solid", color = "red", size=0.5) +
    labs(title=title, subtitle = subtitle) + theme_bw() +
    theme(plot.title = element_text(size=16, face="bold", hjust = 0.5),
          plot.subtitle = element_text(size=14, hjust = 0.5),
          axis.title.x = element_text(size=14, face="bold"),
          axis.title.y = element_text(size=14, face="bold"),
          axis.text.x = element_text(size=13),
          axis.text.y = element_text(size=13),
          legend.position = "bottom", plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
  
  print(p)
  if(!is.na(filename))
    dev.off()
}
# PPI enrichment among a given gene set using STRING REST API
# (https://string-db.org/help/api/)
string_enc = function(genes){
  request_url = 'https://string-db.org/api/tsv-no-header/ppi_enrichment'
  params = list("identifiers" = paste(genes, collapse = "%0d"), "species" = 9606, 
                required_score=900)
  r <- content(POST(request_url, body = params), as = "text")
  r = strsplit(str_replace(r, "\n", ""), split="\t")[[1]]
  cc = r[4] # clustering coefficient
  pval = r[6] # PPI enrichment p-value
  if(pval == "0.0")
    pval = "1.0e-16"
  return(c(cc, pval))
}

workingDir = "/path/to/git/repository/"
setwd(workingDir)

# human samples
GSE147507_data = read.table('input_data/Count Data/GSE147507/Genes_Counts_Matrix_Human_78Samples.txt', 
                            sep = "\t", header = T)
GSE147507_genes = row.names(GSE147507_data) = GSE147507_data$Gene
GSE147507_data$Gene = NULL
GSE147507_data = GSE147507_data[,grepl("Calu3", names(GSE147507_data), fixed = T)]
# ordering samples to make sure the controls come first (for DEG analysis)
GSE147507_data = GSE147507_data[,sort(names(GSE147507_data), decreasing = T)]

# Ad5-hACE2-sensitized mice 
GSE150847_data = read.table('input_data/Count Data/GSE150847/GSE150847_Genes_Counts_Matrix.txt', sep = "\t", header = T)
ortho_info = read.table('input_data/Count Data/GSE150847/MouseToHumanOrthologsFromBiomart.txt', sep = "\t", header = T)
row.names(GSE150847_data) = as.vector(GSE150847_data$Gene)
GSE150847_genes = unique(as.vector(ortho_info[ortho_info$Gene.name %in% as.vector(GSE150847_data$Gene),]$Human.gene.name))
GSE150847_data$Gene = NULL

# African green monkey
GSE153940_data = read.table('input_data/Count Data/GSE153940/GSE153940_Genes_Counts_Matrix.txt', sep = "\t", header = T)
row.names(GSE153940_data) = GSE153940_data$Gene
GSE153940_data$Gene = NULL
GSE153940_gene_info = read.table('input_data/Count Data/GSE153940/GSE153940_Genes_TPM_Matrix_WithHumanOrtho.txt', 
                                 sep = "\t", header = T)
GSE153940_genes = unique(as.vector(GSE153940_gene_info$Gene))
GSE153940_gene_info = as.data.frame(cbind(GSE153940_gene_info$Gene, GSE153940_gene_info$HumanOrtho))
names(GSE153940_gene_info) = c("Symbol", "HumanOrtho")
row.names(GSE153940_gene_info) = GSE153940_gene_info$Symbol

# Conserved transcriptome - Table S2
cons_up = read.table('input_data/SARS-CoV2_DEGs/Conserved_UP.txt')$V1
cons_down = read.table('input_data/SARS-CoV2_DEGs/Conserved_DN.txt')$V1
cons_degs = c(cons_up, cons_down)
PPI_genes = read.table('input_data/SARS-CoV2_DEGs/336 PPI.txt')$V1

GSE147507_up = as.vector(read.table('input_data/Count Data/GSE147507/GSE147507_Up.txt', header = F)$V1)
GSE147507_dn = as.vector(read.table('input_data/Count Data/GSE147507/GSE147507_Dn.txt', header = F)$V1)
GSE147507_degs = union(GSE147507_up, GSE147507_dn)

GSE150847_up = as.vector(read.table('input_data/Count Data/GSE150847/GSE150847_Up.txt', header = F)$V1)
GSE150847_dn = as.vector(read.table('input_data/Count Data/GSE150847/GSE150847_Dn.txt', header = F)$V1)
GSE150847_degs = union(GSE150847_up, GSE150847_dn)

GSE153940_up = as.vector(read.table('input_data/Count Data/GSE153940/GSE153940_Up.txt', header = F)$V1)
GSE153940_dn = as.vector(read.table('input_data/Count Data/GSE153940/GSE153940_Dn.txt', header = F)$V1)
GSE153940_degs = union(GSE153940_up, GSE153940_dn)

numcores = detectCores()
registerDoParallel(numcores-2)

############## random permutations - all 3 studies (Figure S4A) ####################
num_perms = 1000
permutation_results = foreach(i=1:num_perms) %dopar% {
  print(paste("Perm:",i))
  GSE147507_ran_degs = get_degs(randomize_counts(GSE147507_data), 3, 3, 10, 3)
  GSE147507_ran_up = row.names(GSE147507_ran_degs$Up)
  GSE147507_ran_dn = row.names(GSE147507_ran_degs$Dn)
  GSE147507_ran_degs = union(GSE147507_ran_up, GSE147507_ran_dn)
  results_up = get_overlap_pval(GSE147507_up, GSE147507_ran_up)
  results_dn = get_overlap_pval(GSE147507_dn, GSE147507_ran_dn)
  results_GSE147507 = list("Up_cnt"=length(GSE147507_ran_up), "Up_overlap"=results_up[1],
                                    "Up_jscore"=results_up[2], "Up_pval"=results_up[3], 
                                    "Dn_cnt"=length(GSE147507_ran_dn), "Dn_overlap"=results_dn[1],
                                    "Dn_jscore"=results_dn[2], "Dn_pval"=results_dn[3])
  
  GSE150847_ran_degs = get_degs(randomize_counts(GSE150847_data), 3, 3, 10, 3)
  GSE150847_ran_up = row.names(GSE150847_ran_degs$Up)
  GSE150847_ran_up = unique(as.vector(ortho_info[ortho_info$Gene.name %in% 
                                                   GSE150847_ran_up,]$Human.gene.name))
  GSE150847_ran_dn = row.names(GSE150847_ran_degs$Dn)
  GSE150847_ran_dn = unique(as.vector(ortho_info[ortho_info$Gene.name %in% 
                                                   GSE150847_ran_dn,]$Human.gene.name))
  GSE150847_ran_degs = union(GSE150847_ran_up, GSE150847_ran_dn)
  results_up = get_overlap_pval(GSE150847_up, GSE150847_ran_up)
  results_dn = get_overlap_pval(GSE150847_dn, GSE150847_ran_dn)
  results_GSE150847 = list("Up_cnt"=length(GSE150847_ran_up), "Up_overlap"=results_up[1],
                                    "Up_jscore"=results_up[2], "Up_pval"=results_up[3], 
                                    "Dn_cnt"=length(GSE150847_ran_dn), "Dn_overlap"=results_dn[1],
                                    "Dn_jscore"=results_dn[2], "Dn_pval"=results_dn[3])
  
  GSE153940_ran_degs = get_degs(randomize_counts(GSE153940_data), 3, 2, 10, 2)
  GSE153940_ran_up = row.names(GSE153940_ran_degs$Up)
  GSE153940_ran_up = unique(as.vector(GSE153940_gene_info[GSE153940_gene_info$Symbol %in% 
                                                            GSE153940_ran_up,]$HumanOrtho))
  GSE153940_ran_dn = row.names(GSE153940_ran_degs$Dn)
  GSE153940_ran_dn = unique(as.vector(GSE153940_gene_info[GSE153940_gene_info$Symbol %in% 
                                                            GSE153940_ran_dn,]$HumanOrtho))
  GSE153940_ran_degs = union(GSE153940_ran_up, GSE153940_ran_dn)
  results_up = get_overlap_pval(GSE153940_up, GSE153940_ran_up)
  results_dn = get_overlap_pval(GSE153940_dn, GSE153940_ran_dn)
  results_GSE153940 = list("Up_cnt"=length(GSE153940_ran_up), "Up_overlap"=results_up[1],
                                    "Up_jscore"=results_up[2], "Up_pval"=results_up[3], 
                                    "Dn_cnt"=length(GSE153940_ran_dn), "Dn_overlap"=results_dn[1],
                                    "Dn_jscore"=results_dn[2], "Dn_pval"=results_dn[3])
  
  cons_ran_up = unique(c(intersect(GSE147507_ran_up, GSE150847_ran_up), 
              intersect(GSE150847_ran_up, GSE153940_ran_up),
              intersect(GSE147507_ran_up, GSE153940_ran_up)))
  
  cons_ran_dn = unique(c(intersect(GSE147507_ran_dn, GSE150847_ran_dn), 
                     intersect(GSE150847_ran_dn, GSE153940_ran_dn),
                     intersect(GSE147507_ran_dn, GSE153940_ran_dn)))
  
  cons_ran_degs = unique(c(cons_ran_up, cons_ran_dn))
  results_up = get_overlap_pval(cons_up, cons_ran_up)
  results_dn = get_overlap_pval(cons_down, cons_ran_dn)
  results_cons = list("Up_cnt"=length(cons_ran_up), "Up_overlap"=results_up[1],
                                    "Up_jscore"=results_up[2], "Up_pval"=results_up[3], 
                                    "Dn_cnt"=length(cons_ran_dn), "Dn_overlap"=results_dn[1],
                                    "Dn_jscore"=results_dn[2], "Dn_pval"=results_dn[3])
  return(list("GSE147507"=results_GSE147507, "GSE150847" = results_GSE150847, 
              "GSE153940"=results_GSE153940, "cons"=results_cons))
}

# collecting results for each study from all randomized trials
all_results_GSE147507 = lapply(permutation_results, '[[','GSE147507')
all_results_GSE150847 = lapply(permutation_results, '[[','GSE150847')
all_results_GSE153940 = lapply(permutation_results, '[[','GSE153940')
all_results_cons = lapply(permutation_results, '[[','cons')
#save(all_results_GSE147507, all_results_GSE150847, all_results_GSE153940, 
     #all_results_cons, file='RData/benchmarking_individial_degs.RData')

load(file='RData/benchmarking_individial_degs.RData')

df_names = names(all_results_GSE147507[[1]])
up_mean = round(mean(unlist(lapply(all_results_GSE147507, '[[', 'Up_cnt'))),1)
up_std = round(sd(unlist(lapply(all_results_GSE147507, '[[', 'Up_cnt'))),2)
dn_mean = round(mean(unlist(lapply(all_results_GSE147507, '[[', 'Dn_cnt'))),1)
dn_std = round(sd(unlist(lapply(all_results_GSE147507, '[[', 'Dn_cnt'))),2)
count_plots(all_results_GSE147507, num_perms, length(GSE147507_up), length(GSE147507_dn), scale_val = 100,
            subtitle = paste("Observed: ", "Up = ", length(GSE147507_up), "; Down = ", length(GSE147507_dn), 
                             "\n", "random: ", "Up = ", up_mean, ' ± ', up_std, "; Down = ", dn_mean, ' ± ', dn_std, sep=""),
            title="GSE147507", filename="GSE147507_DEG_counts.tiff", width = 700, height=1200)

results_GSE147507 = as.data.frame(t(stri_list2matrix(all_results_GSE147507)))
colnames(results_GSE147507) = df_names
results_GSE147507$Trial = as.vector(paste("Trial", seq(1:num_perms)))
results_GSE147507$Total_cnt = as.numeric(results_GSE147507$Up_cnt)+as.numeric(results_GSE147507$Dn_cnt)
View(results_GSE147507[,c('Trial', 'Up_cnt', 'Dn_cnt', 'Total_cnt')])

up_mean = round(mean(unlist(lapply(all_results_GSE150847, '[[', 'Up_cnt'))),1)
up_std = round(sd(unlist(lapply(all_results_GSE150847, '[[', 'Up_cnt'))),2)
dn_mean = round(mean(unlist(lapply(all_results_GSE150847, '[[', 'Dn_cnt'))),1)
dn_std = round(sd(unlist(lapply(all_results_GSE150847, '[[', 'Dn_cnt'))),2)
count_plots(all_results_GSE150847, num_perms, length(GSE150847_up), length(GSE150847_dn), scale_val = 100,
            subtitle=paste("Observed: ", "Up = ", length(GSE150847_up), "; Down = ", length(GSE150847_dn), 
                           "\n", "random: ", "Up = ", up_mean, ' ± ', up_std, "; Down = ", dn_mean, ' ± ', dn_std, sep=""),
            title="GSE150847", filename="GSE150847_DEG_counts.tiff", width=700, height = 1200)

results_GSE150847 = as.data.frame(t(stri_list2matrix(all_results_GSE150847)))
colnames(results_GSE150847) = df_names
results_GSE150847$Trial = as.vector(paste("Trial", seq(1:num_perms)))
results_GSE150847$Total_cnt = as.numeric(results_GSE150847$Up_cnt)+as.numeric(results_GSE150847$Dn_cnt)
View(results_GSE150847[,c('Trial', 'Up_cnt', 'Dn_cnt', 'Total_cnt')])

up_mean = round(mean(unlist(lapply(all_results_GSE153940, '[[', 'Up_cnt'))),1)
up_std = round(sd(unlist(lapply(all_results_GSE153940, '[[', 'Up_cnt'))),2)
dn_mean = round(mean(unlist(lapply(all_results_GSE153940, '[[', 'Dn_cnt'))),1)
dn_std = round(sd(unlist(lapply(all_results_GSE153940, '[[', 'Dn_cnt'))),2)

count_plots(all_results_GSE153940, num_perms, length(GSE153940_up), length(GSE153940_dn),scale_val = 50,
            subtitle=paste("Observed: ", "Up = ", length(GSE153940_up), "; Down = ", length(GSE153940_dn), 
                           "\n", "random: ", "Up = ", up_mean, ' ± ', up_std, "; Down = ", dn_mean, ' ± ', dn_std, sep=""),
            title="GSE153940", filename= "GSE153940_DEG_counts.tiff", width = 700, height = 1200)

results_GSE153940 = as.data.frame(t(stri_list2matrix(all_results_GSE153940)))
colnames(results_GSE153940) = df_names
results_GSE153940$Trial = as.vector(paste("Trial", seq(1:num_perms)))
results_GSE153940$Total_cnt = as.numeric(results_GSE153940$Up_cnt)+as.numeric(results_GSE153940$Dn_cnt)
View(results_GSE153940[,c('Trial', 'Up_cnt', 'Dn_cnt', 'Total_cnt')])

############## random permutations - keeping one constant and change the other 2 studies (Figure S4B) ############
num_perms = 1000
perm_results1 = foreach(i=1:num_perms) %dopar% {
  GSE150847_ran_degs = get_degs(randomize_counts(GSE150847_data), 3, 3, 10, 3)
  GSE150847_ran_up = row.names(GSE150847_ran_degs$Up)
  GSE150847_ran_up = unique(as.vector(ortho_info[ortho_info$Gene.name %in% GSE150847_ran_up,]$Human.gene.name))
  GSE150847_ran_dn = row.names(GSE150847_ran_degs$Dn)
  GSE150847_ran_dn = unique(as.vector(ortho_info[ortho_info$Gene.name %in% GSE150847_ran_dn,]$Human.gene.name))
  GSE150847_ran_degs = union(GSE150847_ran_up, GSE150847_ran_dn)
  
  GSE153940_ran_degs = get_degs(randomize_counts(GSE153940_data), 3, 2, 10, 2)
  GSE153940_ran_up = row.names(GSE153940_ran_degs$Up)
  GSE153940_ran_up = unique(as.vector(GSE153940_gene_info[GSE153940_gene_info$Symbol %in% GSE153940_ran_up,]$HumanOrtho))
  GSE153940_ran_dn = row.names(GSE153940_ran_degs$Dn)
  GSE153940_ran_dn = unique(as.vector(GSE153940_gene_info[GSE153940_gene_info$Symbol %in% GSE153940_ran_dn,]$HumanOrtho))
  GSE153940_ran_degs = union(GSE153940_ran_up, GSE153940_ran_dn)
  
  cons_ran_up = unique(c(intersect(GSE147507_up, GSE150847_ran_up), 
                         intersect(GSE150847_ran_up, GSE153940_ran_up),
                         intersect(GSE147507_up, GSE153940_ran_up)))
  
  cons_ran_dn = unique(c(intersect(GSE147507_dn, GSE150847_ran_dn), 
                         intersect(GSE150847_ran_dn, GSE153940_ran_dn),
                         intersect(GSE147507_dn, GSE153940_ran_dn)))
  
  cons_ran_degs = unique(c(cons_ran_up, cons_ran_dn))
  print(paste("Iter:",i,length(cons_ran_degs)))
  results_up = get_overlap_pval(cons_up, cons_ran_up)
  results_dn = get_overlap_pval(cons_down, cons_ran_dn)
  list("Up_cnt"=length(cons_ran_up), "Up_overlap"=results_up[1], "Up_jscore"=results_up[2], 
       "Up_pval"=results_up[3], "Dn_cnt"=length(cons_ran_dn), "Dn_overlap"=results_dn[1],
        "Dn_jscore"=results_dn[2], "Dn_pval"=results_dn[3])

}
perm_results2 = foreach(i=1:num_perms) %dopar% {
  
  GSE147507_ran_degs = get_degs(randomize_counts(GSE147507_data), 3, 3, 10, 3)
  GSE147507_ran_up = row.names(GSE147507_ran_degs$Up)
  GSE147507_ran_dn = row.names(GSE147507_ran_degs$Dn)
  GSE147507_ran_degs = union(GSE147507_ran_up, GSE147507_ran_dn)
  
  GSE150847_ran_degs = get_degs(randomize_counts(GSE150847_data), 3, 3, 10, 3)
  GSE150847_ran_up = row.names(GSE150847_ran_degs$Up)
  GSE150847_ran_up = unique(as.vector(ortho_info[ortho_info$Gene.name %in% GSE150847_ran_up,]$Human.gene.name))
  GSE150847_ran_dn = row.names(GSE150847_ran_degs$Dn)
  GSE150847_ran_dn = unique(as.vector(ortho_info[ortho_info$Gene.name %in% GSE150847_ran_dn,]$Human.gene.name))
  GSE150847_ran_degs = union(GSE150847_ran_up, GSE150847_ran_dn)
  
  cons_ran_up = unique(c(intersect(GSE147507_ran_up, GSE150847_ran_up), 
                         intersect(GSE150847_ran_up, GSE153940_up),
                         intersect(GSE147507_ran_up, GSE153940_up)))
  
  cons_ran_dn = unique(c(intersect(GSE147507_ran_dn, GSE150847_ran_dn), 
                         intersect(GSE150847_ran_dn, GSE153940_dn),
                         intersect(GSE147507_ran_dn, GSE153940_dn)))
  
  cons_ran_degs = unique(c(cons_ran_up, cons_ran_dn))
  print(paste("Iter:",i,length(cons_ran_degs)))
  results_up = get_overlap_pval(cons_up, cons_ran_up)
  results_dn = get_overlap_pval(cons_down, cons_ran_dn)
  list("Up_cnt"=length(cons_ran_up), "Up_overlap"=results_up[1],"Up_jscore"=results_up[2], 
       "Up_pval"=results_up[3], "Dn_cnt"=length(cons_ran_dn), "Dn_overlap"=results_dn[1],
      "Dn_jscore"=results_dn[2], "Dn_pval"=results_dn[3])
  
}
perm_results3 = foreach(i=1:num_perms) %dopar% {
  GSE147507_ran_degs = get_degs(randomize_counts(GSE147507_data), 3, 3, 10, 3)
  GSE147507_ran_up = row.names(GSE147507_ran_degs$Up)
  GSE147507_ran_dn = row.names(GSE147507_ran_degs$Dn)
  GSE147507_ran_degs = union(GSE147507_ran_up, GSE147507_ran_dn)
  
  GSE153940_ran_degs = get_degs(randomize_counts(GSE153940_data), 3, 2, 10, 2)
  GSE153940_ran_up = row.names(GSE153940_ran_degs$Up)
  GSE153940_ran_up = unique(as.vector(GSE153940_gene_info[GSE153940_gene_info$Symbol %in% GSE153940_ran_up,]$HumanOrtho))
  GSE153940_ran_dn = row.names(GSE153940_ran_degs$Dn)
  GSE153940_ran_dn = unique(as.vector(GSE153940_gene_info[GSE153940_gene_info$Symbol %in% GSE153940_ran_dn,]$HumanOrtho))
  GSE153940_ran_degs = union(GSE153940_ran_up, GSE153940_ran_dn)
  
  cons_ran_up = unique(c(intersect(GSE147507_ran_up, GSE150847_up), 
                         intersect(GSE150847_up, GSE153940_ran_up),
                         intersect(GSE147507_ran_up, GSE153940_ran_up)))
  
  cons_ran_dn = unique(c(intersect(GSE147507_ran_dn, GSE150847_dn), 
                         intersect(GSE150847_dn, GSE153940_ran_dn),
                         intersect(GSE147507_ran_dn, GSE153940_ran_dn)))
  
  cons_ran_degs = unique(c(cons_ran_up, cons_ran_dn))
  print(paste("Iter:",i,length(cons_ran_degs)))
  results_up = get_overlap_pval(cons_up, cons_ran_up)
  results_dn = get_overlap_pval(cons_down, cons_ran_dn)
  list("Up_cnt"=length(cons_ran_up), "Up_overlap"=results_up[1],"Up_jscore"=results_up[2], 
       "Up_pval"=results_up[3], "Dn_cnt"=length(cons_ran_dn), "Dn_overlap"=results_dn[1],
      "Dn_jscore"=results_dn[2], "Dn_pval"=results_dn[3])
  
}

num_trials = 3*num_perms
all_perm_results2 = c(perm_results1, perm_results2, perm_results2)
#save(all_perm_results2,file='RData/benchmarking_two_out_of_three.RData')

load(file='RData/benchmarking_two_out_of_three.RData')
up_mean = round(mean(unlist(lapply(all_perm_results2, '[[', 'Up_cnt'))),1)
up_std = round(sd(unlist(lapply(all_perm_results2, '[[', 'Up_cnt'))),2)
dn_mean = round(mean(unlist(lapply(all_perm_results2, '[[', 'Dn_cnt'))),1)
dn_std = round(sd(unlist(lapply(all_perm_results2, '[[', 'Dn_cnt'))),2)

count_plots(all_perm_results2, num_trials, length(cons_up), length(cons_down), scale_val=50, title= "Consensus-DEGs", 
            subtitle=paste("Observed: ", "Up = ", length(cons_up), "; Down = ", length(cons_down), 
                           "\n", "random: ", "Up = ", up_mean, ' ± ', up_std, "; Down = ", dn_mean, ' ± ', dn_std, sep=""),
            filename="Cons_DEG_counts_v3.tiff", width = 800, height = 1200)

results_cons = as.data.frame(t(stri_list2matrix(all_perm_results2)))
colnames(results_cons) = df_names
results_cons$Trial = as.vector(paste("Trial", seq(1:num_trials)))
results_cons$Total_cnt = as.numeric(results_cons$Up_cnt)+as.numeric(results_cons$Dn_cnt)
View(results_cons[,c('Trial', 'Up_cnt', 'Dn_cnt', 'Total_cnt')])

############## randomized experiments - consensus DEGs ###################
# reading individual SARS-CoV-2 DEGS from two in vitro (Calu-3 and Vero E6 cells)
# and one in vivo (Ad5-hACE2-sensitized mice) study. Also included are DEGs from
# nasopharyngeal swabs of human COVID-19 patients (GSE152075)
# Supplemental Table S1
calu3_up = as.vector(read.table('input_data/SARS-CoV2_DEGS/calu3_CoV_2_Up.txt', 
                                     sep = "\t", header = F)$V1)
vero_up = as.vector(read.table('input_data/SARS-CoV2_DEGS//Vero_E6_CoV_2_Up.txt', 
                                    sep = "\t", header = F)$V1)
mm_up = as.vector(read.table('input_data/SARS-CoV2_DEGS//Mm_SARS_CoV_2_up.txt', 
                                  sep = "\t", header = F)$V1)
calu3_dn = as.vector(read.table('input_data/SARS-CoV2_DEGS//calu3_CoV_2_DN.txt', 
                                     sep = "\t", header = F)$V1)
vero_dn = as.vector(read.table('input_data/SARS-CoV2_DEGS//Vero_E6_CoV_2_DN.txt', 
                                    sep = "\t", header = F)$V1)
mm_dn = as.vector(read.table('input_data/SARS-CoV2_DEGS//Mm_SARS_CoV_2_DN.txt', 
                                  sep = "\t", header = F)$V1)

num_perms=1000
all_results_cons = foreach(i=1:num_perms) %dopar% {
  # randomly picking DEGs for each study
  calu3_ran_up = sample(GSE147507_genes, length(calu3_up))
  calu3_ran_dn = sample(GSE147507_genes, length(calu3_dn))
  
  vero_ran_up = sample(GSE153940_genes, length(vero_up))
  vero_ran_dn = sample(GSE153940_genes, length(vero_dn))
  
  mm_ran_up = sample(GSE150847_genes, length(mm_up))
  mm_ran_dn = sample(GSE150847_genes, length(mm_dn))
  
  # identifying consensus DEGs - at least 2 out 3
  cons_ran_up = unique(c(intersect(calu3_ran_up, vero_ran_up), 
                         intersect(vero_ran_up, mm_ran_up),
                         intersect(calu3_ran_up, mm_ran_up)))
  
  cons_ran_dn = unique(c(intersect(calu3_ran_dn, vero_ran_dn), 
                         intersect(vero_ran_dn, mm_ran_dn),
                         intersect(calu3_ran_dn, mm_ran_dn)))
  
  cons_ran_degs = unique(c(cons_ran_up, cons_ran_dn))
  # overlaps between random and actual DEGs
  results_up = get_overlap_pval(cons_up, cons_ran_up)
  results_dn = get_overlap_pval(cons_down, cons_ran_dn)
  results_cons = get_overlap_pval(cons_degs, cons_ran_degs)
  # PPI enrichments among random consensus DEGs
  results_str_enc = string_enc(c(cons_ran_degs, PPI_genes))
  list("Up_cnt"=length(cons_ran_up), "Up_overlap"=results_up[1], "Up_jscore"=results_up[2], 
       "Up_pval"=results_up[3],"Dn_cnt"=length(cons_ran_dn), 'Dn_overlap'=results_dn[1],
      'Dn_jscore'=results_dn[2], 'Dn_pval'=results_dn[3], 'Cons_deg_cnt'=length(cons_ran_degs), 
      "Cons_overlap"=results_cons[1], 'Cons_jscore'=results_cons[2], 'Cons_overlap'=results_cons[3],
      "Clustering_coef"=results_str_enc[1], "PPI_pval"=results_str_enc[2])
}

#save(all_results_cons, file='RData/benchmarking_consensus_degs.RData')

load('RData/benchmarking_consensus_degs.RData')

up_mean = round(mean(unlist(lapply(all_results_cons, '[[', 'Up_cnt'))),1)
up_std = round(sd(unlist(lapply(all_results_cons, '[[', 'Up_cnt'))),2)
dn_mean = round(mean(unlist(lapply(all_results_cons, '[[', 'Dn_cnt'))),1)
dn_std = round(sd(unlist(lapply(all_results_cons, '[[', 'Dn_cnt'))),2)

count_plots(all_results_cons, num_perms, 833, 634, scale_val = 50, title="Consensus-DEGs", 
            subtitle=paste("Observed: ", "Up = ", length(cons_up), "; Down = ", length(cons_down), 
                           "\n", "random: ", "Up = ", up_mean, ' ± ', up_std, "; Down = ", dn_mean, ' ± ', dn_std, sep=""),
            legend_bool = T, filename = "Cons_DEG_counts.tiff", width = 800, height = 1200)
other_plot(all_results_cons, 'PPI_pval', num_perms, -log10(1.0e-16), 
           "Trials","logP", 'PPI Enrichments', "Observed: < 1.0e-16", 
           log_transform = TRUE, filename="PPI_enc.tiff")

results_cons = as.data.frame(t(stri_list2matrix(all_results_cons)))
colnames(results_cons) = names(all_results_cons[[1]])
results_cons$Trial = as.vector(paste("Trial", seq(1:num_perms)))
results_cons$Total_cnt = as.numeric(results_cons$Up_cnt)+as.numeric(results_cons$Dn_cnt)
View(results_cons[,c('Trial', 'Up_cnt', 'Dn_cnt', 'Total_cnt', 'PPI_pval')])

############## randomized experiments - consensus DEGs and PPI enrichments ##############
# these experiments need additional libraries 'igraph', 'MCL'
library(igraph)
library(MCL)
library(mclust)

# PPI links from STRING(v11) filtered based on 
# total score >= 900 or experimental score >= 700
PPI_data = read.table("input_data/other data/filtered_PPI.txt",sep="\t", header = T)

num_perms = 1000
all_results_cons2 = foreach(i=1:num_perms) %dopar%{
  # generating random consensus DEGs (from human gene symbols)
  cons_ran_degs = sample(GSE147507_genes, 1803)
  results_str_enc = string_enc(cons_ran_degs)
  
  # filtering PPI data based on random DEGs
  PPI_data_ran = PPI_data[(PPI_data$Gene1 %in% cons_ran_degs)
                          & (PPI_data$Gene2 %in% cons_ran_degs),]
  PPI_data_ran = PPI_data_ran[!duplicated(t(apply(PPI_data_ran,1,sort))),]
  PPI_data_ran = PPI_data_ran[!PPI_data_ran$Gene1==PPI_data_ran$Gene2,]
  
  # igraph for random interactome
  graph_ran = graph_from_data_frame(PPI_data_ran, directed = F)
  nodes_ran = names(V(graph_ran))
  # MCL clustering
  adj_ran = as_adjacency_matrix(graph_ran)
  mcl_results_ran = mcl(adj_ran, inflation = 2.5, addLoops = F)
  if(grepl("Error",mcl_results_ran))
    modularity_ran = 0.0
  else{
    clustering_ran = mcl_results_ran$Cluster+1
    modularity_ran = modularity(graph_ran, as.vector(clustering_ran))
  }
  list("Clustering_coef"=results_str_enc[1], "PPI_pval"=results_str_enc[2],
        "Modularity"=modularity_ran)
}

#save(all_results_cons2, file='RData/benchmarking_consensus_and_ppi.RData')

load(file='RData/benchmarking_consensus_and_ppi.RData')
other_plot(all_results_cons2, 'PPI_pval', num_perms,-log10(1.0e-16), 
           "Trials","logP", 'PPI Enrichments', "Observed: < 1.0e-16",
           log_transform = TRUE, filename="PPI_enc2.tiff")
other_plot(all_results_cons2, 'Clustering_coef', num_perms, 0.42, 
           "Trials","CC", 'Clustering Coefficient', "Observed ≈ 0.43",
           log_transform = FALSE, filename="PPI_cc2.tiff")

results_cons2 = as.data.frame(t(stri_list2matrix(all_results_cons2)))
colnames(results_cons2) = names(all_results_cons2[[1]])
results_cons2$Trial = as.vector(paste("Trial", seq(1:num_perms)))
results_cons2$Clustering_coef = round(as.numeric(as.vector(results_cons2$Clustering_coef)), 2)
results_cons2$Modularity = round(as.numeric(as.vector(results_cons2$Modularity)), 2)
View(results_cons2[,c('Trial', 'Clustering_coef', 'PPI_pval')])
