library(WGCNA)
workingDir = '/path/to/git/repository/' # set the working directory
setwd(workingDir)
source('Scripts/Utils.R')

###################### Reading SARS-CoV-2 DEGs #####################

# reading individual SARS-CoV-2 DEGS from two in vitro (Calu-3 and Vero E6 cells)
# and one in vivo (Ad5-hACE2-sensitized mice) study. Also included are DEGs from
# nasopharyngeal swabs of human COVID-19 patients (GSE152075)
# Supplemental Table S1
calu3_cov2_up = as.vector(read.table('input_data/SARS-CoV2_DEGS/calu3_CoV_2_Up.txt', 
                                                 sep = "\t", header = F)$V1)
vero_cov2_up = as.vector(read.table('input_data/SARS-CoV2_DEGS//Vero_E6_CoV_2_Up.txt', 
                                                sep = "\t", header = F)$V1)
mm_cov2_up = as.vector(read.table('input_data/SARS-CoV2_DEGS//Mm_SARS_CoV_2_up.txt', 
                                              sep = "\t", header = F)$V1)
human_cov2_up = as.vector(read.table('input_data/SARS-CoV2_DEGS//COVID_human_UP.txt', 
                                                 sep = "\t", header = F)$V1)
calu3_cov2_dn = as.vector(read.table('input_data/SARS-CoV2_DEGS//calu3_CoV_2_DN.txt', 
                                                 sep = "\t", header = F)$V1)
vero_cov2_dn = as.vector(read.table('input_data/SARS-CoV2_DEGS//Vero_E6_CoV_2_DN.txt', 
                                                sep = "\t", header = F)$V1)
mm_cov2_dn = as.vector(read.table('input_data/SARS-CoV2_DEGS//Mm_SARS_CoV_2_DN.txt', 
                                               sep = "\t", header = F)$V1)
human_cov2_dn = as.vector(read.table('input_data/SARS-CoV2_DEGS//COVID_human_Dn.txt', 
                                                  sep = "\t", header = F)$V1)

########### Computing overlaps among DEGs from the three input SARS-CoV-2 studies - Figure 2B ##########
degs = list(calu3_cov2_up, vero_cov2_up, mm_cov2_up, 
            calu3_cov2_dn, vero_cov2_dn, mm_cov2_dn)
deg_names = c("Calu-3_CoV-2_UP", "Vero_E6_CoV-2_UP", "Mm_SARS-CoV-2_UP", 
              "Calu-3_CoV-2_DN",  "Vero_E6_CoV-2_DN", "Mm_SARS-CoV-2_DN")

results = deg_overlaps(degs, degs, deg_names, deg_names, alt = 'greater')
pTable = results$pTable
countTable = results$CountTable
# truncating negative logP value matrix
pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)]);
pTable[pTable>50 ] = 50

# Figure 2B
plot_enrichments(pTable, countTable, c(), c(), deg_names, deg_names,
                 mar=c(8,11,1,1), xangle = 45, cex_text = 1.1, cex_lab = 1.1)

results$enrichments = results$enrichments[as.numeric(as.vector(results$enrichments$Overlap_Count)) > 0,]
results$enrichments = results$enrichments[!results$enrichments$DEG1 == results$enrichments$DEG2,]

#write.table(results$enrichments, 'Enc_results/DEG_Overlaps.txt', sep = "\t", row.names = F, quote = F)
#write.table(pTable, 'Enc_results/DEG_Overlaps_logp.txt', sep = "\t", row.names = T, quote = F)
#write.table(countTable, 'Enc_results/DEG_Overlaps_counts.txt', sep = "\t", row.names = T, quote = F)

########### DEG overlaps including DEGs from nasopharyngeal swabs - Figure S1A ##########
degs = list(calu3_cov2_up, vero_cov2_up, mm_cov2_up, human_cov2_up,
            calu3_cov2_dn, vero_cov2_dn, mm_cov2_dn, human_cov2_dn)
deg_names = c("Calu-3_CoV-2_UP", "Vero_E6_CoV-2_UP", "Mm_SARS-CoV-2_UP", 'Human_CoV-2_UP',
              "Calu-3_CoV-2_DN",  "Vero_E6_CoV-2_DN", "Mm_SARS-CoV-2_DN", 'Human_CoV-2_DN')

results = deg_overlaps(degs, degs, deg_names, deg_names, alt = 'greater')
pTable = results$pTable
countTable = results$CountTable
# truncating negative logP value matrix
pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)]);
pTable[pTable>50 ] = 50

# Figure S1A
plot_enrichments(pTable, countTable, c(), c(), deg_names, deg_names,
                 mar=c(8,11,1,1), xangle = 45, cex_text = 1.1, cex_lab = 1.1)

results$enrichments = results$enrichments[as.numeric(as.vector(results$enrichments$Overlap_Count)) > 0,]
results$enrichments = results$enrichments[!results$enrichments$DEG1 == results$enrichments$DEG2,]

#write.table(results$enrichments, 'Enc_results/DEG_Overlaps.txt', sep = "\t", row.names = F, quote = F)
#write.table(pTable, 'Enc_results/DEG_Overlaps_logp.txt', sep = "\t", row.names = T, quote = F)
#write.table(countTable, 'Enc_results/DEG_Overlaps_counts.txt', sep = "\t", row.names = T, quote = F)

########## Lung scRNA-seq marker enrichment among individual DEGs - Figures S2A and S2B ########
degs = list(calu3_cov2_up, vero_cov2_up, mm_cov2_up, 
            calu3_cov2_dn, vero_cov2_dn, mm_cov2_dn)
deg_names = c("Calu-3_CoV-2_UP", "Vero_E6_CoV-2_UP", "Mm_SARS-CoV-2_UP", 
              "Calu-3_CoV-2_DN",  "Vero_E6_CoV-2_DN", "Mm_SARS-CoV-2_DN")

# reading lung scRNA seq markers from 3 different studies
all_markers = read_markers()
lung_markers_krasnow = all_markers[[1]]
lung_markers_kaminski = all_markers[[2]] 
lung_markers_banovich = all_markers[[3]]

enc_results_krasnow = data.frame()
enc_results_kaminski = data.frame()
enc_results_banovich = data.frame()

for(i in 1:length(degs)){
  enc_results_krasnow = rbind(enc_results_krasnow, 
                              module_enrichment(lung_markers_krasnow, degs[[i]], deg_names[i])$enrichments)
  enc_results_kaminski = rbind(enc_results_kaminski,
                               module_enrichment(lung_markers_kaminski, degs[[i]], deg_names[i])$enrichments)
  enc_results_banovich = rbind(enc_results_banovich,
                               module_enrichment(lung_markers_banovich, degs[[i]], deg_names[i])$enrichments)
}
enc_results_krasnow = enc_results_krasnow[as.numeric(as.vector(enc_results_krasnow$Overlap_Count)) > 0,]
enc_results_krasnow$`Cell Type` = paste(as.vector(enc_results_krasnow$`Cell Type`),"(Krasnow)")

enc_results_kaminski = enc_results_kaminski[as.numeric(as.vector(enc_results_kaminski$Overlap_Count)) > 0,]
enc_results_kaminski$`Cell Type` = paste(as.vector(enc_results_kaminski$`Cell Type`),"(Kaminski)")

enc_results_banovich = enc_results_banovich[as.numeric(as.vector(enc_results_banovich$Overlap_Count)) > 0,]
enc_results_banovich$`Cell Type` = paste(as.vector(enc_results_banovich$`Cell Type`),"(Banovich)")

enc_results_combined = rbind(enc_results_krasnow, enc_results_kaminski, enc_results_banovich)
View(enc_results_combined)
#write.table(enc_results_combined, 'Enc_results/DEG_Marker_enrichments.txt', sep = "\t", row.names = F, quote = F)

### generate logPvalue matrices for further visualizations
logp_krasnow = generate_logps_degs(lung_markers_krasnow, degs, deg_names)
logp_krasnow = cbind(logp_krasnow, rep("Krasnow", dim(logp_krasnow)[1]))

logp_kaminski = generate_logps_degs(lung_markers_kaminski, degs, deg_names)
logp_kaminski = cbind(logp_kaminski, rep("Kaminski", dim(logp_kaminski)[1]))

logp_banovich = generate_logps_degs(lung_markers_banovich, degs, deg_names)
logp_banovich = cbind(logp_banovich, rep("Banovich", dim(logp_banovich)[1]))

# Supplemental Figures S2A and S2B
logp_lung = rbind(logp_krasnow, logp_kaminski, logp_banovich)
#write.table(logp_lung, "Enc_results/DEG_Marker_logp.txt", sep = "\t", quote = F, row.names = T)

############### Marker enrichment among conserved DEG signature ##############
# reading conserved signatures -> differentially expressed in at least 2 out 3 studies
# Supplemental Table S2
covid19_cons_up = as.vector(read.table("input_data/SARS-CoV2_DEGS/Conserved_UP.txt", header = F)$V1)
covid19_cons_dn = as.vector(read.table("input_data/SARS-CoV2_DEGS/Conserved_DN.txt", header = F)$V1)

degs1 = list()
degs1[[1]] = covid19_cons_up
degs1[[2]] = covid19_cons_dn
deg_names1 = c("SARS-CoV-2-Up(Conserved)", "SARS-CoV-2-Down(Conserved)")

enc_results_krasnow = data.frame()
enc_results_kaminski = data.frame()
enc_results_banovich = data.frame()

for(i in 1:length(degs1)){
  enc_results_krasnow = rbind(enc_results_krasnow, 
                              module_enrichment(lung_markers_krasnow, degs1[[i]], deg_names1[i])$enrichments)
  enc_results_kaminski = rbind(enc_results_kaminski,
                               module_enrichment(lung_markers_kaminski, degs1[[i]], deg_names1[i])$enrichments)
  enc_results_banovich = rbind(enc_results_banovich,
                               module_enrichment(lung_markers_banovich, degs1[[i]], deg_names1[i])$enrichments)
}

enc_results_krasnow$`Cell Type` = paste(as.vector(enc_results_krasnow$`Cell Type`),"(Krasnow)")
enc_results_kaminski$`Cell Type` = paste(as.vector(enc_results_kaminski$`Cell Type`),"(Kaminski)")
enc_results_banovich$`Cell Type` = paste(as.vector(enc_results_banovich$`Cell Type`),"(Banovich)")

# Supplemental Table S6 (also used for figures S2C - S2G)
enc_results_combined = rbind(enc_results_krasnow, enc_results_kaminski, enc_results_banovich)
enc_results_combined = enc_results_combined[as.numeric(as.vector(enc_results_combined$Overlap_Count))>0,]
#write.table(enc_results_combined, 'Enc_results/cons_DEG_Marker_enrichments.txt', sep = "\t", row.names = F, quote = F)

# generating negative logP-value matrices for further visualizations
logp_krasnow = generate_logps_degs(lung_markers_krasnow, degs1, deg_names1)
logp_krasnow = cbind(logp_krasnow, rep("Krasnow", dim(logp_krasnow)[1]))

logp_kaminski = generate_logps_degs(lung_markers_kaminski, degs1, deg_names1)
logp_kaminski = cbind(logp_kaminski, rep("Kaminski(Lung)", dim(logp_kaminski)[1]))

logp_banovich = generate_logps_degs(lung_markers_banovich, degs1, deg_names1)
logp_banovich = cbind(logp_banovich, rep("Banovich(Lung)", dim(logp_banovich)[1]))

logp_lung = rbind(logp_krasnow, logp_kaminski, logp_banovich)
#write.table(logp_lung, "Enc_results/cons_DEG_Marker_logp.txt", sep = "\t", quote = F, row.names = T)

############# COVID19 Core conserved DEGS #######
# reading conserved signatures -> differentially expressed in all 3 studies
# Supplemental Table S2
core_cons_up = as.vector(read.table("input_data/SARS-CoV2_DEGS/Core_conserved_UP.txt", header = F)$V1)
core_cons_dn = as.vector(read.table("input_data/SARS-CoV2_DEGS/Core_conserved_DN.txt", header = F)$V1)

degs1 = list()
degs1[[1]] = core_cons_up
degs1[[2]] = core_cons_dn
deg_names1 = c("SARS-CoV-2-Up(Core)", "SARS-CoV-2-Down(Core)")

enc_results_krasnow = data.frame()
enc_results_kaminski = data.frame()
enc_results_banovich = data.frame()

for(i in 1:length(degs1)){
  enc_results_krasnow = rbind(enc_results_krasnow, 
                              module_enrichment(lung_markers_krasnow, degs1[[i]], deg_names1[i])$enrichments)
  enc_results_kaminski = rbind(enc_results_kaminski,
                               module_enrichment(lung_markers_kaminski, degs1[[i]], deg_names1[i])$enrichments)
  enc_results_banovich = rbind(enc_results_banovich,
                               module_enrichment(lung_markers_banovich, degs1[[i]], deg_names1[i])$enrichments)
}

enc_results_krasnow$`Cell Type` = paste(as.vector(enc_results_krasnow$`Cell Type`),"(Krasnow)")
enc_results_kaminski$`Cell Type` = paste(as.vector(enc_results_kaminski$`Cell Type`),"(Kaminski)")
enc_results_banovich$`Cell Type` = paste(as.vector(enc_results_banovich$`Cell Type`),"(Banovich)")

# Supplemental Table S3
enc_results_combined = rbind(enc_results_krasnow, enc_results_kaminski,enc_results_banovich)
enc_results_combined = enc_results_combined[as.numeric(as.vector(enc_results_combined$Overlap_Count))>0,]
#write.table(enc_results_combined, 'Enc_results/Core_cons_DEG_Marker_enrichments.txt', 
            #sep = "\t", row.names = F, quote = F)

########## Marker enrichment in gene modules from SARS-CoV-2 DEG interactome - Figure 4 ########

# reading gene clusters from conserver DEG-PPI interactome in SARS-CoV-2.
# Steps followed to obtain the gene modules:
  # 1. conserved DEGs are combined with SARS-CoV2-Human PPI genes (336 genes)
  # 2. The combined gene set was uploaded to STRING (https://string-db.org/) to identify SARS-CoV-2 interactome
  # 3. The combined interactome is analysed using the Cytoscape (https://cytoscape.org/) tool.
  # 4. clusterMaker2 (http://www.rbvi.ucsf.edu/cytoscape/clusterMaker2/) plugin was used to implement MCL clustering

# Supplemental Table S9 (used for Figure 3)
node_stats = read.csv('input_data/SARS-CoV-2-Cons_MCL_Clusters.txt', sep = "\t")
# removing unclustered geness
node_stats = node_stats[!node_stats$MCL.cluster == '',]
# candidate modules => at least 5 genes
selected_clusters = names(table(node_stats$MCL.cluster)[table(node_stats$MCL.cluster)>4])

cluster_genes = list()
for(i in selected_clusters){
  genes = as.vector(node_stats[node_stats$MCL.cluster==i,]$name)
  cluster_genes[[i]] = genes
}

enc_results_krasnow = data.frame()
enc_results_kaminski = data.frame()
enc_results_banovich = data.frame()

for(i in 1:length(cluster_genes)){
  enc_results_krasnow = rbind(enc_results_krasnow, 
                              module_enrichment(lung_markers_krasnow, cluster_genes[[i]], selected_clusters[i])$enrichments)
  enc_results_kaminski = rbind(enc_results_kaminski,
                               module_enrichment(lung_markers_kaminski, cluster_genes[[i]], selected_clusters[i])$enrichments)
  enc_results_banovich = rbind(enc_results_banovich,
                               module_enrichment(lung_markers_banovich, cluster_genes[[i]], selected_clusters[i])$enrichments)
}

enc_results_krasnow$`Cell Type` = paste(as.vector(enc_results_krasnow$`Cell Type`),"(Krasnow)")
enc_results_kaminski$`Cell Type` = paste(as.vector(enc_results_kaminski$`Cell Type`),"(Kaminski)")
enc_results_banovich$`Cell Type` = paste(as.vector(enc_results_banovich$`Cell Type`),"(Banovich)")

# Supplemental Table S10
enc_results_combined = rbind(enc_results_krasnow, enc_results_kaminski, enc_results_banovich)
enc_results_combined = enc_results_combined[as.numeric(as.vector(enc_results_combined$Overlap_Count))>0,]
#write.table(enc_results_combined, 'Enc_results/SARS-COV-2_Cons_PPI_module_Marker_enrichments.txt', 
            #sep = "\t", row.names = F, quote = F)
########## Phenotype enrichments (PheGenI) among conserved DEGs - Figure S3 ###########
phegen_info = read.csv('input_data/other data/PheGenI_Association_latest.txt', sep = "\t", header = T)
phegen_info = phegen_info[!phegen_info$Context=='intergenic',]
phegen_info = phegen_info[phegen_info$P.Value<0.00001,]

bg_cnt = length(union(unique(phegen_info$Gene), unique(phegen_info$Gene.2)))
# reading conserved signatures -> differentially expressed in at least 2 out 3 studies
# Supplemental Table S2
covid19_cons_up = as.vector(read.table("input_data/SARS-CoV2_DEGS/Conserved_UP.txt", header = F)$V1)
covid19_cons_dn = as.vector(read.table("input_data/SARS-CoV2_DEGS/Conserved_DN.txt", header = F)$V1)
covid19_cons_degs = union(covid19_cons_up, covid19_cons_dn)

cons_trait_enrichments = rbind(trait_enrichment(phegen_info, covid19_cons_up, "COVID19_Up(Conserved)", bg_cnt = bg_cnt)$enrichments,
                               trait_enrichment(phegen_info, covid19_cons_dn, "COVID19_DN(Conserved)", bg_cnt = bg_cnt)$enrichments,
                               trait_enrichment(phegen_info, covid19_cons_degs, "COVID19_DEG(Both)", bg_cnt = bg_cnt)$enrichments)
cons_trait_enrichments = cons_trait_enrichments[as.numeric(as.vector(cons_trait_enrichments$Overlap_Count))>0,]
View(cons_trait_enrichments)

# Supplemental Table S7
#write.table(cons_trait_enrichments, "Enc_results/COVID19_Cons_PheGenI_overlap.txt", sep = "\t")

########## Phenotype enrichments (PheGenI) among SARS-CoV-2 DEG interactome #######

# reading conserved signatures -> differentially expressed in at least 2 out 3 studies
# Supplemental Table S2
covid19_cons_up = as.vector(read.table("input_data/SARS-CoV2_DEGS/Conserved_UP.txt", header = F)$V1)
covid19_cons_dn = as.vector(read.table("input_data/SARS-CoV2_DEGS/Conserved_DN.txt", header = F)$V1)
covid19_cons_degs = union(covid19_cons_up, covid19_cons_dn)

cons_trait_enrichments = rbind(trait_enrichment(phegen_info, covid19_cons_up, "COVID19_Up(Conserved)", bg_cnt = bg_cnt)$enrichments,
                               trait_enrichment(phegen_info, covid19_cons_dn, "COVID19_DN(Conserved)", bg_cnt = bg_cnt)$enrichments,
                               trait_enrichment(phegen_info, covid19_cons_degs, "COVID19_DEG(Both)", bg_cnt = bg_cnt)$enrichments)
cons_trait_enrichments = cons_trait_enrichments[as.numeric(as.vector(cons_trait_enrichments$Overlap_Count))>0,]
# Supplemental Table S7
#write.table(cons_trait_enrichments, "Enc_results/COVID19_Cons_PheGenI_overlap.txt", sep = "\t")
# Supplemental Table S9 (used for Figure 3)
node_stats = read.csv('input_data/SARS-CoV-2-Cons_MCL_Clusters.txt', sep = "\t")
# removing unclustered geness
node_stats = node_stats[!node_stats$MCL.cluster == '',]
# candidate modules => at least 5 genes
selected_clusters = names(table(node_stats$MCL.cluster)[table(node_stats$MCL.cluster)>4])

cluster_trait_enrichments = data.frame()
for(cid in selected_clusters){
  print(cid)
  genes = as.vector(node_stats[node_stats$MCL.cluster==cid,]$name)
  enrichments = trait_enrichment(phegen_info, genes, cid, bg_cnt = bg_cnt)$enrichments
  enrichments = enrichments[as.numeric(as.vector(enrichments$Overlap_Count))>0,]
  cluster_trait_enrichments = rbind(cluster_trait_enrichments, enrichments)
}
cluster_trait_enrichments = cluster_trait_enrichments[as.numeric(as.vector(cluster_trait_enrichments$Overlap_Count))>0,]

# Supplemental Table S12
# write.table(cluster_trait_enrichments, "Enc_results/SARS_COV_2_PPI_module_PheGenI_overlap.txt", sep = "\t")