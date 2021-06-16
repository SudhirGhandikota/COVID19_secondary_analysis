##### Installing Package RUVSeq ############
options(warn=-1)
if (!require("devtools")) install.packages("devtools",repos="http://cran.us.r-project.org")
library(withr)
library(httr)
library(curl)
library(devtools)
if (!require("RUVSeq")) BiocManager::install("RUVSeq")
library(RUVSeq)

####### Reading Command Line Arguments ########
args = commandArgs(trailingOnly=TRUE)
infile = args[1]
num_controls = as.numeric(args[2])
num_cases <- as.numeric(args[3])
counts_threshold <- as.numeric(args[4])
num_samples <- as.numeric(args[5])

print(paste("Number of controls:", num_controls, "Number of cases:", num_cases,
            "Count threshold:", counts_threshold, "Sample threshold:", num_samples))

# returns differentially expressed genes
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
  result <- topTags(lrt, n=dim(y)[1]+1, adjust.method="BH", sort.by="logFC")
  return(result$table)
}

# countData = as.matrix(read.table('../input_data/Count Data/GSE147507/Genes_Counts_Matrix_Human_Series7.txt', 
#                                   #sep="\t", header=TRUE, row.names=1, check.names=F))
countData = as.matrix(read.table(infile, sep="\t", header=TRUE, row.names=1, check.names=F))
result_DE = get_degs(countData, num_controls, num_cases, counts_threshold, num_samples)

# saving DE results
outdir = dirname(infile)
outfile = paste(outdir, "/", "DE_results.txt", sep="")
print(paste("Outpath:", outfile))
write.table(result_DE, outfile, sep="\t", row.names = T, quote = F, col.names = T)