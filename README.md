# COVID19_secondary_analysis
## Secondary analysis of transcriptomes of SARS-CoV-2 infection models to characterize COVID-19
<br/>
<p align="center"><img src="graphical_abstract.png" style="vertical-align:middle" width="600" height="600"></p>

This repository contains the various input files, generated outputs and scripts associated with our paper titled above. These scripts can be used to generate or reproduce the data published in our main text and other supplemental items. We also share the R objects and cytoscape (https://cytoscape.org/) session files associated with all the figures published in our work. The folder structure of this repository is described below:

* <u><b>input_data/</b></u>: This directory contains the various files used as inputs to our study

  * <u><b>Count data/</b></u> - Raw counts from the three SARS-CoV-2 studies (two <i>in vitro</i> models and one <i>in vivo</i> model) used in our research. In case of the two animal models, the corresponding human orthologs used are also available
  
  * <u><b>Lung Markers/</b></u> - Lung scRNA-seq markers from three different human lung studies utilized in our work.
  
  * <u><b>SARS-CoV-2 DEGs/</b></u> - Individual differentially expressed gene (DEG) lists identified from the three input studies. along with the consensus transcriptomic signature.Also included are the DEGs from nasopharyngeal swabs from human COVID-19 patients (GSE152075) and SARS-CoV-2 human interactants.
  
  * <u><b>other data/</b></u> - this folder contains gene-phenotype/trait associations (compressed files) from both GWAS Catalog (https://www.ebi.ac.uk/gwas/) and PheGenI (https://www.ncbi.nlm.nih.gov/gap/phegeni) used in our study.
  
* <u><b>Scripts/</b></u>: This directory includes the script files used in our analysis.

   * <u>COVID_enrichments.R</u> - For identifying enriched cell types and traits in the SARS-CoV-2 consensus transcriptome and candidate protein modules from the integrated SARS-CoV-2 interactome. Also useful for producing the Supplemental tables from our work.
   
   * <u>Utils.R</u> - Contains helper functions used in our enrichment analysis script.
   
   * <u>COVID_benchmarking.R</u> - Randomized trials conducted to test the robustness of individual DEGs and the consensus transcriptome from the three input disease models used in our framework. Also included are the experiments used to validate the level of connectivity observed among the consensus signature along with their interactions with the SARS-CoV-2 virus-host interactants.
   
   * <u>MCL_Clustering.R</u> - Optional script provided to identify protein modules by implementing MCL clustering algorithm. This scrippt takes the consensus signature and the SARS-CoV-2 virus-host interactants and constructs an interactome network using filtered protein-protein interactions from STRING (https://string-db.org/). File containing the filtered interactions (total_score >= 900; experimental_score >= 700) can be found at <i>input_data/other data/</i>
   
   * <u>phenotype_enrichments.ipynb</u> - iPython notebook that we used to parse child traits/terms from experimental factor ontology (EFO) hierarchy which are then used to compute enrichment of phenotype associations (including the child traits) from GWAS Catalog.

* <u><b>RData/</b></u> - R objects to reproduce the results from our benchmarking experiments.

* <u><b>Figures_data/</b></u> - Contains the cytoscape session files to generate the network figures (both main and Supplemental) presented in our work. These session objects also contain the input networks that were used to generate the visualizations


  