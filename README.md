# Predicting natural killer cell infiltration in metastatic melanoma
This repository contains code for the manuscript:
* Cursons *et al*. (*in press*). *Cancer Immunology Research*. A gene signature predicting natural killer cell infiltration and improved survival in metastatic melanoma patients. **DOI**: not-yet-known.

For further information please contact:
* Dr Joe Cursons: cursons.j (at) wehi.edu.au
* Dr Momeneh (Sepideh) Foroutan: momeneh.foroutan (at) unimelb.edu.au
* Prof. Nicholas D Huntington: nicholas.huntington (at) monash.edu.au
* Dr Melissa J Davis: davis.m (at) wehi.edu.au

## Overview
Code in this repository includes contributions from:
* Joe Cursons (cursons.j (at) wehi.edu.au)
   - analysis.py :: a python script to perform the majority of analyses and produces figures for the associated manuscript. A more detailed overview of the classes and functions within this package is given below

* Momeneh (Sepideh) Foroutan (momeneh.foroutan (at) unimelb.edu.au) 
   - data :: a folder containing extra data used in this analysis
      - bldCells_DGEList.RData :: gene expression/abundance data from blood cells used for the differential expression analysis as part of NK signature curation (original data are from Linsley et al. (2014). *PLoS One*.; differential expression analysis script is noted below)
	  - Cursons_Guimaraes_NKsignature_CIR_2019.csv :: genes from the NK signature created for this report using Supplementary Table S1 in the paper.
	  - Foroutan2016_TGFb_EMT_signature_upDown.txt :: genes for the TGF-B EMT signature from Foroutan et al. (2016). *Mol. Canc. Res.*
	  - Homo_sapiens.gene_info :: gene mapping information used
	  - Thiery_EMTsignature_both_tumour_cellLine_EntrezIDs.txt :: Epithelial and mesenchymal gene sets from Tan et al. [JP Thierry] (2014). *EMBO Mol. Med.*
   - reports :: Rmarkdown reports
      - MANIFEST.txt :: The NIH/NCI Genomic Data Commons manifest file for TCGA SKCM samples used in this analysis (this is automatically downloaded when using TCGAbiolinks package to download the data) 
	  - Human_genes__GRCh38_p12_.rda :: RData for human gene annotation (this is automatically downloaded when using TCGAbiolinks package to download the data)
	  - NK_scoring_survival.html/.Rmd :: Reports showing an alternative implementation of the NK scoring survival analysis. 
	    - **NB:** some of the parameters used for thresholding scores as well as the format of the TCGA data vary slightly in this analysis and it is not an exact reproduction of the python script used for the associated manuscript, although sample scores etc are consistent
	  - NK_singleCell_SadeFeldman.html/.Rmd :: the single cell RNA-seq analysis of melanoma data from Sade-Feldman et al. (2018). *Cell*. - included as part of Fig. 3 within the manuscript
	  - NK_singleCell_Tirosh.html/.Rmd :: the single cell RNA-seq analysis of melanoma data from Tirosh et al. (2015). *Science* - included as part of Fig. 2 & Fig. 3 within the manuscript
   - script :: R scripts
      - boxplot_vertical_ggplot.R :: A custom fuction used to generate boxplot of Tirosh data included as Fig. 3A in NK_singleCell_Tirosh.html/.Rmd report
	  - DE_NK_Genes_BloodBulk.R :: differential expression analysis to identify genes expressed at higher levels in NK cells relative to other blood populations (using data from Linsley et al noted above).
	  - Survival_analysis.R :: A custom function for R-based reimplementation of survival analysis for genes investigated in the NK_scoring_survival.html/.Rmd report
	    - **NB:** this is not the version used in the manuscript, it is an alternative implementation showing consistent results
   
   
## Data
The authors would like to acknowledge the authors of previous works who have made their data freely available for further analyses, as well as patients who generously donated tissue for these data. Data analysed here include:
* The Cancer Genome Atlas skin cutaneous melanoma (TCGA-SKCM) data

   - The TCGA SKCM data represent a large cohort of matched RNA-seq, patient outcome, as well as genomic sequencing data. RNA-seq and patient outcome data were used here.
   - Please refer to the original TCGA SKCM manuscript for further details: Genomic Classification of Cutaneous Melanoma. *Cell* (2015). 161(7): pp. 1681-1696. **DOI**: 10.1016/j.cell.2015.05.044
   - Data can be downloaded directly from the NIH/NCI genomic data commons (GDC): https://portal.gdc.cancer.gov/projects/TCGA-SKCM
   
* The LM-MEL cell line panel:

   - Behren A, *et al* (2013). The Ludwig Institute for Cancer Research Melbourne Melanoma Cell Line Panel. *Pigment Cell & Melanoma Research*. 26: 597-600. **DOI**: 10.1111/pcmr.12097
   - Data were downloaded directly from ArrayExpress with the accession number: E-MTAB-1496

* Tirosh *et al*. single cell RNA-seq data from dissociated melanoma samples

   - Tirosh I, *et al*. (2016). Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq. *Science*. 352(6282): 189-96. **DOI**: 10.1126/science.aad0501 
   - Data were downloaded directly from the NCBI gene expression omnibus with accession number: GSE72056

* Sade-Feldman *et al*. single cell RNA-seq data from dissociated melanoma samples 

   - Sade-Feldman M, *et al*. (2018). Defining T Cell States Associated with Response to Checkpoint Immunotherapy in Melanoma. *Cell*. 175(4): 998-1013.e20. **DOI**: 10.1016/j.cell.2018.10.038
   - Data were downloaded directly from the NCBI gene expression omnibus with accession number: GSE120575
     - **NB**: only pre-treatment samples were used for this analysis

* Linsley *et al*. RNA-seq data from sorted immune cell populations: 

   - Linsley PS, *et al*. (2014). Copy number loss of the interferon gene cluster in melanomas is linked to reduced T cell infiltrate and poor patient prognosis. *PLoS One*. 9(10):e109760. **DOI**: 10.1371/journal.pone.0109760
   - Data were downloaded directly from the NCBI gene expression omnibus with accession number: GSE60424
   
* Novershtern *et al*. microarray data from sorted immune cell population:

   - Novershtern N, *et al*. (2011). Densely interconnected transcriptional circuits control cell states in human hematopoiesis. *Cell*. 144(2): 296-309. **DOI**: 10.1016/j.cell.2011.01.004
   - Data were downloaded directly from the NCBI gene expression omnibus with accession number: GSE24759

## Python analysis script

### Overview
The majority of analyses included within the manuscript were performed using a python script (analysis.py). This script also contains the code to produce the published figures.

### Dependencies
The authors would like to convey their gratitude to the developers of python libraries used in this analysis. In particular:
* [lifelines](https://lifelines.readthedocs.io/en/latest/)
   - [GitHub](https://github.com/CamDavidsonPilon/lifelines/)
* [matplotlib](https://matplotlib.org/)
* [networkx](https://networkx.github.io/)
* [numpy](https://www.numpy.org/)
* [pandas](https://pandas.pydata.org/)
* [scipy](https://www.scipy.org/)

### Script classes
* PreProc: a series of functions which perform pre-processing on different data sets to make sample accession easier for corresponding subsets of samples etc
   - split_tcga_met_vs_pri()
   - tcga_skcm_rna_data()
   - tcga_skcm_data()
   - tcga_skcm_classifications()
   - tcga_skcm_met_sites()
   - lm_mel_data()
   - gse60424_data()
   - gse24759_data()
   - gse24759_subsets()
   - tcga_histology()
   - density_scatters()
   - refine_NK_signature()
* Analyse:
   - 
* Plot: