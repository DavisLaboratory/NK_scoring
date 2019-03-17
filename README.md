# Predicting natural killer cell infiltration in metastatic melanoma
This repository contains code for the manuscript:
* Cursons *et al*. (*in press*). *Cancer Immunology Research*. A gene signature predicting natural killer cell infiltration and improved survival in metastatic melanoma patients. DOI: not-yet-known.

For further information please contact:
* Dr Joe Cursons: cursons.j (at) wehi.edu.au
* Prof. Nicholas D Huntington: nicholas.huntington (at) monash.edu.au
* Dr Melissa J Davis: davis.m (at) wehi.edu.au

## Overview
Code in this repository includes:
* analysis.py :: a python script to perform the majority of analyses and produces figures for the manuscript listed above

   This script was written by Joe Cursons (cursons.j (at) wehi.edu.au)
   
* singlecell.Rmd :: scripts for analysing single cell RNA-seq data
* R implementation/NK_analysis_FPKM :: an independent R implementation of the survival analysis which generates similar results to the python script above

   These scripts were written by Momeneh (Sepideh) Foroutan (momeneh.foroutan (at) unimelb.edu.au) 
   
## Data
The authors would like to acknowledge the authors of previous works who have made their data freely available for further analyses. Data analysed here include:
* The Cancer Genome Atlas skin cutaneous melanoma (TCGA-SKCM) data

   The TCGA SKCM data represent a large cohort of matched RNA-seq, patient outcome, as well as genomic sequencing data. RNA-seq and patient outcome data were used here.
   Please refer to the original TCGA SKCM manuscript for further details: Genomic Classification of Cutaneous Melanoma. *Cell* (2015). 161(7): pp. 1681-1696. DOI: 10.1016/j.cell.2015.05.044
   Data can be downloaded directly from the NIH/NCI genomic data commons (GDC): https://portal.gdc.cancer.gov/projects/TCGA-SKCM
   The authors would like to thank the patients who generously donated tissue for these data.
   
* The LM-MEL cell line panel

* Tirosh *et al*. single cell RNA-seq data from dissociated melanoma samples

   Tirosh I, *et al*. (2016). Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq. *Science*. 352(6282): 189-96. DOI: 10.1126/science.aad0501 
   Data were downloaded directly from the NCBI gene expression omnibus with accession number: GSE72056

* Sade-Feldman *et al*. single cell RNA-seq data from dissociated melanoma samples 

   Sade-Feldman M, *et al*. (2018). Defining T Cell States Associated with Response to Checkpoint Immunotherapy in Melanoma. *Cell*. 175(4): 998-1013.e20. DOI: 10.1016/j.cell.2018.10.038
   Data were downloaded directly from the NCBI gene expression omnibus with accession number: GSE120575
   **NB**: only pre-treatment samples were used for this analysis

* Linsley *et al*. RNA-seq data from sorted immune cell populations: 

   Linsley PS, *et al*. (2014). Copy number loss of the interferon gene cluster in melanomas is linked to reduced T cell infiltrate and poor patient prognosis. *PLoS One*. 9(10):e109760. DOI: 10.1371/journal.pone.0109760
   Data were downloaded directly from the NCBI gene expression omnibus with accession number: GSE60424
   
* Novershtern *et al*. microarray data from sorted immune cell population:

   Novershtern N, *et al*. (2011). Densely interconnected transcriptional circuits control cell states in human hematopoiesis. *Cell*. 144(2): 296-309. DOI: 10.1016/j.cell.2011.01.004