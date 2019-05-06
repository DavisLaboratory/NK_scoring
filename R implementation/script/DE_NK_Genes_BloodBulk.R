
##======================= Author =======================
## ------------ Momeneh (Sepideh) Foroutan -------------
## ------------ Last updated: 03 May 2019 --------------
##======================================================
## In this script, we performs DE analysis on the bulk RNA-seq data from a blood data set.
## The SRA files for the blood data have been downloaded in Sep 2016 from GSE60424, and analysed using Rsubread, featureCounts, and edgeR package. Here we read the DGEList for this data which includes the count output from featureCounts, as well as sample annotations.
## After performing DE analysis, we subset the test statistics for the genes in the combined NK signature.
## The results of this script contribute to the Supplementary Table S1 in Cursons et al paper (A gene signature predicting natural killer cell infiltration and improved survival in melanoma patients).

#########################################################

##---- Load required libraries and set the paths
library(edgeR)
library(Homo.sapiens) 
dataPath <- "./data/"

##---- Read combined NK signature genes
nkSig <- read.csv(paste0(dataPath, "Cursons_Guimaraes_NKsignature_CIR_2019.csv"),
  stringsAsFactors = F)

nk_signature <- nkSig$HGNC.Symbol

## Make sure that we have updated gene symbols
nkIDs <- alias2SymbolUsingNCBI(nk_signature, 
  paste0(dataPath, "Homo_sapiens.gene_info"))

nkIDs <- nkIDs[complete.cases(nkIDs$GeneID), ]

##----- Read DGEList object for the expression
load(paste0(dataPath, "bldCells_DGEList.RData"))
dge <- bldCells_DGEList

##---- Divide samples into Nk and non_NK samples
dge$samples$celltype_s <- as.character(dge$samples$celltype_s)
dge$samples$groups[dge$samples$celltype_s == "NK"] <- "NK"
dge$samples$groups[ ! dge$samples$celltype_s == "NK"] <- "non-NK"

##---- remove the whole blood samples
dge <- dge[, ! dge$samples$celltype_s == "Whole Blood", ]

##---- Filter genes
kp <- rowSums(cpm(dge$counts) > 4) >= 4 | row.names(dge$counts) %in% nkIDs$GeneID
summary(kp)
dgeFilt <- dge[kp, ]

##---- Perform TMM normalisation
dgeFiltNorm <- calcNormFactors(dgeFilt)


##----- Perform DE analysis (we use voom-limma + TREAT)
design <- model.matrix(~ 0 + dge$samples$groups)
colnames(design) <- c("NK", "nonNK")

contr.matrix <- makeContrasts(
  NK_nonNK = NK - nonNK,
  levels = colnames(design))

v <- voom(dgeFiltNorm, design, plot=TRUE)

## Two genes from the combined NK signature gene list do not present in the expression data:
nkIDs$GeneID[ ! nkIDs$GeneID %in% row.names(dgeFiltNorm$counts)]
## "3803"  "28639"

vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)

## First we use eBayes but it gives many DEGs
efit <- eBayes(vfit)
plotSA(efit)
summary(decideTests(efit))

## Use TREAT to lower the number of DEGs
tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit)
summary(dt)

##----- Map gene IDs to symbols
geneid <- rownames(dgeFiltNorm) 

genes <-
  select(
    Homo.sapiens,
    keys = geneid,
    columns = c("SYMBOL", "TXCHROM"),
    keytype = "ENTREZID"
  ) 
head(genes)

genes <- genes[!duplicated(genes$ENTREZID), ]

DEstats   <- topTreat(tfit, coef = 1, n = Inf)
DEstatsID <- merge(genes, DEstats, by.x = "ENTREZID", by.y = "GeneID")

## Get the stats for the all NK genes in the combined signature.
sigCombIDs_stat <- DEstatsID[DEstatsID$ENTREZID %in% nkIDs$GeneID, ]
sigCombIDs_stat <- sigCombIDs_stat[, -c(3, 4, 5, 6, 7)]

## Export the results
write.table(sigCombIDs_stat, paste0(outPath, "DEstat_sigCombined.txt"), row.names = F, sep = "\t")





