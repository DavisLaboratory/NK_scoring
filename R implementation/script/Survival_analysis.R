
##======================= Author =======================
## ------------ Momeneh (Sepideh) Foroutan -------------
## ------------ Last updated: 02 May 2019 --------------
##======================================================

# This function depends on the below libraries 
# library(survival)
# library(ggfortify)
# library(RColorBrewer)

#=================== Survival analysis =================
# This function looks at the associations between different variables and survival outcome. These variables can be one of the below options that stratifies samples for survival analysis:
  
#  **expr** : expression of a gene, will be split based on 33%-tile and 66%-tile (e.g. low, medium, high)\n
# **score** : score of a single signatre, will be split based on 33%-tile and 66%-tile (low, medium, high)\n
# **covariate** : A continouse covariate (e.g. age), will be split based on 33%-tile and 66%-tile (low, medium, high)\n
# **score_expr** : stratifies samples based on scores from a signature (high and low) and expression of a gene (high and low)\n
# **covariate_expr** : startifies samples according to covariate (age; high and low) and expression of a gene (high and low)\n
# **score_covariate** : stratifies samples according to scores from a single signature (high and low) and covariate (age; high and low)\n
# **expr_expr** : stratifies samples according to expression of two genes (high gene1/high gene2, high gene1/low gene2, etc)\n
# **score_score** : stratifies samples according to scores obtained from two signatures (high score1/high score2, high score1/low score2, etc)

##---------------------- INPUT ARGUMENTS ------------------------

# 1. **data**: a `SummarizedExperiment` object; for example, here we use the TCGA SKCL data that we downloaded using TCGAbiolink package (see above). This is a specific data object in R that stores expression data as well as several meta-data. Therefore, this function can not take a data frame as input at this stage. The `SummarizedExperiment` object needs to further have an `assay` slot called "logFPKM". There must be an "external_gene_name" column in the rowData that has gene names in the same format as the genes that we provide in the gene arguments of the function. To learn more about SummarizedExperiment object and how to construct it, please see https://bioconductor.org/packages/release/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html.
# 2. **stratify**: A character value of one of the above listed options for stratification (e.g. "expr", "score_expr").
# 3. **scores**: A data frame with maximum of three columns: one column needs to be called "sample" which has the sample names consistent with sample names in expression data (first argument of teh function), and minimum one or maximum two columns of signature scores, which have "score" as part of their column names. An error will be given if the data has more than two score columns. This argument can be `NULL` if you are not inetersted in the relationship between scores and survival outcome.
# 4. **gene**: A character vector containing the names of maximum of two genes. This argument can be `NULL` if you are not inetersted in the relationship between genes and survival outcome.
# 5. **covariate**: Name of the column for a covariate; This is the age factor by default. This argument can be `NULL` if you are not inetersted in the relationship between the covariate and survival outcome.
# 6. **timeCol**: The name of time column to be used in the survival analysis
# 7. **eventCol**: The name of the event column to be used in survival analysis (e.g. vital_status)
# 8. **nGroup**: The number of groups for each stratification. Can be 2 or 3. For example, a value of 2 (default) generates two groups of samples with high and low expression for a desired gene/score/covariate, while 3 would stratify samples into three groups of high, low, and medium.
# 9. **confInt**: Boolean; if TRUE, the confidence intervals of survival curves are plotted.

##------------------------- FUNCTION ---------------------------

survival_analysis <- function(data,
  stratify = "expr_expr",
  scores = NULL,
  gene = c("XCL1", "CCL5"),
  covariate = NULL,
  timeCol = "os.delay.ch1",
  eventCol = "os.event.ch1",
  nGroup = 2, 
  confInt = T) {

  ## Check data class
  if( ! class(data) == "RangedSummarizedExperiment" & 
      ! class(data) == "SummarizedExperiment"){
    stop(paste0("The data argument needs to be a 'SummarizedExperiment' object, which has assay, colData and rowData"))
  }
  
  ## Check stratify argument
  strata <- c("expr", "score", "covariate",
    "score_expr", "covariate_expr", "score_covariate",
    "expr_expr", "score_score")
  
  if(! stratify %in% strata) {
    stop(paste0("The 'stratify' argument must be one of the ", strata))
  }
  
  ## Check covariate
  if (! is.null(covariate)){
    if(! covariate %in% colnames(colData(data))) {
      stop(paste0("The covariate column specified does not exist in the colData(data); check colnames(colData(data)) to select a valid column name for covariate"))
    }
  }
  
  ## Check time and event columns
  if (! timeCol %in% colnames(colData(data))){
    stop(paste0("The timeCol specified does not exist in the colData(data); check colnames(colData(data)) to select a valid column name for survival time"))
  }
  
  if (! eventCol %in% colnames(colData(data))){
    stop(paste0("The eventCol specified does not exist in the colData(data); check colnames(colData(data)) to select a valid column name for survival event"))
  }
  
  ##---------------------------------------- Scores
  if (!is.null(scores)) {
    col_names <- colnames(assay(data))
    newAnnot <- merge(
      colData(data),
      scores,
      by = "sample",
      sort = F,
      all.x = T
    )
    
    row.names(newAnnot) <- col_names
    
    scoreCol <- grep("score", colnames(scores), value = T)
    
    if (length(scoreCol) > 2) {
      stop(paste0("score data must have maximum of 2 score columns"))
    }
    
    ## add high/low/medium score:
    currentscoreCol <- scoreCol[1]
    median_score <- median(newAnnot[, currentscoreCol])
    lowQ_score <- quantile(newAnnot[, currentscoreCol], prob = 0.33)
    upQ_score <- quantile(newAnnot[, currentscoreCol], prob = 0.66)
    
    newAnnot$scores_2status <-
      ifelse(
        newAnnot[, currentscoreCol] >= median_score,
        paste("High ", currentscoreCol),
        paste("Low ", currentscoreCol)
      )

    newAnnot$scores_3status[newAnnot[, currentscoreCol] >= upQ_score] <-
      paste("High ", currentscoreCol)
    newAnnot$scores_3status[newAnnot[, currentscoreCol] <= lowQ_score] <-
      paste("Low ", currentscoreCol)
    newAnnot$scores_3status[newAnnot[, currentscoreCol] < upQ_score &
                              newAnnot[, currentscoreCol] > lowQ_score] <-
      paste("Medium ", currentscoreCol)
    
    
    if (length(scoreCol) == 2) {
      currentscoreCol <- scoreCol[2]
      median_score <- median(newAnnot[, currentscoreCol])
      lowQ_score <-
        quantile(newAnnot[, currentscoreCol], prob = 0.33)
      upQ_score <-
        quantile(newAnnot[, currentscoreCol], prob = 0.66)
      
      newAnnot$scores2_2status <-
        ifelse(
          newAnnot[, currentscoreCol] >= median_score,
          paste("High ", currentscoreCol),
          paste("Low ", currentscoreCol)
        )
      
      newAnnot$scores2_3status[newAnnot[, currentscoreCol] >= upQ_score] <-
        paste("High ", currentscoreCol)
      newAnnot$scores2_3status[newAnnot[, currentscoreCol] <= lowQ_score] <-
        paste("Low ", currentscoreCol)
      newAnnot$scores2_3status[newAnnot[, currentscoreCol] < upQ_score &
                                 newAnnot[, currentscoreCol] > lowQ_score] <-
        paste("Medium ", currentscoreCol)
    }

    ## save this new annotation data as sample annotation for the data
    colData(data) <- newAnnot
  }
  
  ##-------------------------------------- Covariate
  
  if (!is.null(covariate)) {
    newAnnot <- colData(data)
    median_cov <- median(newAnnot[, covariate], na.rm = T)
    lowQ_cov <-
      quantile(newAnnot[, covariate], prob = 0.33, na.rm = T)
    upQ_cov <-
      quantile(newAnnot[, covariate], prob = 0.66, na.rm = T)
    
    newAnnot$cov_2status <-
      ifelse(newAnnot[, covariate] >= median_cov,
             "High covariate",
             "Low covariate")
    
    newAnnot$cov_3status[newAnnot[, covariate] >= upQ_cov] <-
      "High covariate"
    newAnnot$cov_3status[newAnnot[, covariate] <= lowQ_cov] <-
      "Low covariate"
    newAnnot$cov_3status[newAnnot[, covariate] < upQ_cov &
                           newAnnot[, covariate] > lowQ_cov] <-
      "Medium covariate"
    
    ## save this new annotation data as sample annotation for the data
    colData(data) <- newAnnot
  }
  
  currentData <- data

  ##------------------------------------- Gene expression
    
  if (!is.null(gene)) {
    if (sum(rowData(data)$external_gene_name %in% gene) < 1) {
      stop(paste0(gene, " : this gene does not present in the expression data"))
    }
    if (length(gene) > 2) {
      stop(paste0("Please provide maximum of 2 genes at a time"))
    }
    
    currentData <-
      currentData[rowData(currentData)$external_gene_name %in% gene, ]
    
    currentGene <- gene[1]
    
    newAnnot <- colData(currentData)
    
    ## calculate median and 33%-tile and 66%-tile of gene expression
    currentData1 <-
      currentData[rowData(currentData)$external_gene_name == currentGene, ]
    
    median_expr <- median(assay(currentData1, "logFPKM"))
    lowQ_expr <- quantile(assay(currentData1, "logFPKM"), prob = 0.33)
    upQ_expr <- quantile(assay(currentData1, "logFPKM"), prob = 0.66)
    
    newAnnot$expr1_2status <-
      ifelse(
        assay(currentData1, "logFPKM") >= median_expr,
        paste0("High ", currentGene),
        paste0("Low ", currentGene)
      )

    newAnnot$expr1_3status[assay(currentData1, "logFPKM") >= upQ_expr] <-
      paste0("High ", currentGene)
    newAnnot$expr1_3status[assay(currentData1, "logFPKM") <= lowQ_expr] <-
      paste0("Low ", currentGene)
    newAnnot$expr1_3status[assay(currentData1, "logFPKM") < upQ_expr &
                             assay(currentData1, "logFPKM") > lowQ_expr] <-
      paste0("Medium ", currentGene)
    
    ## save this new annotation as sample annotation for the data
    colData(currentData) <- newAnnot
    
    if (length(gene) > 1) {
      currentGene <- gene[2]
      newAnnot <- colData(currentData)
      
      ## calculate median and 33%-tile and 66%-tile of gene expression
      median_expr <-
        median(assay(currentData[rowData(currentData)$external_gene_name == currentGene, ], "logFPKM"))
      lowQ_expr <-
        quantile(assay(currentData[rowData(currentData)$external_gene_name == currentGene, ], "logFPKM"), prob = 0.33)
      upQ_expr <-
        quantile(assay(currentData[rowData(currentData)$external_gene_name == currentGene, ], "logFPKM"), prob = 0.66)
      
      newAnnot$expr2_2status <-
        ifelse(
          assay(currentData[rowData(currentData)$external_gene_name == currentGene, ], "logFPKM") >= median_expr,
          paste0("High ", currentGene),
          paste0("Low ", currentGene)
        )
      
      newAnnot$expr2_3status[assay(currentData[rowData(currentData)$external_gene_name == currentGene, ], "logFPKM") >= upQ_expr] <-
        paste0("High ", currentGene)
      newAnnot$expr2_3status[assay(currentData[rowData(currentData)$external_gene_name == currentGene, ], "logFPKM") <= lowQ_expr] <-
        paste0("Low ", currentGene)
      newAnnot$expr2_3status[assay(currentData[rowData(currentData)$external_gene_name == currentGene, ], "logFPKM") < upQ_expr &
                               assay(currentData[rowData(currentData)$external_gene_name == currentGene, ], "logFPKM") > lowQ_expr] <-
        paste0("Medium ", currentGene)
    }
    
    colData(currentData) <- newAnnot
  }
  
  currentData <-
       currentData[, complete.cases(colData(currentData)[, timeCol])]
  
  ##------------------------------------- Check for stratification type
  if (stratify == "expr") {
    currentStrata <- paste0("expr1_", nGroup, "status")
    mainTitle <- gene[1]
    }
  if(stratify == "score"){
    currentStrata <- paste0("scores_", nGroup, "status")
    mainTitle <- scoreCol[1]
  }
  if(stratify == "covariate"){
    currentStrata <- paste0("cov_", nGroup, "status")
    mainTitle <- covariate
  }
  if (stratify == "score_expr") {
    currentSt_score <- paste0("scores_", nGroup, "status")
    currentSt_expr  <- paste0("expr1_", nGroup, "status")
    colData(currentData)$score_expr <-
      paste0(colData(currentData)[, currentSt_score],
        " / ",
        colData(currentData)[, currentSt_expr])
    currentStrata <- "score_expr"
    mainTitle <- paste(scoreCol[1], " & ", gene)
  }
  if (stratify == "covariate_expr") {
    if(is.null(covariate) | is.null(gene)){
      stop("Make sure both covriate and gene are provided")
    }
    currentSt_cov  <- paste0("cov_", nGroup, "status")
    currentSt_expr <- paste0("expr1_", nGroup, "status")
    ## Remove samples with NA annotation as covariate:
    currentData <- currentData[ , complete.cases(colData(currentData)[, currentSt_cov])]
    colData(currentData)$cov_expr <-
      paste0(colData(currentData)[, currentSt_cov],
        " / ",
        colData(currentData)[, currentSt_expr])
    currentStrata <- "cov_expr"
    mainTitle <- paste(covariate, " & ", gene[1])
  }
  if (stratify == "score_covariate") {
    currentSt_score <- paste0("scores_", nGroup, "status")
    currentSt_cov   <- paste0("cov_", nGroup, "status")
    colData(currentData)$score_cov <-
      paste0(colData(currentData)[, currentSt_score],
        " / ",
        colData(currentData)[, currentSt_cov])
    currentStrata <- "score_cov"
    mainTitle <- paste(scoreCol[1], " & ", covariate)
  }
  if (stratify == "expr_expr") {
    currentSt_expr1 <- paste0("expr1_", nGroup, "status")
    currentSt_expr2 <- paste0("expr2_", nGroup, "status")
    colData(currentData)$expr_expr <-
      paste0(colData(currentData)[, currentSt_expr1],
        " / ",
        colData(currentData)[, currentSt_expr2])
    currentStrata <- "expr_expr"
    mainTitle <- paste(gene[1], " & ", gene[2])
  }
  if (stratify == "score_score") {
    currentSt_score1 <- paste0("scores_", nGroup, "status")
    currentSt_score2 <- paste0("scores2_", nGroup, "status")
    colData(currentData)$score_score <-
      paste0(
        colData(currentData)[, currentSt_score1],
        " / ",
        colData(currentData)[, currentSt_score2])
    currentStrata <- "score_score"
    mainTitle <- paste(scoreCol[1], " & ", scoreCol[2])
  }
  
  ##----------------------------- Fit survival model
    tt <- data.frame(table(colData(currentData)[, currentStrata]))
    tt$Var1 <- as.character(tt$Var1)
    tt$Freq <- as.character(tt$Freq)
    
    ## Add number of samples in each group (e.g. high vs low expression)
    for (i in 1:nrow(tt)) {
      colData(currentData)$currentStrata_n[colData(currentData)[, currentStrata] == tt$Var1[i]] <-
        paste0(tt$Var1[i], " (", tt$Freq[i], ")")
    }
  
    fitValues <- survfit(Surv(colData(currentData)[, timeCol],
      as.numeric(as.factor(
        colData(currentData)[, eventCol]
      )) - 1) ~
        colData(currentData)$currentStrata_n)
    
    ss <- survdiff(Surv(colData(currentData)[, timeCol],
      as.numeric(as.factor(
        colData(currentData)[, eventCol]
      )) - 1) ~
        colData(currentData)$currentStrata_n)
    
  ##------------------------------- Calculate p-value
  ## This does not adjust for any covariates
    
  pval <- ifelse (is.na(ss), next, (round(1 - pchisq(
    ss$chisq, length(ss$n) - 1
  ), 6)))[[1]]
  
  
  ##------------------------------ Plot survival curve
  cols <- c(brewer.pal(9, "Set1")[c(2, 3, 4, 5, 7, 8)],
    brewer.pal(8, "Dark2")[c(8, 1, 4, 6)])
  
  p <- autoplot(fitValues, surv.size = 1.5, conf.int = confInt) + 
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) +
    ggtitle(paste0(mainTitle,
                   # " (Chisq = ", round(ss$chisq, 3),
                   " (p = ", pval, ")")) +
    ylab("Survival") +
    xlab("Time") +
    theme_bw()
  
  p
  
}










