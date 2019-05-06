
##======================= Author =======================
## ------------ Momeneh (Sepideh) Foroutan -------------
## ------------ Last updated: 03 May 2019 --------------
##======================================================

## This function takes a long format data, with columns for genes, normalised expression values and sample/cell annotations. Then, it generates a vertical plot for boxplots comparing expression of genes across sample types, coloured by specified annotations.

##----------------- INPUT ARGUMENTS --------------------

## data: a data frame contatining at least three columns for gene identifiers, gene expression values and annotation.
## x.axis.column: Name of the column containing expression values, this will be the x axis in the plot
## y.axis.column: Name of the column containing gene identifiers, this will be the y axis in the plot
## annot.column: Name of the column containing sample annotations, this will be used to colour boxplots
## title: Character; title of plot.
## cols: colours to be used to colour boxplots
## textSize: size of texts in the plot

##--------------------- FUNCTION -----------------------

boxplot_vertical_ggplot <- function(data,
  x.axis.column = "LogTPM",
  y.axis.column = "Genes",
  annot.column = "Non.malignant",
  title = "",
  cols = brewer.pal(8, "Set1"),
  textSize = 1.5) {
  
  legend_title_size <- 8 * textSize
  legend_text_size  <- 7 * textSize
  
  p <-  ggplot(data = data) +
    geom_boxplot(
      aes_string(x = y.axis.column,
        y = x.axis.column,
        col = annot.column), 
      size = 0.5,
      position = position_dodge(width = 0.6),
      width = 1.1,
      alpha = 0.8,
      outlier.colour = NULL,
      outlier.size = 0.2 # not actually needed
      ) +
    # geom_point(data = data[data[, x.axis.column] > data$upper.limit |
    #     data[, x.axis.column]  < data$lower.limit, ],
    #   aes_string(
    #     x = y.axis.column,
    #     y = x.axis.column,
    #     color = annot.column
    #   )) +
    
    theme_bw() +
    # scale_y_reverse()+
    coord_flip() +
    scale_colour_manual(values = cols, name = "Cell types") +
    
    ylab ("log(TPM)") +
    xlab("") +
    ggtitle(title) +
    theme(
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = rel(textSize)),
      axis.text.x = element_text(angle = 0, size = rel(textSize)),
      axis.text.y = element_text(
        angle = 0,
        size = rel(textSize),
        face = "italic"
      ),
      axis.line = element_line(colour = "black"),
      axis.ticks = element_line(),
      legend.margin = margin(unit(0, "cm")),
      legend.title = element_text(face = "italic", size = legend_title_size),
      legend.text = element_text(size = legend_text_size),
      legend.key.size = unit(3, 'lines'),
      plot.title = element_text(
        face = "bold",
        size = rel(textSize),
        hjust = 0.5
      )
    ) +
    guides(fill = guide_legend(reverse = TRUE))
  
  return(p)
}


