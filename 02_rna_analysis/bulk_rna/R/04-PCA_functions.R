#PCA_functions.R

### Bundle ggplot layers in a pca Theme

pca_personal_theme <- function(){
  theme_classic() +
  theme(panel.background = element_rect(colour = "black", linewidth = 1),
        plot.title = element_text(size=rel(1.5), hjust = 0.5),
        axis.title = element_text(size=rel(1.25)),
        axis.text = element_text(size=rel(1.15)),
        legend.title = element_text(size=rel(1.25)),
        legend.text= element_text(size=rel(1.15))) 
}

calculate_PCA <- function(
  myDGEList,  
  numberTopGenes = 2000){
  
  lcpm <- cpm(myDGEList,log = TRUE)
  topVarGenes <- head(order(rowVars(lcpm), decreasing = TRUE),n = numberTopGenes)
  edata.PCA <- lcpm[topVarGenes,]

  # Use base prcomp for PCA calculation
  res.pca <- prcomp(t(edata.PCA), scale = TRUE)
}

pca_axis_label <- function(
  pca_result,
  pca_axis,
  sig_digits=1
  ){
    #percent of variation on each PC axis
    eigs <- pca_result$sdev^2
    
    propVar <- 100*(eigs / sum(eigs))
    
    return(paste0("PC",pca_axis,": ", round(propVar[pca_axis],sig_digits), "% variance"))
}



