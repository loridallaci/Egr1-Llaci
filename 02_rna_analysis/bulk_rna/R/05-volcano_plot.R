# Updated 10/24/2023 to plot -log(p val) instead of -log(FDR)


volcano_plot_from_Ebayes <- function(
    myEBayes,
    myContrastCoef,
    myTitle,
    myMinLFC=1,
    myMaxPVal=0.05
){
  # basic volcano plot
  cexValue=0.8
  
  limma_pvals = topTable(myEBayes, coef = myContrastCoef, number = nrow(myEBayes),
                             sort.by = "none")$P.Val
  ebayesResults <- data.frame("Log2FC"= myEBayes$coefficients[,myContrastCoef],
                              "LogPVal"= -log(limma_pvals,base=10),
                              "SYMBOL"= myEBayes$genes$SYMBOL)
  
  # Basic Volcano Plot
  with(ebayesResults, plot(Log2FC, LogPVal, pch=20, main=myTitle, xlab=expression('Log'[2]*' fold change'),ylab=expression('-Log'[10]*' P-value')))
  
  # Add colored points: red if upregulated, blue if downregulated
  with(subset(ebayesResults, Log2FC>=myMinLFC & LogPVal >=-log(myMaxPVal,base=10)), points(Log2FC, LogPVal, pch=20, col="red"))
  with(subset(ebayesResults, Log2FC<=-myMinLFC & LogPVal >=-log(myMaxPVal,base=10)), points(Log2FC, LogPVal, pch=20, col="dodgerblue"))
}