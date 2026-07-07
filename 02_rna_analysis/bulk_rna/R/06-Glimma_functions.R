# makeVolanoPlot_from_eBayes
#
# Updated 10/24/23 to generate volcano plot from -log10(p value) instead of FDR. 
#   Using argument yAxisPlot to replace plotLogFDR
#
#
# Purpose: given a gene expression set tested with eBayes, and optional subsetting, 
#   the subroutine will generate a volcano plot of the desired genes
#
# Input:
#  fittedObj: the output of eBayes in limma package. Either array or voom(rna seq)
#  contrastCoef: the index of the contrast of interest. 1 by default (the value if only 1 contrast)
#  yAxisPlot: String. Default PVal. Allowed "PVal", "FDR", "BStat"
#      plot -log(P-value,base=10) or -log(FDR,base=10) on y axis. 
#      Alternatively, if BStat, will plot the limma B-statistic ($lods) on y axis. 
#  countMatrix: appropriately normalized count matrix. e.g exprs(eset) of RMA normalized arrays, v$E of voom treated RNA-seq
#  contrastMatrix: the contrastMatrix used to call eBayes. Alternatively, define keepGroups
#  keepGroups: a vector with the names of the elements in groupVector which should be retained. 
#      Define this instead of contrastMatrix if for some reason there is no contrastMatrix
#  folderTitle: desired full foldername for output folder (will create a subdirectory in working directory)
#  DEResults: the output of decideTests with desired thresholds
#  groupFactor: contains group assignments used to generate the design matrix.
#  subsetVector: either NULL (default) or a TRUE/FALSE vector of rows to keep. If NULL, all rows plotted. This is intended to capture
#    functional suites, e.g. genes annotated as transcription factors
#  annotationDF: either NULL (default) or a passed data frame. If NULL, will create from fittedObj$genes.
#  launchFLAG: either FALSE (default) or TRUE to indicate if Glimma should open the created volcano plot in a browser
#
# Output:
#
# Usage:
#  makeVolcanoPlot_from_eBayes(subsetVector=subsetVector,annotationDF=efit$genes,DEList=decideTests(efit),fittedObj=efit,
#    contrastCoef=1,folderTitle="Out",myVoom=v,contrastMatrix=contr.matrix)
makeVolcanoPlot_from_eBayes <- function(
    fittedObj,
    contrastCoef=1,
    yAxisPlot = "PVal", # accepted PVal, FDR, BStat
    countMatrix,
    contrastMatrix,
    keepGroups=NULL,
    folderTitle,
    DEresults,
    groupFactor,
    subsetVector=NULL,
    annotationDF=NULL,
    launchFlag = FALSE
)
{
  require(Glimma)
  
  # Grab and format the DE testing results
  ## get the DE results (1,0,-1 from the decideTest object)
  myDt <- DEresults[,contrastCoef]
  ## Translate DE status to words
  myDE <- c("Downregulated", "NoDifference", "Upregulated")[as.factor(myDt)]
  
  # Make the annotation dataframe, including DE status
  ## If annotationDF is null, use the fittedObject genes slot
  if(is.null(annotationDF)){
    cat("annotationDF was not defined, generating annotation data frame from fitted object\n\n")
    annotationDF <- fittedObj$genes
  }
  
  ## Check that SYMBOL is included as an entry for side figure plotting
  if(!("SYMBOL"%in% colnames(annotationDF))){
    stop(cat("\n\nFatal error in makeVolcanoPlot_from_eBayes: the annotationDF must ",
             "contain a column titled SYMBOL (or the $genes slot of the ebayes output ",
             "must contain such a column).\n\n"))
  }
  ##Build the anno table for plotting. Glimma requires a column titled GeneID
  anno <- cbind("GeneID"=row.names(annotationDF),annotationDF,"DE"=myDE)
  
  # Subset the values, if desired
  ## If subsetVector is null, keep all the rows
  if(is.null(subsetVector)){
    subsetVector <- rep(TRUE,dim(fittedObj)[1])
  }
  
  # Subset the fitted object, DE factor, annotation DF, eset
  ## identify the groups to keep
  if(is.null(keepGroups)){
    keepGroups <- row.names(contrastMatrix)[contrastMatrix[,contrastCoef]!=0]
  }
  subsetFit <- fittedObj[subsetVector,]
  subsetAnno <- anno[subsetVector,]
  # commented this out on 12/20/2020
  #  subsetDt <- myDt[subsetVector]  
  subsetDt <- myDt[subsetVector,]
  
  # subset counts to only keep the rows in subsetVector and samples from contrast of interest
  subsetCount <- countMatrix[subsetVector,]
  subsetCount <- subsetCount[,groupFactor %in% keepGroups]
  
  # Limit the groupFactor to only those groups in the contrast of interest
  myGroupFactor <- groupFactor[groupFactor %in% keepGroups]
  myGroupFactor <- droplevels.factor(myGroupFactor)
  
  # generate the plot
  if(!yAxisPlot %in% c("PVal","FDR","BStat")){
    stop(cat("\n\nFatal error in makeVolcanoPlot_from_eBayes: the yAxisPlot value ",
             "must be PVal, FDR, or BStat. PVal is default and recommended.\n\n"))
  } else if (yAxisPlot == "PVal"){
    yAxis = topTable(subsetFit,coef = contrastCoef,number = dim(subsetFit)[1],
                               sort.by = "none")$P.Val
    yAxisLabel = "-log10(P-value)"
    
    glXYPlot(x=subsetFit$coefficients[,contrastCoef],y=-log(yAxis,base = 10),
             status=subsetDt,anno = subsetAnno,main=colnames(subsetFit)[contrastCoef],
             xlab = "log2FC", ylab= yAxisLabel, folder=folderTitle, counts=subsetCount,
             groups=myGroupFactor,samples =colnames(subsetCount), 
             side.main="SYMBOL", side.ylab="Log2(cpm)",launch = launchFlag,
             display.columns = colnames(subsetAnno))
  } else if (yAxisPlot == "FDR") {
    yAxis = topTable(subsetFit,coef = contrastCoef,number = dim(subsetFit)[1],
                               sort.by = "none")$adj.P.Val
    
    yAxisLabel = "-log10(FDR)"
    
    glXYPlot(x=subsetFit$coefficients[,contrastCoef],y=-log(yAxis,base = 10),
             status=subsetDt,anno = subsetAnno,main=colnames(subsetFit)[contrastCoef],
             xlab = "log2FC", ylab= yAxisLabel, folder=folderTitle, counts=subsetCount,
             groups=myGroupFactor,samples =colnames(subsetCount), 
             side.main="SYMBOL", side.ylab="Log2(cpm)",launch = launchFlag,
             display.columns = colnames(subsetAnno))
    
  } else {
    # plot the B statistic in lods
    glXYPlot(x=subsetFit$coefficients[,contrastCoef],y=subsetFit$lods[,contrastCoef],
             status=subsetDt,anno = subsetAnno,main=colnames(subsetFit)[contrastCoef],
             xlab="logFC", ylab="logodds",folder=folderTitle, counts=subsetCount,
             groups=myGroupFactor,samples =colnames(subsetCount), 
             side.main="SYMBOL",side.ylab="Log2(cpm)",launch = launchFlag,
             display.columns = colnames(subsetAnno))
  }
  
}