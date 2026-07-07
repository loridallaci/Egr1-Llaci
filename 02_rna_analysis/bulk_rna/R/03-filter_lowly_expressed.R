#############
# plot_density_countData
#############
#
# Reads that are rarely detected are less reliable, because the presence or 
#    absence of the gene count may be due more to chance or sequencing effort 
#    than real biology.  Removing these genes from subsequent analysis runs the 
#    risk of false negative results, but at the benefit of fewer false positives 
#    and increased statistical power to make claims on more highly sequenced genes.
#
# This is script to visualize read distributions of DGELists before and after
#    removing lowly expressed genes
#
# Usage:
#
#
# Arguments:
#
#
# Return

plot_density_count_data <- function(
    myDGEList,
    myTitle="Unfiltered counts",
    noYLim=TRUE,
    plotZoom=FALSE,
    yUpper=0.10
){
  par(mar=c(4.1,4.1,2.1,1))
  
  # convert the counts to log base 2
  nsamples <- ncol(myDGEList)
  lcpm <- cpm(myDGEList,log = TRUE)
  col <- rainbow(nsamples)
  
  # determine the lcpm cutoff used by limma filterByExpr
  L <- mean(myDGEList$samples$lib.size) * 1e-6
  M <- median(myDGEList$samples$lib.size) * 1e-6
  
  # Calculate yUpper if noYLim flag passed. Else use passed yUpper
  if(noYLim){
    # Get a maximum ylim
    temp <- vector(mode="numeric",length=nsamples)
    for(i in 1:nsamples){temp[i] <- max(density(lcpm[,i])$y)}
    yUpper <- range(temp)[2] + 0.025
  }
  
  if(plotZoom){
    plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,yUpper), xlim=c(-1,10),las=2, 
         main="", xlab="")
  }
  else{
    plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,yUpper), las=2, 
         main="", xlab="")
  }
  
  title(main=myTitle,xlab="log2(CPM)",ylab = "Density")
  abline(v=log2(10/M + 2/L), lty=3)
  
  for (i in 2:nsamples){
    den <- density(lcpm[,i])
    lines(den$x, den$y, col=col[i], lwd=2)
  }
}