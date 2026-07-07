# 07-TF_Ligand_Receptor_functions.R
#
# Updated 10/24/23 to reflect yAxisPlot change in volcano plotting (see 06-Glimma_functions)


# Note - THIS CODE BASE ONLY WORKS ON BRIAN's DESKTOP. 
# NEED TO THINK OF WAY TO MAKE THIS MORE GENERAL

###############
# Subsetting functions
###############

# This function takes as input 2 files, an annotation data frame for all genes in the fitted object, 
#    and the data frame for genes of interest
# Both dataframes must have a column titled "ENTREZID"
# It will return a logical vector of length==number of genes to indicate if that probe should be retained for plotting
# @Usage: subsetVector <- make_subsetVector(allGenes=efit$genes,matchGenes=AnimalFT,verbose=FALSE)
# includes check to remove NA values. Not implemented before Taka analysis Dec 4
make_subsetVector <- function(
    allGenes,
    matchGenes,
    verbose=FALSE
)
{
  if(!("ENTREZID"%in% colnames(allGenes))){
    stop(cat("The object allGenes",allGenes,"passed to make_subsetProbe needs a column named 'ENTREZID'\n\n"))
  }
  if(!("ENTREZID"%in% colnames(matchGenes))){
    stop(cat("The object matchGenes",matchGenes,"passed to make_subsetProbe needs a column named 'ENTREZID'\n\n"))
  }
  
  entrezMatch <- allGenes$ENTREZID %in% matchGenes$ENTREZID
  # in case NA was passed and there was an NA in matchGenes
  entrezMatch[is.na(allGenes$ENTREZID)] <- FALSE
  
  if(verbose){
    cat("Number of genes in allGenes:",dim(allGenes)[1],"\n")
    cat("Number of unique Entrez in matchGenes:",length(unique(matchGenes$ENTREZID)),"\n")
    cat("Number of genes from allGenes with match in matchGenes:",sum(entrezMatch),"\n\n")
  }
  
  return(entrezMatch)
}




############
# Animal Transcription Factor Database
############

# genesTable is, for instance, the value of fit2$genes. Must contain column ENTREZID
subset_by_AnimalTF <- function(
    genesTable
)
{
  # load the transcription factor flat file
  animalTF.df <- read.delim(file="~/Data/TranscriptionFactors/AnimalTFDB/AnimalTfDb_annotated.txt",header=TRUE,sep="\t",stringsAsFactors = FALSE)
  subsetVector <- make_subsetVector(allGenes = genesTable,matchGenes = animalTF.df)
  return(subsetVector)
}

############
# Toronto Transcription Factor Database
############

# This is a human TF database
# genesTable is, for instance, the value of fit2$genes. Must contain column ENTREZID
subset_by_TorontoTF <- function(
    genesTable
)
{
  # load the transcription factor flat file
  torontoTF.df <- read.delim(file="~/Data/Human_GeneSets/TranscriptionFactors/Toronto_HumanTF_Formatted.txt",header=TRUE,sep="\t",stringsAsFactors = FALSE)
  subsetVector <- make_subsetVector(allGenes = genesTable,matchGenes = torontoTF.df)
  return(subsetVector)
}

############
# Fantom Mouse Ligand/Receptor
############

subset_by_FantomMouseLigand <- function(
    genesTable
)
{
  # load the mouse Ligand Flat Files from Fantom
  mouseLigand <- read.delim(file="~/Data/GeneSets/Fantom/Fantom_MouseLigand_Annotated.txt",header=TRUE,sep="\t",stringsAsFactors = FALSE)
  # change the column title from entrezgene to ENTREZID for downstream
  colnames(mouseLigand)[4] <- "ENTREZID"
  subsetVector <- make_subsetVector(allGenes = genesTable,matchGenes = mouseLigand)
  return(subsetVector)
}

subset_by_FantomMouseReceptor <- function(
    genesTable
)
{
  # load the mouse Receptor Flat Files from Fantom
  mouseReceptor <- read.delim(file="~/Data/GeneSets/Fantom/Fantom_MouseReceptor_Annotated.txt",header=TRUE,sep="\t",stringsAsFactors = FALSE)
  # change the column title from entrezgene to ENTREZID for downstream
  colnames(mouseReceptor)[4] <- "ENTREZID"
  subsetVector <- make_subsetVector(allGenes = genesTable,matchGenes = mouseReceptor)
  return(subsetVector)
}

############
# Fantom Human Ligand/Receptor
############

subset_by_FantomHumanLigand <- function(
    genesTable
)
{
  # load the Human Ligand Flat Files from Fantom
  humanLigand <- read.delim(file="~/Data/GeneSets/Fantom/Fantom_HumanLigand_Annotated.txt",header=TRUE,sep="\t",stringsAsFactors = FALSE)
  # change the column title from entrezgene to ENTREZID for downstream
  colnames(humanLigand)[4] <- "ENTREZID"
  subsetVector <- make_subsetVector(allGenes = genesTable,matchGenes = humanLigand)
  return(subsetVector)
}

subset_by_FantomHumanReceptor <- function(
    genesTable
)
{
  # load the Human Receptor Flat Files from Fantom
  humanReceptor <- read.delim(file="~/Data/GeneSets/Fantom/Fantom_HumanReceptor_Annotated.txt",header=TRUE,sep="\t",stringsAsFactors = FALSE)
  # change the column title from entrezgene to ENTREZID for downstream
  colnames(humanReceptor)[4] <- "ENTREZID"
  subsetVector <- make_subsetVector(allGenes = genesTable,matchGenes = humanReceptor)
  return(subsetVector)
}

############
# Amigo Cell Surface, for RC
############

subset_by_AmigoCellSurface <- function(
    genesTable
)
{
  # load the mouse Cell Surface Flat Files from Amigo
  cellSurface <- read.delim(file="~/Data/GeneSets/Amigo/Amigo_MouseCellSurface_Annotated.txt",header=TRUE,sep="\t",stringsAsFactors = FALSE)
  # change the column title from entrezgene to ENTREZID for downstream
  colnames(cellSurface)[4] <- "ENTREZID"
  subsetVector <- make_subsetVector(allGenes = genesTable,matchGenes = cellSurface)
  return(subsetVector)
}  

############
# annotate topTable, mouse
############

geneTable_to_TF_Ligand_Receptor_Table <- function(
    dataframe # must include a column ENTREZID. could be topTable or fit2$genes
)
{
  tf.subset <- subset_by_AnimalTF(dataframe)
  ligand.subset <- subset_by_FantomMouseLigand(dataframe)
  receptor.subset <- subset_by_FantomMouseReceptor(dataframe)
  returnDF <- data.frame("TranscriptionFactor"=tf.subset,"Ligand"=ligand.subset,"Receptor"=receptor.subset)
  row.names(returnDF) <- row.names(dataframe)
  return(returnDF)
}

############
# annotate topTable, human
############

humanGeneTable_to_TF_Ligand_Receptor_Table <- function(
    dataframe # must include a column ENTREZID. could be topTable or fit2$genes
)
{
  tf.subset <- subset_by_TorontoTF(dataframe)
  ligand.subset <- subset_by_FantomHumanLigand(dataframe)
  receptor.subset <- subset_by_FantomHumanReceptor(dataframe)
  returnDF <- data.frame("TranscriptionFactor"=tf.subset,"Ligand"=ligand.subset,"Receptor"=receptor.subset)
  row.names(returnDF) <- row.names(dataframe)
  return(returnDF)
}

############
# Generate subset Glimma plots by coefficient number
############

# makeSubsetGlimmaPlots_by_Coefficient
# Purpose: Generate subset Glimma volcano plots of TF's, receptors, and ligands
#
# Input:
#  fittedObject: the output of eBayes in limma package. Either array or voom(rna seq)
#  contrastCoef: the contract coefficient to plot
#  countMatrix: appropriately normalized count matrix. e.g exprs(eset) of RMA normalized arrays, v$E of voom treated RNA-seq
#  keepGroups: a 2 element vector with names of groups that should be retained / plotted
#  folderTitle: desired full foldername for output folder (will create a subdirectory in working directory)
#  DEResults: the output of decideTests with desired thresholds
#  groupFactor: contains group assignments used to generate the design matrix.
#  annotationDF: either NULL (default) or a passed data frame. If NULL, will create from fittedObj$genes.
#  launchFLAG: either FALSE (default) or TRUE to indicate if Glimma should open the created volcano plot in a browser
#  species: either mouse (default) or human.
#
# Output:
#
# Usage:
#  makeSubsetGlimmaPlots_by_Coefficient(fittedObject=efit,contrastCoef=1,countMatrix=v$E,
#    folderTitle=file.path(deOutput,"Colon_vs_Ileum"),annotationDF=myAnno,
#    DEresults=dt,groupFactor=myGroup)
makeSubsetGlimmaPlots_by_Coefficient <- function(
    fittedObj,
    countMatrix,
    contrastCoef,
    folderTitle,
    keepGroups,
    DEresults,
    groupFactor,
    annotationDF=NULL,
    launchFlag = FALSE,
    species="mouse" # or "human"
)
{
  if(species != "mouse" && species != "human"){
    return(cat("Exiting makeSubsetGlimmaPlots_by_Coefficient subroutine.",
               "Passed argument \'species\' must be \"mouse\" or \"human\".\n\n"))
  }
  
  myStartingDirectory <- getwd()
  dir.create(folderTitle)
  setwd(folderTitle)
  
  if(is.null(annotationDF)){
    annotationDF<- fittedObj$genes[,c("SYMBOL","GENENAME","ENTREZID","BIOTYPE")]
  }
  
  # plot all genes. (there is no passed SubsetVector)
  makeVolcanoPlot_from_eBayes(fittedObj = fittedObj, countMatrix = countMatrix,
                              keepGroups=keepGroups ,DEresults = DEresults,
                              groupFactor = groupFactor,annotationDF = annotationDF,
                              contrastCoef = contrastCoef, folderTitle="AllGenes_Volcano")
  
  if(species=="mouse"){
    # Mouse TF
    subsetVector <- subset_by_AnimalTF(fittedObj$genes)
    makeVolcanoPlot_from_eBayes(fittedObj = fittedObj , countMatrix = countMatrix,
                                keepGroups= keepGroups,DEresults = DEresults,
                                groupFactor = groupFactor, annotationDF = annotationDF,
                                contrastCoef = contrastCoef, 
                                subsetVector = subsetVector,folderTitle="Mouse_TranscriptionFactors_Volcano")
    
    ## Mouse ligands
    subsetVector <- subset_by_FantomMouseLigand(fittedObj$genes)
    makeVolcanoPlot_from_eBayes(fittedObj = fittedObj , countMatrix = countMatrix,
                                keepGroups= keepGroups,DEresults = DEresults,
                                groupFactor = groupFactor, annotationDF = annotationDF,
                                contrastCoef = contrastCoef, 
                                subsetVector = subsetVector,folderTitle="Mouse_Ligands_Volcano")
    
    ## Mouse Receptors
    subsetVector <- subset_by_FantomMouseReceptor(fittedObj$genes)
    makeVolcanoPlot_from_eBayes(fittedObj = fittedObj , countMatrix = countMatrix,
                                keepGroups= keepGroups,DEresults = DEresults,
                                groupFactor = groupFactor, annotationDF = annotationDF,
                                contrastCoef = contrastCoef, 
                                subsetVector = subsetVector,folderTitle="Mouse_Receptors_Volcano")
  }
  if(species=="human"){
    # Human TF
    subsetVector <- subset_by_TorontoTF(fittedObj$genes)
    makeVolcanoPlot_from_eBayes(fittedObj = fittedObj , countMatrix = countMatrix,
                                keepGroups= keepGroups,DEresults = DEresults,
                                groupFactor = groupFactor, annotationDF = annotationDF,
                                contrastCoef = contrastCoef, 
                                subsetVector = subsetVector,folderTitle="Human_TranscriptionFactors_Volcano")
    
    ## Human ligands
    subsetVector <- subset_by_FantomHumanLigand(fittedObj$genes)
    makeVolcanoPlot_from_eBayes(fittedObj = fittedObj , countMatrix = countMatrix,
                                keepGroups= keepGroups,DEresults = DEresults,
                                groupFactor = groupFactor, annotationDF = annotationDF,
                                contrastCoef = contrastCoef, 
                                subsetVector = subsetVector,folderTitle="Human_Ligands_Volcano")
    
    ## Human Receptors
    subsetVector <- subset_by_FantomHumanReceptor(fittedObj$genes)
    makeVolcanoPlot_from_eBayes(fittedObj = fittedObj , countMatrix = countMatrix,
                                keepGroups= keepGroups,DEresults = DEresults,
                                groupFactor = groupFactor, annotationDF = annotationDF,
                                contrastCoef = contrastCoef, 
                                subsetVector = subsetVector,folderTitle="Human_Receptors_Volcano")  
  }
  setwd(as.character(myStartingDirectory))
}