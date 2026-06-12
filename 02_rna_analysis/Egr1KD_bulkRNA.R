********************************** NEW ANALYSIS *********************************************
  091625
# Keep guides for Egr1 separate. How does that change overlap with CC data?
AGAIN
1. Determine which samples are outliers first:
  ```{r parameters}
# Save for another time example of how to set working directory with relative path in RMD
# original (author's machine): "C:/Users/loril/Documents/Egr1/Egr1_bulk_RNAseq/Egr1_KD_030525/"
myWorkingDirectory = "."
figureOutputDirectory = file.path(myWorkingDirectory,"figures_Egr1KD_again_keepEgr1GuidesSeparate_only0filteredout_removeOutliers")
analysisOutputDirectory = file.path(myWorkingDirectory,"analysis_Egr1KD_again_keepEgr1GuidesSeparate_only0filteredout_removeOutliers")
metadataDirectory = file.path(myWorkingDirectory,"metadata")
dataDirectory = file.path(myWorkingDirectory,"data")
deOutputDirectory = file.path(myWorkingDirectory,"diffexp_Egr1KD_again_keepEgr1GuidesSeparate_only0filteredout_removeOutliers")

## Files
countDataFile <- file.path(dataDirectory,"Dedup_Counts.txt")
geneDataFile <- file.path(dataDirectory,"GeneInfo_Egr1KD.csv")
sampleMetadataFile <- file.path(metadataDirectory,"Meta_data_Egr1KD.csv")
```

### Color Palettes and Graphics
```{r graphics}

# http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
palette_cb_gray <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
palette_cb_black <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Tableau color blind 10
# https://public.tableau.com/views/TableauColors/ColorPaletteswithRGBValues?%3Aembed=y&%3AshowVizHome=no&%3Adisplay_count=y&%3Adisplay_static_image=y
# https://jrnold.github.io/ggthemes/reference/tableau_color_pal.html
palette_tableau_cb10 <- c("#1170aa","#fc7d0b","#a3acb9","#57606c","#5fa2ce","#c85200","#7b848f","#a3cce9","#ffbc79","#c8d0d9")

# palette
palette_ALI_only <- pal_nejm()(5)
palette_All_Data <- pal_nejm()(8)

```

### Values for Differential Expression Testing
```{r}
maxFDR <- 0.05
minLFC <- 0.5
```

## Block 2: Load the count data, sample data, and gene data file from our dataset

```{r, fileLoad}
countData <- read.delim(countDataFile, header = TRUE, 
                        stringsAsFactors = FALSE, row.names = 1)
countData

geneData <- read.delim(geneDataFile,header = TRUE,stringsAsFactors = FALSE,row.names=1, sep = ',')
geneData <- geneData[!duplicated(geneData$ensembl_gene_id), ] 
#geneData4 <- geneData3[,-1]
#rownames(geneData4) <- geneData3[,1]
#geneData <- geneData4
head(geneData)

sampleData <- read.delim(sampleMetadataFile,header=TRUE,row.names=1, sep = ',')
sampleData <- as.data.frame(t(sampleData))
#rownames(sampleData) <- paste0(rownames(sampleData), '.')
colnames(sampleData)[3] <- 'group'
sampleData

countData <- countData[, rownames(sampleData)]
countData

```
Remove outliers, No Treatment Neg2, and combine Egr1KD guides per sex

```{r}
library(dplyr)

# Vector of samples to remove
samples_to_remove <- c(
  "F6_Egr1_gRNA2_Rep2",
  "F6_Egr1_gRNA2_Rep3",
  "F6_Egr1_gRNA2_Rep4"
)

sampleData <- sampleData %>%
  # Step 0: Remove specific unwanted samples
  filter(!(rownames(.) %in% samples_to_remove)) %>%
  
  # Step 1: Remove NoTreatment_gRNA2 samples
  filter(Treatment != "NoTreatment_gRNA2")

# Step 2: Combine Egr1KD_gRNA2 and Egr1KD_gRNA3 into one group per sex
#mutate(
#  Treatment = case_when(
#    Treatment %in% c("Egr1KD_gRNA2", "Egr1KD_gRNA3") ~ "Egr1KD_gRNA",
#    TRUE ~ Treatment
#  ),
#  group = paste(Treatment, Sex, sep = "_")
#)

# Make sure countData matches updated sampleData
countData <- countData[, colnames(countData) %in% rownames(sampleData)]

```


filter countData to only keep my samples that are in the metadata
```{r}
library(dplyr)
countData <- countData[, colnames(countData) %in% rownames(sampleData)] 
```

filter geneData to only keep genes that are in the countData
```{r}
geneData <- geneData[geneData$ensembl_gene_id %in% rownames(countData), ] 
```

```{r}
sampleData[,c('Treatment', 'Sex', 'group')] <- lapply(sampleData[,c('Treatment', 'Sex', 'group')] , factor)
```



```{r}
gene_data <- geneData
count_data <- countData
sample_data <- sampleData
```

```{r}
sample_data$Group <- sample_data$group
```

```{r raw_data_check}
# Ensure that sample_data and count_data have same sampleID's
sample_data <- sample_data[order(row.names(sample_data), decreasing = FALSE),]
count_data <- count_data[,order(colnames(count_data), decreasing = FALSE)]

if(!(all.equal(row.names(sample_data),colnames(count_data)))){
  stop("FATAL ERROR: The sampleID's are not identical in the loaded count and metadata files. Please inspect the raw data.")
}

# check if sample_data contains a Group column
if(!("Group" %in% colnames(sample_data))){
  stop(paste0("FATAL ERROR: A category \"Group\" must be in the Sample Metadata. Please inspect the sample metadata file and try again."))
}

# check if count data begins with ENSEMBL ids
if(sum(grep("^ENS", row.names(count_data), invert = TRUE))!=0){
  stop("FATAL ERROR: Gene identifiers should by Ensembl ID's. One or more of your genes did not beging with \"ENS\". Please check your data.")
}
```

The sample metadata can be inspected in this interactive table.

```{r interactive_sample_table}
datatable(sample_data, class = "compact", options = list(pageLength = 25))
```

## Build Gene Annotation

The counts files are indexed by Ensembl id. We use the bioMart package to expand the annotation to include gene Symbols, Names, and other useful details.

```{r build_gene_anno}
#gene_data <- build_gene_table_from_ensembl(ensembl_ids = row.names(count_data), mySpecies = species)
```

# Data Pre-Processing and Quality Control

We begin by exploring and pre-processing the provided counts data. We use the bioinformatics tools [EdgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) and [Limma](https://bioconductor.org/packages/release/bioc/html/limma.html) to organize and process the data. 

```{r create_limma}
# Use remove.zeros= FALSE to keep congruence with genes file
y <- DGEList(counts = count_data, samples = sample_data,genes = gene_data,remove.zeros = FALSE,group = sample_data$group)
```

<!-- Save the Limma Elist object for future analysis. -->
  
  ```{r save_limma}
#saveRDS(y, file.path(output_directory,"MandF_filteredOutliers.RDS"))
```

## Sequencing Effort and Evenness {.tabset}

If samples aren't sequenced deeply enough, we may lose power to detect lowly expressed genes. In addition, large variability in sequencing depth between different samples can introduce bias that might cloud clustering and differential expression. 

We inspect the sequencing effort for each sample and main biological group.

The target depth of sequencing for BRB-seq 3' sequencing is 2-10 million counts per sample.

### Sequence Effort Boxplot

This plot shows the distribution of counts per sample, organized by Group. Hover over the small circles to identify the sample name.

all samples
```{r sequenceEffortBoxplot}
# Generate the stripchart
stripchart = stripchart_sequence_effort_by_group(y, y$samples$Group)

# print the interactive stripchart to the report
stripchart %>% ggplotly(tooltip = c("text","y"))

# Use the %+% operator to add components to the plot, here stripping the label so we can print as a static image

pdf(file=file.path(figureOutputDirectory,"Sequence_Effort_Boxplot_MandF.pdf"),width = 4, height = 4, useDingbats = TRUE)

stripchart %+% aes(text=NULL)

invisible(dev.off())
```

### Sequence Effort Table

This table shows the distribution of counts per sample, organized by Group.

```{r sequenceEffort}
# Assess heterogeneity of sequencing effort by group
sequence_effort <- calculate_sequence_effort_by_group_table(myDGEList = y,myGroup = y$samples$Group)

datatable(sequence_effort, class = "display", 
          caption = htmltools::tags$caption(
            style = 'caption-side: bottom; text-align: left; color: red;',
            'Sequencing Depth by Group (millions of reads)'))
```

# Filter Lowly Expressed Genes and Normalize {.tabset}

Not all genes are expressed in all cells and tissues. Further, extremely low levels of gene expression might not be biologically meaningful. Removing these lowly expressed genes helps to focus the analysis on genes that are more likely to have functional consequence, and improves the ability to statistically model the variance in the data.

## Visualize the Count Distribution of Raw Data

The following charts show how often a gene with given expression is detected in the raw data. The data is transformed to log2-counts per million to reflect different sequencing depths for different samples. To avoid division by zero, a small "pseudocount" is added to each count.

This data shows a typical distribution, including a large number of genes that are not detected at all or only at very low levels (log2-cpm < 0).

```{r plot_density_unfiltered, fig.width=8}
par(mfrow=c(1,2))

plot_density_count_data(myDGEList = y,myTitle = "A. Raw Data", noYLim = TRUE, plotZoom = FALSE)

plot_density_count_data(myDGEList = y,myTitle = "B. Raw Data, Zoomed", noYLim = FALSE, plotZoom = TRUE, yUpper = 0.06)
```

The Vertical dashed line in these plot shows the log2-CPM cutoff used by Limma's*filterByExpr* function.

## Visualize Data Filtered using Limma Function

```{r filter_and_plot_filterByExpr}
keep.exprs <- filterByExpr(y, group = y$samples$Group)
x <- y[keep.exprs,, keep.lib.sizes = FALSE]

# Make another density plot using filtered counts
plot_density_count_data(myDGEList=x, myTitle = "C. Filtered Data", noYLim = TRUE, plotZoom = FALSE)
```

## Visualize Data Filtered using CPM Thresholds

Here we manually set a CPM threshold and number of samples. Trying this since it is a shallowly sequenced project

```{r filter_and_plot_byThreshold}

#minCpmThreshold <- 2
#minNumSamp <- 4

#keep.exprs.2 <- rowSums(cpm(y,log=TRUE)>=minCpmThreshold)>=minNumSamp
#sum(keep.exprs.2)
  
#x.2 <- y[keep.exprs.2,keep.lib.sizes=FALSE]
#x.2 <- calcNormFactors(x.2)

# Make another density plot using filtered counts
#plot_density_count_data(myDGEList=x.2, myTitle = "C. Filtered Data, Direct", noYLim = TRUE, plotZoom = FALSE)
```

Conclusion: use the manual filtering

## Normalize Gene distributions

Gene expression is normalized to sequencing depth using the trimmed Mean of M-Values (TMM) algorithm from **EdgeR**.  
********************************************* MAKE SURE TO RUN THIS CODE BELOW AS X.2, NOT JUST X *********************************************
IN THIS CASE I'M USING COUNTS FILTERED BY EDGEr ONLY
```{r normalize}
x <- calcNormFactors(x, method = "TMM")
```

<!-- Save the filtered DGE list -->
  ```{r save_filtered_dge}
#saveRDS(x,file=file.path(output_directory,"DGElist_filtered_MandF_filteredOutliers.RDS"))
```

# Sample Heatmap

Heatmaps are an important tool to visualize transcriptomic data. These plots can visually highlight outliers that don't cluster with replicates, and suggest how similar or dissimilar the expression of groups are.

This plot shows a heatmap generated by calculating the Euclidean distance between every pair of log2-cpm normalized samples.
```{r allSampleHeatmap, fig.height=7, fig.width=7}
# Calculate the distance matrix and convert to data.frame
sample_dists <- dist(t(cpm(x, log=TRUE)), method="euclidean")
sample_dist_DF <- as.data.frame(as.matrix(sample_dists))

# Prepare annotation data frame
annotation_DF <- x$samples[, colnames(x$samples) %in% c("Group", "Sex")]
rownames(annotation_DF) <- rownames(sample_dist_DF)

# Drop unused levels to avoid color mismatch
annotation_DF$Group <- droplevels(factor(annotation_DF$Group))
annotation_DF$Sex   <- droplevels(factor(annotation_DF$Sex))

# Define colors for annotations based on *exact* present levels
anno_colors <- list(
  Sex = c(
    Female = "#F08080",     # light pink
    Male = "darkgreen"      # green
  ),
  Group = c(
    Egr1KD_gRNA2_Female = "#6A0DAD",
    Egr1KD_gRNA3_Female = "#B266FF",
    Egr1KD_gRNA2_Male = "#FF8C00",
    Egr1KD_gRNA3_Male = "#FFB84D",
    #NoTreatment_gRNA2_Female = "#FF69B4",
    #NoTreatment_gRNA2_Male = "#1E90FF",
    NoTreatment_gRNA1_Male = "#1E90FF",
    NoTreatment_gRNA1_Female = "#FF69B4"
  )
)
# Generate heatmap
p <- pheatmap(
  sample_dist_DF,
  col = colorRampPalette(c("blue", "white", "red"))(20),
  border_color = "black",
  fontsize = 10,
  fontsize_col = 7,
  show_rownames = FALSE,
  show_colnames = TRUE,
  annotation_col = annotation_DF,
  annotation_colors = anno_colors,
  clustering_distance_rows = sample_dists,
  clustering_distance_cols = sample_dists
)

# Save to PDF
pdf(file = file.path(figureOutputDirectory, 
                     "SampleSample_MandF_filteredOutliers_Heatmap.pdf"), 
    width = 10, height = 8)
p
invisible(dev.off())
p
```
```{r}
top_number <- 2000
pca_results <- calculate_PCA(x,numberTopGenes= top_number)

# Make a data.frame with PC1, PC2, and metadata 
pca_data <- data.frame("PC1"=pca_results$x[,1],
                       "PC2"=pca_results$x[,2],
                       "Group"=x$samples$Group,
                       "Sample"=colnames(x))

# Factor to keep group order congruent with heatmaps
pca_data$Group <- factor(pca_data$Group, levels = levels(x$samples$Group))
group_colors <- c(
    Egr1KD_gRNA2_Female = "#6A0DAD",
    Egr1KD_gRNA3_Female = "#B266FF",
    Egr1KD_gRNA2_Male = "#FF8C00",
    Egr1KD_gRNA3_Male = "#FFB84D",
    #NoTreatment_gRNA2_Female = "#FF69B4",
    #NoTreatment_gRNA2_Male = "#1E90FF",
    NoTreatment_gRNA1_Male = "#1E90FF",
    NoTreatment_gRNA1_Female = "#FF69B4"
)



levels(factor(pca_data$Group))
names(group_colors)

ggPretty <- ggplot(pca_data, aes(PC1, PC2, text = Sample)) + 
  geom_point(aes(fill = Group), shape = 21, size = 4) +
  scale_fill_manual(values = group_colors) +  # 👈 custom colors here
  xlab(pca_axis_label(pca_results, pca_axis = 1)) +
  ylab(pca_axis_label(pca_results, pca_axis = 2)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.75) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.75) +
  ggtitle(label = "RNA-Seq", subtitle = paste0(top_number, " most variable probes")) +
  pca_personal_theme() +
  coord_fixed()

pdf(file=file.path(figureOutputDirectory,"PCA_MandF_filteredOutliers.pdf"),width = 7,height = 7,onefile = TRUE,useDingbats = FALSE)
ggPretty
invisible(dev.off())


ggPretty
```
```{r}
library(ggrepel)

# Female samples
pca_data_female <- pca_data[grepl("Female", pca_data$Group), ]

# Male samples
pca_data_male <- pca_data[grepl("Male", pca_data$Group), ]

# Female PCA plot
ggFemale <- ggplot(pca_data_female, aes(PC1, PC2, text = Sample)) + 
  geom_point(aes(fill = Group), shape = 21, size = 4) +
   geom_text_repel(aes(label = Sample), size = 2.5) +
  scale_fill_manual(values = group_colors) +
  xlab(pca_axis_label(pca_results, pca_axis = 1)) +
  ylab(pca_axis_label(pca_results, pca_axis = 2)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.75) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.75) +
  ggtitle("RNA-Seq PCA (Females)", subtitle = paste0(top_number, " most variable probes")) +
  pca_personal_theme() +
  coord_fixed()

# Male PCA plot
ggMale <- ggplot(pca_data_male, aes(PC1, PC2, text = Sample)) + 
  geom_point(aes(fill = Group), shape = 21, size = 4) +
   geom_text_repel(aes(label = Sample), size = 2.5) +
  scale_fill_manual(values = group_colors) +
  xlab(pca_axis_label(pca_results, pca_axis = 1)) +
  ylab(pca_axis_label(pca_results, pca_axis = 2)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.75) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.75) +
  ggtitle("RNA-Seq PCA (Males)", subtitle = paste0(top_number, " most variable probes")) +
  pca_personal_theme() +
  coord_fixed()

pdf(file = file.path(figureOutputDirectory, "PCA_Females_filteredOutliers.pdf"), width = 7, height = 7)
print(ggFemale)
dev.off()

pdf(file = file.path(figureOutputDirectory, "PCA_Males_filteredOutliers.pdf"), width = 7, height = 7)
print(ggMale)
dev.off()




```


Plot males in Pc2 and PC3
```{r}
top_number <- 2000
pca_results <- calculate_PCA(x,numberTopGenes= top_number)

# Make a data.frame with PC1, PC2, and metadata 
pca_data <- data.frame("PC1"=pca_results$x[,1],
                       "PC2"=pca_results$x[,2],
                       "PC3"=pca_results$x[,3],
                       "Group"=x$samples$Group,
                       "Sample"=colnames(x))

# Factor to keep group order congruent with heatmaps
pca_data$Group <- factor(pca_data$Group, levels = levels(x$samples$Group))


levels(factor(pca_data$Group))
names(group_colors)

ggPretty <- ggplot(pca_data, aes(PC2, PC3, text = Sample)) + 
  geom_point(aes(fill = Group), shape = 21, size = 4) +
  scale_fill_manual(values = group_colors) +  # 👈 custom colors here
  xlab(pca_axis_label(pca_results, pca_axis = 2)) +
  ylab(pca_axis_label(pca_results, pca_axis = 3)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.75) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.75) +
  ggtitle(label = "RNA-Seq", subtitle = paste0(top_number, " most variable probes")) +
  pca_personal_theme() +
  coord_fixed()

pdf(file=file.path(figureOutputDirectory,"PCA_MandF_PC2and3_filteredOutliers.pdf"),width = 7,height = 7,onefile = TRUE,useDingbats = FALSE)
ggPretty
invisible(dev.off())


ggPretty
```
```{r}
# Female samples
pca_data_female <- pca_data[grepl("Female", pca_data$Group), ]

# Male samples
pca_data_male <- pca_data[grepl("Male", pca_data$Group), ]

# Female PCA plot
ggFemale <- ggplot(pca_data_female, aes(PC2, PC3, text = Sample)) + 
  geom_point(aes(fill = Group), shape = 21, size = 4) +
   geom_text_repel(aes(label = Sample), size = 2.5) +
  scale_fill_manual(values = group_colors) +
  xlab(pca_axis_label(pca_results, pca_axis = 2)) +
  ylab(pca_axis_label(pca_results, pca_axis = 3)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.75) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.75) +
  ggtitle("RNA-Seq PCA (Females)", subtitle = paste0(top_number, " most variable probes")) +
  pca_personal_theme() +
  coord_fixed()

# Male PCA plot
ggMale <- ggplot(pca_data_male, aes(PC2, PC3, text = Sample)) + 
  geom_point(aes(fill = Group), shape = 21, size = 4) +
   geom_text_repel(aes(label = Sample), size = 2.5) +
  scale_fill_manual(values = group_colors) +
  xlab(pca_axis_label(pca_results, pca_axis = 2)) +
  ylab(pca_axis_label(pca_results, pca_axis = 3)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.75) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.75) +
  ggtitle("RNA-Seq PCA (Males)", subtitle = paste0(top_number, " most variable probes")) +
  pca_personal_theme() +
  coord_fixed()

pdf(file = file.path(figureOutputDirectory, "PCA2_Females_filteredOutliers.pdf"), width = 7, height = 7)
print(ggFemale)
dev.off()

pdf(file = file.path(figureOutputDirectory, "PCA2_Males_filteredOutliers.pdf"), width = 7, height = 7)
print(ggMale)
dev.off()


```



## Block 3: Create a DESeq2 object and perform basic pre-filtering

Use the three files from our minicourse to make a DESEq object called “dds”. 

The design should be “~0+Group” to standardize our ability to compare groups later.

### Create a DESeq2 object 

Now feed the CPM filtered output to DESEQ2

```{r}
library("DESeq2")

# This is actually all it takes to make  basic object
#dds <- DESeqDataSetFromMatrix(countData = countData,
#                              colData = sampleData,
#                              design=~0+group)

# *************************************** CHANGING ALL THIS BELOW TO MATCH THE ORIGINAL METHOD *********************************************

# Step 3: Extract filtered raw count matrix
#filtered_counts <- x$counts

# Step 4: Ensure sampleData rownames match filtered_counts colnames
# (if not already done)
#sampleData <- sampleData[colnames(filtered_counts), ]

# Step 5: Build DESeq2 object with filtered data
#dds <- DESeqDataSetFromMatrix(countData = filtered_counts,
#                              colData = sampleData,
#                              design = ~ group)  # Replace with your design

# Step 6: Run DESeq2
#dds <- DESeq(dds)

# Step 7: Get results
#res <- results(dds)

# View results
#head(res[order(res$pvalue), ])

library("DESeq2")

# This is actually all it takes to make  basic object
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = sampleData,
                              design=~0+group)
```

### Store the gene Data as a metadata feature

Here is a line of code that will store the values of geneData into the metadata of dds

```{r}
# First: make sure geneData rownames are gene IDs (so they match rownames(dds))
geneData
rownames(geneData) <- geneData$ensembl_gene_id
geneData
# Then: subset and reorder geneData to match dds
geneData <- geneData[rownames(dds), ]

# Now assign it safely to mcols
mcols(dds) <- DataFrame(mcols(dds), geneData)
```

### Minimal prefiltering

Here we will remove any genes that have values in only a single sample (Michael Love's approach in the workflow).

```{r}
dds2 <- dds[rowSums(counts(dds)) > 1, ] # select non zero
# scale by Sequencing effort, because we have removed some rows of data
dds2 <- estimateSizeFactors(dds2)
dds2
dim(dds2)
```

In practice, I often use a more stringent criteria that requires a minimum count in at least a certain number of samples (where the number of samples is proportional to the smallest group category). We will build up these approaches later with real data.

### Calculate the regularized logarithmic transformation of the data

There are two baked in methods to transform or normalize the data; regularized logarithms (rlog) and variance stabilizing transformations (vst). You can use these instead of log2 transforms to normalize data before plotting PCA, heatmap, etc.

Note that when you actually run the differential expression test, you will do that on the plain counts, and the differential expression function will do the appropriate transformations and model fitting.

```{r}
rld <- vst(dds2, blind = FALSE)
rld
```

## Block 4: Exploratory Data Analysis and QC Plots

### Assess sequencing effort by group

We want to know if any samples are sequenced at a much different depth than the rest of the samples

To understand this code, look up the aggregate() function from the base "stats" package.

```{r sequenceEffort}
# Assess heterogeneity of sequencing effort by group (currently each collection of three replicates)
depthBygroup <- data.frame(
  "Min"=aggregate(colSums(assay(dds2)),by=list(dds2$group),min)[,2],
  "Median"=aggregate(colSums(assay(dds2)),by=list(dds2$group),median)[,2],
  "Max"=aggregate(colSums(assay(dds2)),by=list(dds2$group),max)[,2],
  "Mean"=aggregate(colSums(assay(dds2)),by=list(dds2$group),mean)[,2])

row.names(depthBygroup) <- levels(factor(dds2$group))
depthBygroup <- depthBygroup / 1E6

kable(depthBygroup)
```

### Assess sequencing effort with a boxplot

```{r sequenceEffortBoxplot}
boxplot(colSums(assay(dds2))~factor(dds2$group),col=palette_All_Data,ylab="Reads per sample (Millions)")
stripchart(colSums(assay(dds2))~factor(dds2$group),method="jitter",vertical=TRUE,add=TRUE,col="black")
```

## Block 5: Make sample distance heatmaps


```{r}
# ---- Pretty sample distance heatmap with rlog ----
sampleDistanceHeatmap <- function(
    rld,                               # DESeq2 rlog object
    method = "euclidean", 
    myTitle = "Sample Distance Heatmap", 
    figureOutputDirectory = NULL
) {
  require(pheatmap)
  require(RColorBrewer)
  
  # Calculate distances
  sample_dists <- dist(t(assay(rld)), method = method)
  sample_dist_DF <- as.data.frame(as.matrix(sample_dists))
  
  # Pull annotation data (use group + Sex from colData)
  annotation_DF <- as.data.frame(colData(rld)[, c("group", "Sex")])
  rownames(annotation_DF) <- rownames(colData(rld))
  
  # Drop unused levels
  annotation_DF$group <- droplevels(factor(annotation_DF$group))
  annotation_DF$Sex   <- droplevels(factor(annotation_DF$Sex))
  
  # Define colors for annotations (your original palette)
  anno_colors <- list(
    Sex = c(
      Female = "#F08080",     # light pink
      Male   = "darkgreen"    # green
    ),
    group = c(
      Egr1KD_gRNA2_Female = "#6A0DAD",
      Egr1KD_gRNA3_Female = "#B266FF",
      Egr1KD_gRNA2_Male   = "#FF8C00",
      Egr1KD_gRNA3_Male   = "#FFB84D",
      NoTreatment_gRNA1_Male   = "#1E90FF",
      NoTreatment_gRNA1_Female = "#FF69B4"
    )
  )
  
  # Generate heatmap
  p <- pheatmap(
    sample_dist_DF,
    col = colorRampPalette(c("blue", "white", "red"))(20),
    border_color = "black",
    fontsize = 10,
    fontsize_col = 7,
    show_rownames = FALSE,
    show_colnames = TRUE,
    annotation_col = annotation_DF,
    annotation_colors = anno_colors,
    clustering_distance_rows = sample_dists,
    clustering_distance_cols = sample_dists,
    main = myTitle
  )
  
  # Optionally save to PDF
  if (!is.null(figureOutputDirectory)) {
    pdf(file = file.path(
      figureOutputDirectory, 
      paste0(gsub(" ", "_", myTitle), ".pdf")
    ), width = 10, height = 8)
    p
    invisible(dev.off())
  }
  
  return(p)
}

# ---- Run it ----
pdf(file = file.path(figureOutputDirectory, "MandF_filteredOutliers_vstTransformed.pdf"), width = 7, height = 7)

p1 <- sampleDistanceHeatmap(
  rld, 
  method = "euclidean",
  myTitle = "Euclidean Distance on vst transformed counts",
  figureOutputDirectory = figureOutputDirectory
)
print(p1)
dev.off()
p1
```

Make a sample-distance heatmap using rlog transformed values

```{r}
#sampleDistanceHeatmap(counts = assay(rld), method = "euclidean",myTitle = "Euclidean Distance on Rlog transformed counts")
```

## Block 7: MDS Plots with Gliima

I use the package "Glimma" a lot. It creates searchable, interactive web pages to visualize differential expression results. It is a great tool to explore the data and is very easy to share with collaborators.

You can read the package vignette if you choose to work with this on your own.

I don't always make interactive PCA plots like this, but I'm including this example to show you the power of the package.

```{r}
# If you have never installed it before, first install the package
library(Glimma)

# interactive MDS plot
# don't inlude "sizeFactor" or "Barcode"
keepColumns <- ("Barcode" != colnames(colData(rld))) & ("sizeFactor" != colnames(colData(rld)))
group.df <- data.frame(colData(rld)[,keepColumns])

# Glimma is odd about folders
setwd(analysisOutputDirectory)

glMDSPlot(as.matrix(assay(rld)), groups = group.df, labels = row.names(group.df), 
          main = paste0("MDS of Rlog counts"),
          folder = "Glimma_MandF_filtered_MDS", launch=FALSE)


# make the plot with ALI samples only
#rld.lot6only <- rld[,colData(rld)$Lot=="L6"]
#keepColumns <- ("Barcode" != colnames(colData(rld.lot6only))) & ("sizeFactor" != colnames(colData(rld.lot6only)))
#group.df <- data.frame(colData(rld.lot6only)[,keepColumns])

#glMDSPlot(as.matrix(assay(rld.lot6only)), groups = group.df, labels = #row.names(group.df), 
#main = paste0("MDS of Rlog counts, lot6 only"),
#folder = "Glimma_MDS_Lot6_FEMALE_allDoses_Only_Demo", launch = #FALSE)

setwd(myWorkingDirectory)
```

## Block 8: Differential Expression Testing

Since we already put design into our dds object, we need just one line of code to perform DE testing

```{r}
dds2 <- DESeq(dds2)
```

```{r}
res <- results(dds2,contrast = c("group", "NoTreatment_gRNA1_Male","NoTreatment_gRNA1_Female"))
res$SYMBOL <- mcols(dds2)[,colnames(mcols(dds2)) %in% "external_gene_name"]
#res$ENTREZID <- mcols(dds2)[,colnames(mcols(dds2)) %in% "ENTREZID"]
#res$GENENAME <- mcols(dds2)[,colnames(mcols(dds2)) %in% "GENENAME"]

res[1:5,]

write.table(res,file=file.path(deOutputDirectory, "Male_NoTreatg1_vs_Female_NoTreatg1_DE_vst_filtered_091625.txt"),quote = FALSE,sep="\t",col.names = NA)
```


```{r}
res <- results(dds2,contrast = c("group", "Egr1KD_gRNA2_Male","NoTreatment_gRNA1_Male"))
res$SYMBOL <- mcols(dds2)[,colnames(mcols(dds2)) %in% "external_gene_name"]
#res$ENTREZID <- mcols(dds2)[,colnames(mcols(dds2)) %in% "ENTREZID"]
#res$GENENAME <- mcols(dds2)[,colnames(mcols(dds2)) %in% "GENENAME"]

res[1:5,]

write.table(res,file=file.path(deOutputDirectory, "Male_Egr1KDg2_vs_Male_NoTreatg1_DE_vst_filtered_091625.txt"),quote = FALSE,sep="\t",col.names = NA)
```

```{r}
res <- results(dds2,contrast = c("group", "Egr1KD_gRNA3_Male","NoTreatment_gRNA1_Male"))
res$SYMBOL <- mcols(dds2)[,colnames(mcols(dds2)) %in% "external_gene_name"]
#res$ENTREZID <- mcols(dds2)[,colnames(mcols(dds2)) %in% "ENTREZID"]
#res$GENENAME <- mcols(dds2)[,colnames(mcols(dds2)) %in% "GENENAME"]

res[1:5,]

write.table(res,file=file.path(deOutputDirectory, "Male_Egr1KDg3_vs_Male_NoTreatg1_DE_vst_filtered_091625.txt"),quote = FALSE,sep="\t",col.names = NA)
```

```{r}
res <- results(dds2,contrast = c("group", "Egr1KD_gRNA2_Female","NoTreatment_gRNA1_Female"))
res$SYMBOL <- mcols(dds2)[,colnames(mcols(dds2)) %in% "external_gene_name"]
#res$ENTREZID <- mcols(dds2)[,colnames(mcols(dds2)) %in% "ENTREZID"]
#res$GENENAME <- mcols(dds2)[,colnames(mcols(dds2)) %in% "GENENAME"]

res[1:5,]

write.table(res,file=file.path(deOutputDirectory, "Female_Egr1KDg2_vs_Female_NoTreatg1_DE_vst_filtered_091625.txt"),quote = FALSE,sep="\t",col.names = NA)
```
```{r}
res <- results(dds2,contrast = c("group", "Egr1KD_gRNA3_Female","NoTreatment_gRNA1_Female"))
res$SYMBOL <- mcols(dds2)[,colnames(mcols(dds2)) %in% "external_gene_name"]
#res$ENTREZID <- mcols(dds2)[,colnames(mcols(dds2)) %in% "ENTREZID"]
#res$GENENAME <- mcols(dds2)[,colnames(mcols(dds2)) %in% "GENENAME"]

res[1:5,]

write.table(res,file=file.path(deOutputDirectory, "Female_Egr1KDg3_vs_Female_NoTreatg1_DE_vst_filtered_091625.txt"),quote = FALSE,sep="\t",col.names = NA)
```

```{r}
res <- results(dds2,contrast = c("group", "Egr1KD_gRNA2_Male","Egr1KD_gRNA2_Female"))
res$SYMBOL <- mcols(dds2)[,colnames(mcols(dds2)) %in% "external_gene_name"]
#res$ENTREZID <- mcols(dds2)[,colnames(mcols(dds2)) %in% "ENTREZID"]
#res$GENENAME <- mcols(dds2)[,colnames(mcols(dds2)) %in% "GENENAME"]

res[1:5,]

write.table(res,file=file.path(deOutputDirectory, "Male_Egr1KDg2_vs_Male_Egr1KDg2_DE_vst_filtered_091625.txt"),quote = FALSE,sep="\t",col.names = NA)
```
```{r}
res <- results(dds2,contrast = c("group", "Egr1KD_gRNA3_Male","Egr1KD_gRNA3_Female"))
res$SYMBOL <- mcols(dds2)[,colnames(mcols(dds2)) %in% "external_gene_name"]
#res$ENTREZID <- mcols(dds2)[,colnames(mcols(dds2)) %in% "ENTREZID"]
#res$GENENAME <- mcols(dds2)[,colnames(mcols(dds2)) %in% "GENENAME"]

res[1:5,]

write.table(res,file=file.path(deOutputDirectory, "Male_Egr1KDg3_vs_Female_Egr1KDg3_DE_vst_filtered_091625.txt"),quote = FALSE,sep="\t",col.names = NA)
```