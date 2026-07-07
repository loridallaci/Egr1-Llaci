library(limma)
library(ggplot2)
library(dplyr)

#############
# calculate_sequence_effort_by_group_table
#############
#
# Usage
# calculate_sequence_effort_by_group_table(myDGEList=y,myGroup=y$samples$Group)
#
# Arguments
#   myDGEList: a DGE list object returned by call to Limma::DGEList
#   myGroup: a vector, usually from the samples table, with the sample grouping to compare
#
# Return
#   a data frame with minimum sequencing depth, maximum sequencing depth, and 
#    median sequencing depth for the sample groups (in millions of counts)

calculate_sequence_effort_by_group_table <- function(
    myDGEList,
    myGroup
)
{
  depth_by_group <- data.frame("Minimum_Depth"=aggregate(colSums(myDGEList$counts),by=list(myGroup),min)[,2], 
                             "Maximum_Depth"=aggregate(colSums(myDGEList$counts),by=list(myGroup),max)[,2],
                             "Median_Depth"=aggregate(colSums(myDGEList$counts),by=list(myGroup),median)[,2],
                             "Mean_Depth"=aggregate(colSums(myDGEList$counts),by=list(myGroup),mean)[,2])
  # rename with the experimental groups
  row.names(depth_by_group) <- levels(factor(myGroup))
  
  # express as millions of counts
  depth_by_group <- round(x=depth_by_group / 1E6,digits = 2)
}

#############
# stripchart_sequence_effort_by_group
#############
#
# Usage
# stripchart_sequence_effort_by_group(myDGEList = y,myGroup = y$samples$Group)
#
# Arguments
#   myDGEList: a DGE list object returned by call to Limma::DGEList
#   myGroup: a vector, usually from the samples table, with the sample grouping to compare
#
# Return
#   a ggplot2 object. The plot is a stripchart (boxplot) of counts per sample,
#     aggregated by metadata group. The plot includes an aesthetic "text" which
#     can be used to generate an interactive plot with plotly functions
#
# Example: 
#  stripchart <- stripchart_sequence_effort_by_group(y,y$samples$Group)
#  stripchart %>% ggplotly(tooltip = c("text","y"))

stripchart_sequence_effort_by_group <- function(
    myDGEList,
    myGroupFactor
){
  data <- data.frame("Counts"=colSums(myDGEList$counts),
                     "Group"=myGroupFactor,
                     "Sample"=colnames(myDGEList))
  
  data <- data %>% mutate(MillionCounts = Counts / 1E6)
  
  stripchart <- ggplot(data, aes(x = Group, y = MillionCounts, fill = Group, text = Sample)) +
    geom_boxplot() +
    geom_jitter() +
    theme_bw() +
    ggtitle("Sequencing Effort") +
    ylab("Counts (Millions)")
  
  return(stripchart)
}
