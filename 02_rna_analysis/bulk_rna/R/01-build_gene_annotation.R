library(biomaRt)

# Usage
# gene_table <- build_gene_table_from_ensembl(ensembl_ids=row.names(count_data), species="mouse")
# 
# Arguments
#   ensemble_ids: a vector containing a list of Ensembl ID's. Intended use is to 
#       pass the row.names used to build the current count table
#   species: "mouse" or "human", to choose the appropriate mart
#
# Return
#   A data.frame, with a row for every element of the ensembl_id vector. 
#   For every row, there are 5 columns: "GENENAME","ENTREZID","BIOTYPE","SYMBOL","ENSEMBLID")
#   If no value is available, it returns NA

build_gene_table_from_ensembl <- function(
  ensembl_ids,
  mySpecies="mouse" # or "human"
){
  # Building on code shared by Kevin Blighe, https://www.biostars.org/p/388949/#388960
  
  mart <- useMart('ENSEMBL_MART_ENSEMBL')
  
  if(mySpecies == "mouse"){
    mart <- useDataset('mmusculus_gene_ensembl', mart)
  }
  else(
    if(mySpecies == "human"){
      mart <- useDataset('hsapiens_gene_ensembl', mart)
    }
    else{
      return("Species must be 'mouse' or 'human' in function build_gene_table_from_ensemble\n\n")
    }
  )
  
  # searchAttributes(mart,pattern = "symbol")
  annot_lookup <- getBM(
    mart = mart,
    attributes = c(
      'wikigene_description',
      'ensembl_gene_id',
      'entrezgene_id',
      'gene_biotype',
      #    'mgi_symbol',
      'external_gene_name'))
  
  # Use the ensembl keys to build the annotation
  annot_lookup <- annot_lookup[which(annot_lookup$ensembl_gene_id %in% ensembl_ids),]
  annot_lookup <- annot_lookup[match(ensembl_ids, annot_lookup$ensembl_gene_id),]
  all.equal(ensembl_ids, annot_lookup$ensembl_gene_id) # check that annots are aligned
  
  gene_data <- data.frame(matrix(nrow =length(ensembl_ids),ncol=5))
  row.names(gene_data) <- ensembl_ids
  colnames(gene_data) <- c("GENENAME","ENTREZID","BIOTYPE","SYMBOL","ENSEMBLID")
  gene_data$GENENAME <- annot_lookup$wikigene_description
  gene_data$ENTREZID <- annot_lookup$entrezgene_id
  gene_data$BIOTYPE <- annot_lookup$gene_biotype
  gene_data$SYMBOL <- annot_lookup$external_gene_name
  gene_data$ENSEMBLID <- ensembl_ids
  
  # replace blanks with NA, and <NA> with "NA"
  gene_data[gene_data==""] <- "NA"
  gene_data[is.na(gene_data)] <- "NA"
  
  return(gene_data)
}