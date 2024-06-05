# ---- load requirements ----

library(edgeR) # TMM (calcNormFactors, cpm)

ALPHA <- 0.05

gene_map <- read.csv('data/expression/gene_name_map.csv')
gene_map <- setNames(gene_map$gene_name, nm = gene_map$ensembl_id)

tissue_map <- read.csv('data/expression/tissue_name_map.csv')
tissue_map <- setNames(tissue_map$long_name, nm = tissue_map$short_name)

y_groups <- read.csv('results/groups/Y_groups.csv', row.names = 'rfid')
mt_groups <- read.csv('results/groups/MT_groups.csv', row.names = 'rfid')

# ---- helper functions ----

# load haplotype groups for given rats
# chrom: name of chromosome to load haplotypes for
# rfids: list of RFIDs to read haplotypes groups for
# returns: list of haplotype groups in order of RFID
chrom_groups <- function(chrom, rfids) {
  if (chrom == 'Y') { 
    y_groups[rfids, 'Y_group'] 
  } else { 
    mt_groups[rfids, 'MT_group']
  }
}

# load RSEM from RatGTex
# tissue: short name of the tissue on RatGTex
# type: what type of expression data (log2, tpm, iqn) to load; default 'log2'
# returns: data.frame with genes x samples raw RSEM values; genes are row names
#          first three columns are metadata instead of samples
load_rsem <- function(tissue, type = 'log2') {
  cat('Loading', tissue, 'data', '\n')
  rsem <- read.csv(paste0('data/expression/', tissue, '.rn7.expr.', 
                          type, '.bed.gz'),
                   # check.names = FALSE keeps #chr from becoming X.chr, etc.
                   sep = '\t', check.names = FALSE)
  
  metadata <- data.frame(rsem[ , c(1, 4)])
  colnames(metadata) <- c('chr', 'ensembl_id')
  metadata$tissue <- tissue
  
  # data.table can't handle row names
  rsem <- data.frame(rsem[ , -c(1:3)], row.names = 'gene_id', check.names = F)
  if (type == 'log2') rsem <- round((2 ** rsem) - 1, 2)
  
  cbind(metadata, rsem)
}

# filter RSEM samples and genes for low expression or lack of haplotype
# rsem: data.frame with genes x samples raw RSEM values; genes are row names
# chrom: name of chromosome for logging purposes
# returns: genes x samples raw RSEM values with some rows and columns removed
filter_rsem <- function(rsem, chrom) {
  no_haplotype <- is.na(chrom_groups(chrom, colnames(rsem)))
  cat(sum(no_haplotype), 'rats removed for lacking a', chrom, 'haplotype', '\n')
  rsem <- rsem[ , !no_haplotype]
  
  low_expr_genes <- rowSums(rsem > 0) < ncol(rsem) * 0.1
  cat(sum(low_expr_genes), 'genes removed for low expression.',
      '<10% of samples with a', chrom, 'haplotype have RSEM > 0', '\n')
  rsem <- rsem[!low_expr_genes, ]
  
  cat(nrow(rsem), 'genes and', ncol(rsem), 
      'samples for', chrom, 'haplotype tests', '\n')
  rsem
}

# convert raw RSEM values to TMM-normalized CPM
# rsem: data.frame with genes x samples raw RSEM values; genes are row names
# returns: data.frame of the same size but with values in CPM instead
rsem_to_cpm <- function(rsem) {
  dgelist <- DGEList(rsem)
  dgelist <- calcNormFactors(dgelist, method = 'TMM')
  cpm(dgelist)
}

# convert a chromosome name to a gene expression result filename
# chrom: name of chromosome in filename
# returns: filename as a string
test_filename <- function(chrom)
  paste0('results/associations/', chrom, '_gene_expression_tests.csv')