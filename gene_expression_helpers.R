# ---- load requirements ----

library(edgeR) # TMM (calcNormFactors, cpm)

# map from short to long name for tissues
tissue_map <- c(
  'BLA' = 'Basolateral amygdala', 'Brain' = 'Brain hemisphere', 
  'Eye' = 'Eye', 'IL' = 'Infralimbic cortex', 'LHb' = 'Lateral habenula', 
  'NAcc' = 'Nucleus accumbens core', 'NAcc2' = 'Nucleus accumbens core',
  'OFC' = 'Orbitofrontal cortex', 'PL' = 'Prelimbic cortex', 
  'PL2' = 'Prelimbic cortex'
)

# map from gene IDs to common names
gene_map <- c(
  'ENSRNOG00000057231' = 'Ddx3y', 'ENSRNOG00000055562' = 'Dkc1',
  'ENSRNOG00000033615' = 'Mt-nd3', 'ENSRNOG00000043866' = '16S rRNA', 
  'ENSRNOG00000029971' = 'Mt-nd5', 'ENSRNOG00000029042' = 'Mt-nd6',
  'ENSRNOG00000029707' = 'Mt-nd4', 'ENSRNOG00000030644' = 'Mt-nd1',
  'ENSRNOG00000030478' = '12S rRNA', 'ENSRNOG00000031033' = 'Mt-nd2',
  'ENSRNOG00000031053' = 'Mt-nd4l'
)

y_groups <- read.csv('results/groups/Y_groups.csv', row.names = 'rfid')
mt_groups <- read.csv('results/groups/MT_groups.csv', row.names = 'rfid')

# ---- helper functions ----

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
  
  # process metadata: chromosome and gene ID
  metadata <- data.frame(rsem[, c(1, 4)])
  colnames(metadata) <- c('chr', 'ensembl_id')
  metadata$tissue <- tissue
  
  # data.table can't handle row names
  rsem <- data.frame(rsem[, -c(1:3)], row.names = 'gene_id', check.names = F)
  # now all values are numeric, it is safe to convert log_2(RSEM + 1) -> RSEM
  if (type == 'log2') rsem <- round((2 ** rsem) - 1, 2)
  
  cbind(metadata, rsem)
}

# filter RSEM samples and genes for low expression or lack of haplotype
# rsem: data.frame with genes x samples raw RSEM values; genes are row names
# chrom: name of chromosome for logging purposes
# returns: genes x samples raw RSEM values with some rows and columns removed
filter_rsem <- function(rsem, chrom) {
  # subset to only rats with a haplotype before filtering on expression
  no_haplotype <- is.na(groups[[chrom]][colnames(rsem), 1])
  cat(sum(no_haplotype), 'rats removed for lacking a', chrom, 'haplotype', '\n')
  rsem <- rsem[, !no_haplotype]
  
  # filter out genes with little expression
  low_expr <- rowSums(rsem > 0) < ncol(rsem) * 0.1
  cat(sum(low_expr), 'genes removed for low expression.',
      '<10% of samples with a', chrom, 'haplotype have RSEM > 0', '\n')
  rsem <- rsem[!low_expr,]
  
  cat(nrow(rsem), 'genes and', ncol(rsem), 
      'samples for', chrom, 'haplotype tests', '\n')
  rsem
}

# convert raw RSEM values to TMM-normalized CPM
# rsem: data.frame with genes x samples raw RSEM values; genes are row names
# returns: data.frame of the same size but with values in CPM instead
rsem_to_cpm <- function(rsem) {
  # build edgeR-friendly object for normalization
  dgelist <- DGEList(rsem)
  dgelist <- calcNormFactors(dgelist, method = 'TMM')
  
  # export TMM-normalized CPM
  cpm(dgelist)
}

# convert a chromosome name to a gene expression result filename
# chr: name of chromosome in filename
# returns: filename as a string
test_filename <- function(chrom)
  paste0('results/associations/gene_expression_', chrom, '_tests.csv')