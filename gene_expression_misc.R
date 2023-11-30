source('code/gene_expression_helpers.R')

# ---- summarize gene expression tests performed ----

# summarize association tests for a chromosome
# chrom: name of chromosome in filename
# returns: data.frame with summary of tests for this chromosome
summarize_tests <- function(chrom) {
  # count tests by tissue
  tests <- read.csv(test_filename(chrom))
  summary <- aggregate(ensembl_id ~ tissue + n1 + n2, tests, length)
  
  # rename columns to add chromosome name, for disambiguation
  colnames(summary) <- c('tissue', paste0('n ', chrom, c(1:2, ' tests')))
  summary
}

# put Y and MT side-by-side
summary <- merge.data.frame(summarize_tests('Y'), 
                            summarize_tests('MT'), by = 'tissue')
# add full name of tissue for context
summary$tissue <- paste0(tissue_map[summary$tissue], ' (', summary$tissue, ')')
# table S2
write.csv(summary, 'results/associations/gene_expression_summary.csv', 
          quote = FALSE, row.names = FALSE)

# ---- extracting significant gene expression associations ----

# extract significant associations and save to file
# chrom: name of chromosome in filename
# returns: nothing
extract_sig_tests <- function(chrom) {
  # read sig tests
  tests <- read.csv(test_filename(chrom))
  tests <- tests[tests$adj_p < 0.05,]
  
  # look up gene and tissue name
  tests$gene <- gene_map[tests$ensembl_id]
  tests$tissue <- paste0(tissue_map[tests$tissue], ' (', tests$tissue, ')')
  
  # subset to only necessary columns
  tests <- tests[, c('ensembl_id', 'gene', 'tissue', 'chr', 'adj_p')]
  colnames(tests) <- c('Ensembl ID', 'gene', 'tissue', 'chr', 'q-value')
  
  write.csv(tests, 
            paste0('results/associations/gene_expression_', 
                   chrom, '_sig_tests.csv'), 
            quote = FALSE, row.names = FALSE)
}

# table 1
extract_sig_tests('Y')
# table 2
extract_sig_tests('MT')

# ---- demonstrating why IQN fails for highly-expressed genes ----

# load TPM and IQN files for the same tissue
tpm <- load_rsem('Brain', 'tpm')[, -c(1:3)]
iqn <- load_rsem('Brain', 'iqn')[, -c(1:3)]
# an MT gene which is highly expressed in many samples
nd4 <- 'ENSRNOG00000029707'

# rank() assigns small numbers small ranks; negate so large TPM are ranked first
ranked_tpm <- -tpm
ranked_tpm[] <- lapply(ranked_tpm, rank)

# 296 samples have Mt-nd4 has the 9th-most TPM
table(unlist(ranked_tpm[nd4,]))
# Therefore, 296 samples have identical IQN
table(unlist(iqn[nd4,]))