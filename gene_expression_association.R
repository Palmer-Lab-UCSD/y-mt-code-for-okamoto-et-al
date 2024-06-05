# ----setup ----

source('code/gene_expression_helpers.R')

# run Wilcoxon rank-sums test on TMM-normalized RSEM
# rsem: data.frame with genes x samples raw RSEM values; genes are row names
# chrom: name of chromosome (Y or MT)
# returns: data.frame with gene, p-value, and nonzero counts by group
tmm_wilcox <- function(rsem, chrom) {
  rsem <- filter_rsem(rsem, chrom)
  rfids <- colnames(rsem)
  is_haplo1 <- startsWith(chrom_groups(chrom, rfids), paste0(chrom, '1'))
  
  cat('Converting filtered', chrom, 'RSEM to TMM-normalized CPM', '\n')
  tmm <- rsem_to_cpm(rsem)
  
  pvals <- apply(tmm, 1, function(x)
    wilcox.test(x[is_haplo1], x[!is_haplo1], exact = FALSE)$p.value)
  
  data.frame(test_chr = chrom, ensembl_id = rownames(tmm), 
             n1 = sum(is_haplo1), n2 = sum(!is_haplo1), p = pvals)
}

# run gene expression associations against nonrecombinant chromosome haplotypes
# tissue: short name of the tissue on RatGTex
# returns: data.frame with results of all tests, against both Y and MT
haplogroup_tests <- function(tissue) {
  rsem <- load_rsem(tissue)
  # separate gene metadata from expression values
  metadata <- rsem[ , c(1:3)]
  rsem <- rsem[ , -c(1:3)]
  
  y_pvals <- tmm_wilcox(rsem, 'Y')
  mt_pvals <- tmm_wilcox(rsem, 'MT')
  
  all_pvals <- rbind(merge.data.frame(metadata, y_pvals, by = 'ensembl_id'),
                     merge.data.frame(metadata, mt_pvals, by = 'ensembl_id'))
  all_pvals$tissue <- tissue
  
  all_pvals
}

# ---- run analysis -----

tests <- do.call('rbind', lapply(names(tissue_map), haplogroup_tests))

for (chrom in c('Y', 'MT')) {
  chrom_tests <- tests[tests$test_chr == chrom, 
                       c('ensembl_id', 'chr', 'tissue', 'n1', 'n2', 'p')]
  chrom_tests$adj_p <- p.adjust(chrom_tests$p, method = 'BH')
  chrom_tests <- chrom_tests[order(chrom_tests$p), ]
  
  write.csv(chrom_tests, test_filename(chrom), quote = FALSE, row.names = FALSE)

  sig_tests <- chrom_tests[chrom_tests$adj_p < ALPHA, ]
  sig_tests$gene <- gene_map[sig_tests$ensembl_id]
  long_tissue_names <- tissue_map[sig_tests$tissue]
  sig_tests$tissue <- paste0(long_tissue_names, ' (', sig_tests$tissue, ')')
  
  sig_tests <- sig_tests[ , c('ensembl_id', 'gene', 'tissue', 'chr', 'adj_p')]
  colnames(sig_tests) <- c('Ensembl ID', 'gene', 'tissue', 'chr', 'q-value')
  
  write.csv(sig_tests, gsub('tests.csv', 'tests_sig.csv', test_filename(chrom)), 
            quote = FALSE, row.names = FALSE)
}
