# ---- load requirements ----

source('code/gene_expression_helpers.R')

# remove notes about what founders this chromosome is from, e.g. Y1 (ACI) -> Y1
y_groups$Y_group <- substr(unlist(y_groups$Y_group), start = 1, stop = 2)
mt_groups$MT_group <- substr(unlist(mt_groups$MT_group), start = 1, stop = 3)

groups <- list('Y' = y_groups, 'MT' = mt_groups)

# ---- helper functions ----

# run Wilcoxon rank-sums test on TMM-normalized RSEM
# rsem: data.frame with genes x samples raw RSEM values; genes are row names
# chrom: name of chromosome in haplotype labels (e.g. Y from Y1)
# returns: data.frame with gene, p-value, and nonzero counts by group
tmm_wilcox <- function(rsem, chrom) {
  # preprocess to remove bad genes & samples (see filter_rsem for details)
  rsem <- filter_rsem(rsem, chrom)
  rfids <- colnames(rsem)
  # split by haplotype
  is_group1 <- groups[[chrom]][rfids, 1] == paste0(chrom, '1')
  
  # TMM-normalized CPM is named "tmm" to avoid naming conflicts
  tmm <- rsem_to_cpm(rsem)
  cat('Converting filtered', chrom, 'RSEM to TMM-normalized CPM', '\n')
  
  # run association tests
  pvals <- apply(tmm, 1, function(x)
    wilcox.test(x[is_group1], x[!is_group1], exact = FALSE)$p.value)
  
  data.frame(test_chr = chrom, ensembl_id = rownames(tmm), 
             n1 = sum(is_group1), n2 = sum(!is_group1), p = pvals)
}

# run gene expression associations against nonrecombinant chromosome haplotypes
# tissue: short name of the tissue on RatGTex
# returns: data.frame with results of all tests, against both Y and M
haplogroup_tests <- function(tissue) {
  rsem <- load_rsem(tissue)
  # separate gene metadata from expression values
  metadata <- rsem[, c(1:3)]
  rsem <- rsem[, -c(1:3)]
  
  # run separately for each chromosome
  y_pvals <- tmm_wilcox(rsem, 'Y')
  mt_pvals <- tmm_wilcox(rsem, 'MT')
  
  # add metadata back in
  rbind(merge.data.frame(metadata, y_pvals, by = 'ensembl_id'),
        merge.data.frame(metadata, mt_pvals, by = 'ensembl_id'))
}

# ---- run analysis -----

tests <- do.call('rbind', lapply(names(tissue_map), haplogroup_tests))

# save relevant result columns for a specific chromosome
# chrom: name of chromosome for extracting rows
# returns: nothing
save_chrom_results <- function(chrom) {
  chrom_tests <- tests[tests$test_chr == chrom, 
                       c('ensembl_id', 'chr', 'tissue', 'n1', 'n2', 'p')]
  
  chrom_tests$adj_p <- p.adjust(chrom_tests$p, method = 'BH')
  
  write.csv(chrom_tests[order(chrom_tests$p), ], test_filename(chrom), 
            quote = FALSE, row.names = FALSE)
}

save_chrom_results('Y')
save_chrom_results('MT')