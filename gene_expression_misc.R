source('code/gene_expression_helpers.R')

# ---- summarize gene expression tests performed ----

# summarize association tests for a chromosome
# chrom: name of chromosome in filename
# returns: data.frame with summary of tests for this chromosome
count_tests <- function(chrom) {
  tests <- read.csv(test_filename(chrom))
  # use long name in summary table, for clarity
  tests$tissue <- paste0(tissue_map[tests$tissue], ' (', tests$tissue, ')')
  # count tests performed
  summary <- aggregate(ensembl_id ~ tissue + n1 + n2, tests, length)

  # disambiguate columns by adding chromosome name
  colnames(summary) <- c('tissue', paste0('n ', chrom, c(1:2, ' tests')))
  summary
}

summary <- merge.data.frame(count_tests('Y'), count_tests('MT'), by = 'tissue')
# table S2
write.csv(summary, 'results/associations/gene_expression_summary.csv', 
          quote = FALSE, row.names = FALSE)

# ---- saving TMM-normalized CPM for certain gene-tissue combos ----

# save TMM for one gene-tissue combination to file
# x: short tissue name, then Ensembl ID
# returns: nothing
save_tmm <- function(x) {
  # preprocess RSEM as done during analysis (but always keep females)
  rsem <- filter_rsem(load_rsem(x[1])[ , -c(1:3)], 'MT')
  tmm <- rsem_to_cpm(rsem)
  
  expr <- data.frame(rfid = colnames(tmm), expr = tmm[x[2], ])
  filename <- paste0('data/expression/', x[2], '_', x[1], '_expr.csv')
  write.csv(expr, filename, row.names = FALSE)
}

y_tests <- read.csv(test_filename('Y'))
apply(y_tests[y_tests$adj_p < ALPHA, c('tissue', 'ensembl_id')], 1, save_tmm)

mt_tests <- read.csv(test_filename('MT'))
apply(mt_tests[mt_tests$adj_p < ALPHA, c('tissue', 'ensembl_id')], 1, save_tmm)

# genes around/in a PAV on Y; see File S1
y_pav_genes <- names(gene_map)[gene_map %in% c('Dkc1', 'Med14Y')]
y_pav_to_save <- expand.grid(names(tissue_map), y_pav_genes)
apply(y_pav_to_save, 1, save_tmm)

# ---- saving TMM-normalized CPM for PCA on top Y/MT genes ----

cis_contigs <- read.table('data/expression/Y_MT_cis_contigs.tsv', header = 1)
brain_rsem <- load_rsem('Brain')
# chromosome is the only useful metadata from the first 3 columns
gene_chrs <- setNames(brain_rsem$chr, nm = rownames(brain_rsem))
brain_rsem <- brain_rsem[ , -c(1:3)]

for (chrom in c('Y', 'MT')) {
  tmm <- rsem_to_cpm(filter_rsem(brain_rsem, chrom))
  gene_ids <- rownames(tmm)
  contigs <- strsplit(cis_contigs[cis_contigs$chr == chrom, 'contigs'], ',')

  chrom_genes <- tmm[gene_chrs[gene_ids] %in% contigs[[1]], ]
  half_have_expr <- apply(chrom_genes, 1, function(x) mean(x != 0) >= 0.5 )
  chrom_genes <- t(chrom_genes[half_have_expr, ])

  write.csv(chrom_genes, paste0('data/expression/brain_', chrom, '_genes.csv'), 
            quote = FALSE)
}

# ---- demonstrating why IQN fails for highly-expressed genes ----

tpm <- load_rsem('Brain', 'tpm')[ , -c(1:3)]
iqn <- load_rsem('Brain', 'iqn')[ , -c(1:3)]
nd4 <- names(gene_map)[gene_map == 'Mt-nd4']

# rank() assigns small numbers small ranks; negate so large TPM are ranked first
ranked_tpm <- -tpm
ranked_tpm[] <- lapply(ranked_tpm, rank)

# 296 samples have Mt-nd4 as the 9th-most TPM
table(unlist(ranked_tpm[nd4, ]))
# Therefore, 296 samples have identical IQN
table(unlist(iqn[nd4, ]))
