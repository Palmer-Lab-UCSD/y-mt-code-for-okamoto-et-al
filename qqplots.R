# ---- load requirements ----

library(ggplot2) # QQ plot
library(cowplot) # QQ plot theme and arrangement

source('code/gene_expression_helpers.R')

# base colors for points
BLACK <- '#000000'
GREY <- '#888888'

# https://www.nature.com/articles/nmeth.1618
COLORBLIND_PALETTE <- c(GREY, BLACK, '#E69F00', '#56B4E9', '#009E73', '#F0E442')

# ---- core QQ plot functions ----

# match up p-values (log10-transformed) to expected quantiles, large to small
# pvals: p-values from tests, possibly with NAs
# other_cols: metadata corresponding to each p-value, to be retained in output
# returns: data.frame with ordered p-values and corresponding expected quantiles
#          columns are p (-log10 transformed), quantile, is_sig, then other_cols
expected_pvals <- function(pvals, other_cols) {
  # remove tests lacking a p-value
  na_pval <- is.na(pvals)
  other_cols <- other_cols[!na_pval,]
  pvals <- pvals[!na_pval]
  
  # -log10 of a (null) uniform distribution is exponential with mean log(10)
  quantiles <- (1:length(pvals) - 0.5) / length(pvals)
  quantiles <- qexp(quantiles, rate = log(10))
  
  # p < alpha = 0.05 AFTER adjustment is significant
  is_sig <- p.adjust(pvals, method = 'BH') < 0.05
  
  # p-values and their metadata must be sorted to match quantiles
  sorted_order <- order(pvals, decreasing = TRUE)
  # use -log10 transformed p-values to emphasize significant ones
  cbind(data.frame(p = -log10(pvals[sorted_order]), quantile = quantiles,
                   is_sig = is_sig[sorted_order]),
        other_cols[sorted_order,])
}

# make (and save) a QQ plot for phenotype associations
# chrom: name of chromosome for plot title
# description: description of test set for plot title
# color_priority: priority of colors to plot over each other; increasing order
# returns: completed QQ plot, all points plotted and no shape legend
pval_qqplot <- function(tests, chrom, description, color_priority) {
  # add expected quantiles to data in order to plot
  to_plot <- expected_pvals(tests$p, tests[, c('color', 'text')])
  
  plt <- ggplot(to_plot, aes(quantile, p, color = color, shape = is_sig)) + 
    # plot actual vs. expected and line of expectation
    geom_point() + geom_abline(slope = 1) + guides(shape = 'none') +
    
    geom_text(aes(label = text), color = BLACK, parse = TRUE,
              vjust = 'middle', hjust = 'top', nudge_x = -0.05) +
    
    # increase font size for readability
    theme(title = element_text(size = 15)) + theme_cowplot() +
    
    labs(title = paste0('QQ plot for ', chrom, ' haplotype to\n',
                        description, ' associations'),
         x = 'Expected -log10(pval)', y = 'Actual -log10(pval)')
  
  if (length(color_priority) > 1) sapply(color_priority, function(color) {
    # re-plot dots with color priority so they appear on top
    plt <<- plt + geom_point(data = to_plot[to_plot$color == color,])
  })
  
  plt
}

# ---- specialized QQ plot functions -----

# make a QQ plot for tests run on old GWAS traits with GCTA
# chrom: name of chromosome for plot title and filename
# returns: completed QQ plot, customized for these traits
gcta_qqplot <- function(chrom) {
  tests <- read.csv(paste0('results/associations/gwas_phenotype_', 
                           chrom, '_tests.csv'))
  # label any and all traits which are significant
  tests$text <- ifelse(tests$adj_p < 0.05, tests$name, '')
  tests$color <- ''
  
  # disable color legend since all dots are the same color
  pval_qqplot(tests, chrom, 'GWAS phenotype', c('')) + guides(color = 'none')
}

# make a QQ plot for tests run on gene expression data
# chrom_label: name of chromosome for plot title and filename
# cis_chroms: names of chromosomes to consider as cis to chrom_label
# genes_ids: Ensembl IDs of genes to label (i.e. the top associations)
# returns: completed QQ plot, customized for these traits
eqtl_qqplot <- function(chrom_label, cis_chroms, gene_ids) {
  tests <- read.csv(test_filename(chrom_label))
  # convert tissue names from e.g. BLA -> Basolateral~amygdala for labels
  tests$tissue <- gsub(' ', '~', tissue_map[unlist(tests$tissue)])

  # coloring will be done by gene
  tests$color <- ifelse(
    # convert Ensembl IDs to common names for the top genes
    tests$ensembl_id %in% gene_ids, gene_map[tests$ensembl_id],
    # non-top genes are colored by cis-ness
    ifelse(tests$chr %in% cis_chroms, 'other (cis)', 'other (trans)')
  )
  # need to be a factor to use as a color aesthetic
  color_order <- c('other (trans)', 'other (cis)', gene_map[gene_ids])
  tests$color <- factor(tests$color, levels = color_order)
  
  # label only the top genes (identified by p-value)
  best_p <- head(sort(tests$p), length(gene_ids))
  tests$text <- ifelse(tests$p %in% best_p,
                       paste0('italic("', tests$color, '")~"in"~', 
                              tests$tissue), '')
  
  pval_qqplot(tests, chrom_label, 'gene expression', color_order) +
    # force colors to be from colorblind palette
    scale_color_manual(limits = color_order, values = COLORBLIND_PALETTE) +
    # legend will be at top and will not have a label for itself
    theme(legend.position = 'top', 
          legend.text = element_text(size = 15, face = 'italic')) +
    labs(color = '') + guides(color = guide_legend(nrow = 2))
}

# ---- combined QQ plots ----

# make side-by-side QQ plots for both kinds of tests
# chrom: name of chromosome for plot titles and filenames
# cis_chroms: names of chromosomes to consider as cis to chrom_label
# gene_names: common names of the top genes
# returns: finished plot grid
combined_qqplot <- function(chrom, cis_chroms, gene_names) {
  gene_ids <- sapply(gene_names, function(x) names(gene_map)[gene_map == x])
  
  gcta_plot <- gcta_qqplot(chrom)
  eqtl_plot <- eqtl_qqplot(chrom, cis_chroms, unlist(gene_ids))
  
  # place plots side-by-side and label the panels
  plot_grid(gcta_plot, eqtl_plot, labels = 'AUTO', ncol = 2)
}