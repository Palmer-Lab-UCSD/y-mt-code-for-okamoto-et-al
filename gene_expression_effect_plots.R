# ---- environment set-up ----

# prevent scientific notation for axes
options(scipen = 3)

library(ggplot2) # effect plots
library(cowplot) # plot theme and arrangement

source('code/gene_expression_helpers.R')

# load sex information and add to y_groups (to make Y1, Y2, and F rats)
sample_info <- read.csv('data/sample_info.csv')
y_groups <- ifelse(sample_info$sex == 'F', 'Female', 
                   y_groups[sample_info$rfid, 1])
names(y_groups) <- sample_info$rfid

mt_groups$MT_group <- ifelse(mt_groups$MT_group == 'MT1 (BN)', 
                             mt_groups$MT_group, 'MT2 (BUF etc.)')

groups <- list('Y' = as.data.frame(y_groups), 'MT' = as.data.frame(mt_groups))

# ---- helper functions ----

# read significant gene expression associations from a file
# chrom: name of chromosome in filename
# returns: data.frame of significant hits
read_sig_hits <- function(chrom) {
  tests <- read.csv(test_filename(chrom))
  tests[tests$adj_p < 0.05, ]
}

# make an effect plot
# x: information about the significant hit
# tmm: data.frame with genes x samples TMM-normalized CPM; genes are row names
# groups_df: data.frame with RFIDs as row names and haplotype as first column
# chrom: name of chromosome for the title and p-value column names
# returns: completed effect plot as a ggplot object
effect_plot <- function(x, tmm, groups_df, chrom) {
  rfids <- colnames(tmm)
  
  # pull ID information (for use with gene/tissue_map)
  ensembl_id <- x[['ensembl_id']]
  tissue_id <- x[['tissue']]
  
  # gene expression values
  vals <- tmm[ensembl_id,]
  
  # look up haplotype for each sample, and count n for each group
  type <- groups_df[rfids, 1]
  type <- paste0(type, '\nn = ', table(type)[type], '')
  # make a factor for plotting purposes
  type <- factor(type)
  
  # plot gene expression on the y axis, separating haplotypes along the x
  ggplot(data.frame(val = vals, haplotype = type), aes(haplotype, val)) + 
    
    # violin plot with box-plot quantiles drawn across
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), scale = 'width') + 
    # jitter points only to the side, to keep them at the correct y position
    geom_jitter(height = 0, width = 0.2) + 
    
    # clean up visual appearance, then blow up tick labels for readability
    theme_minimal_hgrid() + theme(axis.text.x = element_text(size = 15)) +
    
    labs(
      title = bquote(.(tissue_map[tissue_id])~'expression of'
                     ~italic(.(gene_map[ensembl_id]))), x = '', y = 'CPM'
    )
}

# generate effect plots for all significant DE genes
# chrom: name of chromosome; used for groups, file names, and labels
# returns: list of effect plots as ggplot objects
all_sig_effect_plots <- function(chrom) {
  # loop over each hit's information
  apply(read_sig_hits(chrom), 1, function(x) {
    # preprocess RSEM for this tissue as done during analysis
    rsem <- load_rsem(x[['tissue']])
    rsem <- filter_rsem(rsem[, -c(1:3)], chrom)
    tmm <- rsem_to_cpm(rsem)
    
    effect_plot(x, tmm, groups[[chrom]], chrom)
  })
}

# arrange effect plots in a grid
# effect_plot_list: list of effect plots as ggplot objects
# start_i: what index in the alphabet to start labelling at
# ncol: how many columns the grid has
# nrow: how many rows the grid has
# returns: finished plot grid
arrange_plots <- function(effect_plot_list, start_i, ncol, nrow)
  plot_grid(plotlist = effect_plot_list, ncol = ncol, nrow = nrow,
            labels = LETTERS[start_i:(start_i + length(effect_plot_list) - 1)])