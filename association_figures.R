source('code/gene_expression_effect_plots.R')
source('code/qqplots.R')

y_effect_plots <- all_sig_effect_plots('Y')
mt_effect_plots <- all_sig_effect_plots('MT')

# Figure 3A-B
y_qqplot <- combined_qqplot('Y', c('Y', 'JACYVU010000493.1', 'MU150190.1', 
                                   'MU150191.1', 'MU150192.1', 'MU150193.1'),
                            c('Ddx3y', 'Dkc1'))
# Figure 4A-B
mt_qqplot <- combined_qqplot('MT', 'MT', c('Mt-nd3', '16S rRNA', 'Mt-nd5'))

# Figure 3C-D
all_y_effect_plots <- plot_grid(plotlist = y_effect_plots, 
                                ncol = 2, nrow = 1, labels = c('C', 'D'))


# Figure 3
y_assoc <- plot_grid(y_qqplot, all_y_effect_plots, ncol = 1, nrow = 2,
                     rel_heights = c(2, 1))
save_plot('results/plots/Y_association.svg', y_assoc, ncol = 1, nrow = 2,
          base_height = 4, base_asp = 3)

# Figure 4C-F
subset_mt_effect_plots <- plot_grid(plotlist = mt_effect_plots[c(1, 2, 3, 7)], 
                                    ncol = 2, nrow = 2, 
                                    labels = c('C', 'D', 'E', 'F'))

# Figure 4
mt_assoc <- plot_grid(mt_qqplot, subset_mt_effect_plots, ncol = 1, nrow = 2)
save_plot('results/plots/MT_association.svg', mt_assoc, ncol = 1, nrow = 2,
           base_height = 5.5, base_asp = 2)

# gene IDs for MT hits
mt_genes <- read_sig_hits('MT')$ensembl_id

# group MT IDs by gene type
nd3_id <- names(gene_map)[gene_map == 'Mt-nd3']
other_nd_id <- names(gene_map)[startsWith(gene_map, 'Mt-nd')]
other_nd_id <- other_nd_id[other_nd_id != nd3_id]
rrna_id <- names(gene_map)[endsWith(gene_map, 'rRNA')]

# order plots as rRNA, then Mt-nd3, then other Mt-nd sub-units
all_mt_order <- c(which(mt_genes %in% rrna_id), which(mt_genes == nd3_id), 
                  which(mt_genes %in% other_nd_id))

# Figure S6
all_mt_effect_plots <- plot_grid(plotlist = mt_effect_plots[all_mt_order], 
                                 ncol = 3, nrow = 6, labels = 'AUTO')
save_plot('results/plots/all_MT_effect_plots.svg', all_mt_effect_plots,
          ncol = 3, nrow = 6)
