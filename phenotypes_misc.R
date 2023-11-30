# ---- load requirements ----

library(ggplot2) # effect plots
library(cowplot) # plot theme

source('code/gene_expression_helpers.R')

pheno <- read.csv('data/phenotypes/gwas_phenotypes_table.csv', row.names = 1)
traits_by_project <- read.csv('data/phenotypes/trait_dictionary.csv')

# haplotype groups
y_groups <- read.csv('results/groups/Y_groups.csv')
mt_groups <- read.csv('results/groups/MT_groups.csv')

kidneys <- read.csv('data/phenotypes/kidneys.csv')

# ---- test for kidney association ----

mt_kidney <- merge.data.frame(mt_groups, kidneys)
# count RFIDs in MT1/one kidney, MT2/one, MT1/two kidneys, MT2/two
mt_kidney <- aggregate(rfid ~ MT_group + n_kidneys, mt_kidney, length)

# convert to nice matrix
mt_kidney <- matrix(mt_kidney$rfid, nrow = 2, dimnames = 
                      list(c('MT1 (BN)', 'MT2 (nearly ACI)'),
                           c('1 kidney', '2 kidneys')))

# table S3
write.csv(mt_kidney, 'results/phenotypes/mt_vs_kidneys.csv', quote = FALSE)

fisher.test(mt_kidney, alternative = 'greater')

# ---- count rats ----

has_haplotype <- union(y_groups$rfid, mt_groups$rfid)
has_gwas_pheno <- rownames(pheno)
has_gene_expr <- sapply(names(tissue_map), function(tissue) 
  colnames(load_rsem(tissue))[-c(1:3)])
has_gene_expr <- unique(unlist(has_gene_expr, use.names = F))
has_pheno <- union(has_gwas_pheno, has_gene_expr)

paste(length(has_haplotype), 'rats have at least one haplotype,',
      length(has_gwas_pheno), 'rats have at least one GWAS phenotype,',
      length(has_gene_expr), 'rats have at least one RNA-seq sample.',
      length(has_pheno), 'rats have at least one phenotype, and',
      length(intersect(has_haplotype, has_pheno)), 
      'have both a phenotype and a haplotype')

# ---- table S1 ----

# summarize trait-project combinations as traits/rats by project
n_traits <- aggregate(trait ~ project, traits_by_project, length)
colnames(n_traits) <- c('project', 'traits')

# count haplotypes in each project
n_rats <- do.call('rbind', lapply(n_traits$project, function(name) {
  traits <- traits_by_project[traits_by_project$project == name, 'trait']
  proj_pheno <- pheno[, traits]
  # remove rats with no phenotypes in this project
  rfids <- rownames(proj_pheno)[rowSums(!is.na(proj_pheno)) > 0]
  
  # will list n Y1, Y2, MT1, MT2
  counts <- c(table(y_groups[y_groups$rfid %in% rfids, 'Y_group']),
              table(mt_groups[mt_groups$rfid %in% rfids, 'MT_group']))
  setNames(counts, c('Y1', 'Y2', 'MT1', 'MT2'))
}))

write.csv(cbind(n_traits, n_rats), 'results/associations/phenotype_summary.csv', 
          quote = FALSE, row.names = FALSE)

# ---- figure S5 ----

sig_name <- 'MZ531.3646417_5.08009'
sig_pheno <- pheno[!is.na(pheno[[sig_name]]), sig_name, drop = FALSE]

# look up haplotype for each sample, and count n for each group
to_plot <- merge.data.frame(sig_pheno, y_groups, by.x = 'row.names', 
                            by.y = 'rfid')
colnames(to_plot)[2] <- 'metabolite'
to_plot$Y_group <- paste0(to_plot$Y_group, '\nn = ', 
                          table(to_plot$Y_group)[to_plot$Y_group])
# make a factor for plotting purposes
to_plot$Y_group <- factor(to_plot$Y_group)

# plot values on the y axis, separating haplotypes along the x
ggplot(to_plot, aes(Y_group, metabolite)) + 
  
  # violin plot with box-plot quantiles drawn across
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), scale = 'width') + 
  # jitter points only to the side, to keep them at the correct y position
  geom_jitter(height = 0, width = 0.2) + 
  
  # clean up visual appearance, then blow up tick labels for readability
  theme_minimal_hgrid() + theme(axis.text.x = element_text(size = 15)) +
  
  labs(title = paste('Levels of', sig_name), x = '', y = 'normalized residuals')
ggsave('results/plots/metabolite.svg')
