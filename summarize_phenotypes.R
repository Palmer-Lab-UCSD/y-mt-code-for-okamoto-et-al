# load test results
mlma_files <- list.files('results/phenotypes/', pattern = '.mlma', 
                         full.names = TRUE)

tests <- do.call('rbind', lapply(mlma_files, function(filename) {
    mlma <- read.csv(filename, sep = '\t')
    mlma$name <- gsub('.mlma', '', basename(filename))
    
    mlma
  }))

# separate by chromosome and save results
invisible(sapply(c('Y', 'MT'), function(chrom) {
  chr_tests <- tests[tests$SNP == paste0('type', chrom), 
                     c('name', 'Freq', 'p')]
  chr_tests$adj_p <- p.adjust(chr_tests$p, method = 'BH')
  
  write.csv(chr_tests[order(chr_tests$p), ], 
            paste0('results/associations/gwas_phenotype_', chrom, '_tests.csv'), 
            quote = FALSE, row.names = FALSE)
}))