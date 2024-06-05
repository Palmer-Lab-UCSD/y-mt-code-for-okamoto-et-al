# ---- load requirements ----

library(genio) # write PLINK-binary genotype files

sample_info <- read.csv('data/sample_info.csv', row.names = 'rfid')

y_groups <- read.csv('results/groups/Y_groups.csv')
y_groups$Y_group <- ifelse(y_groups$Y_group == 'Y1 (ACI)', 0, 2)
mt_groups <- read.csv('results/groups/MT_groups.csv')
mt_groups$MT_group <- ifelse(mt_groups$MT_group == 'MT1 (BN)', 0, 2)

geno <- merge.data.frame(y_groups, mt_groups, all = TRUE)

pheno <- read.csv('data/phenotypes/GWAS_phenotype_table.csv', row.names = 1)
kidneys <- read.csv('data/phenotypes/kidneys.csv')

# ==== Analysis prep ====

# ---- write genotype files ----

# file formats: https://www.cog-genomics.org/plink/1.9/formats

bed <- as.matrix(t(geno[ , c(2, 3)]), nrow = 2)
rownames(bed) <- c('typeY', 'typeMT')

bim <- data.frame(chr = c(22, 23), id = c('typeY', 'typeMT'), posg = 0, 
                  pos = 1, ref = c('Y1', 'MT1'), alt = c('Y2', 'MT2'))

sexes <- sample_info[geno$rfid, 'sex']
fam <- data.frame(fam = 0, id = geno$rfid, pat = 0, mat = 0, 
                  sex = ifelse(sexes == 'M', 1, 2), pheno = -9)

write_plink('results/groups/my_haplotypes', bed, bim = bim, fam = fam)

# ---- write phenotype files ----

sapply(colnames(pheno), function(trait) {
  use_row <- !is.na(pheno[[trait]])
  pheno_df <- data.frame(FID = 0, IID = rownames(pheno)[use_row], 
                         pheno = pheno[use_row, trait])
  
  # file format: https://yanglab.westlake.edu.cn/software/gcta/#GREMLanalysis
  write.table(pheno_df, paste0('data/phenotypes/', trait, '.txt'), sep = '\t',
              quote = FALSE, row.names = FALSE, col.names = FALSE)
})

# ==== Phenotype analysis ====

# run PheWAS using GCTA: for each phenotype (e.g. <name>=`pheno1`), use options
# --mlma --bfile results/groups/my_haplotypes --grm data/genotypes/autosomes
# --pheno data/phenotypes/<name>.txt --out results/phenotypes/<name>

# ---- summarize result files ----

mlma_files <- list.files('results/phenotypes/', pattern = '.mlma',
                         full.names = TRUE)

tests <- do.call('rbind', lapply(mlma_files, function(filename) {
    mlma <- read.csv(filename, sep = '\t')
    mlma$name <- gsub('.mlma', '', basename(filename))
    
    mlma
  }))

for (chrom in c('Y', 'MT')) {
  chr_tests <- tests[tests$SNP == paste0('type', chrom), c('name', 'Freq', 'p')]
  chr_tests$adj_p <- p.adjust(chr_tests$p, method = 'BH')

  write.csv(chr_tests[order(chr_tests$p), ], 
            paste0('results/associations/', chrom, '_GWAS_phenotype_tests.csv'), 
            quote = FALSE, row.names = FALSE)
}

# ---- test for kidney association ----

mt_kidney <- merge.data.frame(mt_groups, kidneys)
mt_kidney <- aggregate(rfid ~ MT_group + n_kidneys, mt_kidney, length)

dim_names <- list(c('MT1 (BN)', 'MT2 (nearly ACI)'), c('1 kidney', '2 kidneys'))
mt_kidney <- matrix(mt_kidney$rfid, nrow = 2, dimnames = dim_names)

# table S5
write.csv(mt_kidney, 'results/associations/MT_vs_kidneys.csv', quote = FALSE)
fisher.test(mt_kidney, alternative = 'greater')
