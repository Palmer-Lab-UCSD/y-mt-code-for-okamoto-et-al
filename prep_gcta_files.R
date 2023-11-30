# ---- load requirements ----

library(genio) # write PLINK-style genotype files

# has RFIDs has sexes for all rats
geno_log <- read.csv('data/sample_info.csv')

# haplotype groups
y_groups <- read.csv('results/groups/Y_groups.csv')
mt_groups <- read.csv('results/groups/MT_groups.csv')

# encode haplotype as ref/alt genotype
y_groups$Y_group <- ifelse(y_groups$Y_group == 'Y1 (ACI)', 0, 2)
mt_groups$MT_group <- ifelse(mt_groups$MT_group == 'MT1 (BN)', 0, 2)
# general genotype table with RFIDs to haplotypes
geno <- merge.data.frame(y_groups, mt_groups, all = TRUE)

# ---- write genotype files ----

# file formats: https://www.cog-genomics.org/plink/1.9/formats

# .bed file is m variants (chrY/chrMT) by n rats
bed <- as.matrix(t(geno[, c(2, 3)]), nrow = 2)
rownames(bed) <- c('typeY', 'typeMT')

# .bim file is m variants by 6; mock up values for pseudo-SNPs
bim <- data.frame(chr = c(22, 23), id = c('typeY', 'typeMT'), posg = 0, 
                  pos = 1, ref = c('Y1', 'MT1'), alt = c('Y2', 'MT2'))

# .fam file is n rats by 6; ignore FID, pat/maternal ID, and phenotype value
sexes <- sapply(geno$rfid, function(iid) geno_log[geno_log$rfid == iid, 'sex'])
fam <- data.frame(fam = 0, id = geno$rfid, pat = 0, mat = 0, 
                  sex = ifelse(sexes == 'M', 1, 2), pheno = -9)

write_plink('results/groups/my_haplotypes', bed, bim = bim, fam = fam)

# ---- write phenotype files ----

pheno <- read.csv('data/phenotypes/gwas_phenotypes_table.csv', row.names = 1)

# write phenotype files for each trait
invisible(sapply(colnames(pheno), function(trait) {
  # use rats with non-NA values for this trait
  use_row <- !is.na(pheno[[trait]])
  
  # file format: https://yanglab.westlake.edu.cn/software/gcta/#GREMLanalysis
  write.table(data.frame(FID = 0, IID = rownames(pheno)[use_row], 
                         pheno = pheno[use_row, trait]), 
              paste0('data/phenotypes/', trait, '.txt'),
              sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
}))