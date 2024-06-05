# ---- setup ----

library(genio) # read PLINK-binary genotype files

# "X" is a matrix corresponding to the BED file
geno <- read_plink('data/genotypes/chrX')$X
# Dkc1 expression in Brain; see Figure 3D
expr <- read.csv('data/expression/ENSRNOG00000055562_Brain_expr.csv')
expr <- setNames(expr$expr, nm = expr$rfid)

sample_info <- read.csv('data/sample_info.csv')
y_groups <- read.csv('results/groups/Y_groups.csv')
sample_info <- merge(sample_info, y_groups, by = 'rfid', all = TRUE)
# convert haplotype group into integer, for use as covariate
sample_info$Y_group <- as.numeric(substr(sample_info$Y_group, 2, 3))
# females need a uniform non-NA Y group to use as covariate
sample_info[sample_info['sex'] == 'F', 'Y_group'] <- 0

# run linear regression of Dkc1 vs. one SNP + Y group covariate
# snp_geno: SNP genotypes encoded as 0/1/2 and labeled by RFIDs
# sex: which sex of RFIDs to analyze
# returns: p-value of the regression
regress_snp_dkc1 <- function(snp_geno, sex) {
  snp_geno <- snp_geno[!is.na(snp_geno)]
  sex_rfids <- sample_info[sample_info$sex == sex, 'rfid']
  ids_to_test <- intersect(names(snp_geno), sex_rfids)
  
  # only bother with the test if there are genotypes
  if(length(ids_to_test) == 0) { return(NA) }
  
  snp_geno <- snp_geno[ids_to_test]
  pheno <- expr[ids_to_test]
  covariates <- sample_info[sample_info$rfid %in% ids_to_test, 'Y_group']
  
  # only bother with the test if genotype varies
  if(max(table(snp_geno)) > length(ids_to_test) - 5) { return(NA) }
  
  model <- lm(pheno ~ snp_geno + covariates)
  unname(summary(model)$coefficients[2, 4])
}

# run and save a chrX-WAS of Dkc1 expression using simple linear regression
# sex: which sex of RFIDs to analyze
# returns: data.frame with SNP name and regression p-value
regress_all_x_snps_dkc1 <- function(sex) {
  pvals_list <- apply(geno, 1, regress_snp_dkc1, sex)
  pvals_df <- data.frame(SNP = names(pvals_list), p = unname(pvals_list))

  long_sex_name <- ifelse(sex == 'F', 'female', 'male')
  filename <- paste0('results/associations/X_WAS_Dkc1_', long_sex_name, '.csv')
  write.csv(pvals_df, filename, row.names = FALSE)
  pvals_df
}

# ---- run analysis -----

male_pvals <- regress_all_x_snps_dkc1('M')
# no need to save female results
invisible(regress_all_x_snps_dkc1('F'))

male_pvals <- male_pvals[!is.na(male_pvals$p),]
top_snp <- male_pvals[male_pvals$p == min(male_pvals$p), 'SNP'][1]

top_geno <- data.frame(colnames(geno), geno[top_snp, ])
colnames(top_geno) <- c('rfid', top_snp)
write.csv(top_geno, 'data/genotypes/X_WAS_Dkc1_male_top_SNP.csv', 
          row.names = FALSE)
