[![DOI](https://zenodo.org/badge/725762461.svg)](https://zenodo.org/doi/10.5281/zenodo.10234037)

This repository contains code to replicate the preprint Okamoto F, Chitre AS, Sanches TM, Chen D, Munro D, NIDA Center for GWAS in Outbred Rats, Polesskaya O, Palmer AA 2023. Y and Mitochondrial Chromosomes in the Heterogeneous Stock Rat Population. *bioRxiv*. doi:[10.1101/2023.11.29.566473](https://doi.org/10.1101/2023.11.29.566473)

# Reproducing this paper's analysis #

- Y and MT haplotyping: all Python notebooks

    * First, run `group_making.ipynb` to create haplotype group files. 
    * The others can be run in any order. They make various panels and figures.
    * Don't run `genotype_helpers.py`; it is just a helper file.

- Gene expression association: `gene_expression_*.R`

    * First, run `gene_expression_association.R` to get p-values.
    * Then, run `gene_expression_misc.R` for a few odds and ends.
    * Don't run any others; they are just helper files.

- Ex-GWAS trait association: 

    * First, run `prep_gcta_files.R` to create GCTA's genotype/phenotype files.
    * Then, run `phenotype_association.sh` to get p-values.
    * Finally, run `summarize_phenotypes.R` to consolidate GCTA's output files.
    * At any point, run `phenotype_misc.R` for some odds and ends.

- Final figures: `association_figures.R`  
  Run after `gene_expression_association.R` and `summarize_phenotypes.R`.

# Directory structure #

The code expects this directory structure within your working directory:

- code
- data
    * expression
    * genotypes
    * phenotypes
- results
    * alignments
    * associations
    * groups
    * phenotypes
    * plots
    * trees

`code/` starts populated with the code included herein.

`data/expression/` starts populated with `log2` BED files for each tissue,
i.e. `<name>.rn7.expr.log2.bed.gz`, where `<name>` is a tissue's short name.
There are also `tpm` and `iqn` files for only "Brain" data. These files may be
found in the `gene_expression.tar` tarball.

`data/genotypes/` starts populated with all input genotypes. These include:
- VCFs (documented in `genotype_helpers.py`, stored in `y_mt_genotypes.tar.gz`)
    * `HS_founders.vcf.gz`
    * `modern_HS_deep_sequenced.vcf.gz`
    * `modern_HS_shallow_sequenced.vcf.gz`
    * `STRs.vcf.gz`
- GRM files (stored in `autosomal_grm.tar.gz`)
    * `autosomes.grm.bin` 
    * `autosomes.grm.id`
    * `autosomes.grm.N.bin`
- `mRatBN_7_2_mt.fasta`: reference MT (also from `y_mt_genotypes.tar.gz`)
- `mt_depth.csv` (read depth along MT for each low-coverage sample)

`data/phenotypes/` starts populated with consolidated phenotype files. These
include (stored in `gwas_phenotypes.tar.gz`):
- `gwas_phenotypes_table.csv` (RFIDs x phenotypes, trait names anonymized)
- `trait_dictionary.csv` (ties trait names to their original projects)
- `project_metadata.csv` (ties project name to metadata)
- `kidneys.csv` (kidney count at birth for all rats tested)

Directly within `data/` is `sample_info.csv` (basic information about each
modern HS rat: RFID, sex, sequencing method, and, if known, birth date)

All scripts write to `results/` subdirectories, with one exception:
`prep_gcta_files.R` populates `data/phenotypes/`

# Core output files #

Many output files are produced in the course of running the analysis scripts.
Some of them are likely to be of more immediate interest for those wishing to
directly use the results here presented. Here is a short list of files which
are not represented by a table within the text or supplements.
The given filenames are relative to/within `results/`.

| File | Description | Format |
|-|-|-|
| `groups/Y_groups.csv` | All rats with a Y haplotype group | CSV file with `rfid` and `Y_group` columns |
| `groups/MT_groups.csv` | All rats with a MT haplotype group | CSV file with `rfid` and `MT_group` columns |
| `associations/gene_expression_Y_tests.csv` | Results of Y haplotype to gene expression associations, sorted by p-value | CSV file with `ensembl_id`, `tissue` [using short names], `chr` [of gene], `n1` [number of Y1 rats], `n2` [number of Y2 rats], `p`, and `adj_p` [from BH] columns | 
| `associations/gene_expression_MT_tests.csv` | Results of MT haplotype to gene expression associations, sorted by p-value | CSV file with `ensembl_id`, `tissue` [using short names],  `chr` [of gene], `n1` [number of MT1 rats], `n2` [number of MT2 rats], `p`, and `adj_p` [from BH] columns | 
| `associations/gwas_phenotype_Y_tests.csv` | Results of Y haplotype to ex-GWAS phenotype associations, sorted by p-value | CSV file with `name` [of phenotype], `Freq`, `p` [both output by GCTA], and `adj_p` [from BH] columns |
| `associations/gwas_phenotype_MT_tests.csv` | Results of MT haplotype to ex-GWAS phenotype associations, sorted by p-value | CSV file with `name` [of phenotype], `Freq`, `p` [both output by GCTA], and `adj_p` [from BH] columns |

Notably, the latter four files are stored in `association_results.tar.gz`.

# Sources of tables/figures #

Locations of output files are relative to/within `results/`. 
All script names are relative to/within `code/`.

Figures and tables are in order of appearance within the text.
Some things not mentioned in the table include:
- The count of haplotyped and phenotyped rats is from `phenotypes_misc.R`
- The mention of "Ranking loses raw abundance information by introducing ties"
  in "Methods -> Gene expression association" is supported by a section of 
  `gene_expression_misc.R`

| Item | Output filename | Description | Script |
|-|-|-|-|
| Figure 1A | `plots/Y_founders_tree.svg` | NJ tree of HS founder Y | `nj_trees_trees.ipynb` |
| Figure 1B | `plots/Y_ref_vs_alt.svg` | Count of reference vs. alternate Y alleles in modern and founder rats | `group_making.ipynb` |
| Figure 1C | `plots/Y_over_time.svg` | Percentage of births with the Y1 haplotype over time | `group_use.ipynb` |
| Figure 1D | `plots/Y_STR_tree.svg` | NJ tree of founder and modern rats, based on STRs | `nj_trees_trees.ipynb` |
| Figure 1E (bars) | `plots/Y1_entropy.svg` | Strength of consensus Y1 haplotype at SNPs where Y1 founders vary | `y_haplotype_id.ipynb` |
| Figure 1E (alignment) | `plots/Y1_align.svg` | Consensus Y1 haplotype at SNPs where Y1 founders vary | `y_haplotype_id.ipynb` |
| Figure 1F | `plots/Y1_mismatches.svg` | Number of Y1 rats deviating from their consensus at each SNP | `y_haplotype_id.ipynb` |
| Figure 1G (bars) | `plots/Y2_entropy.svg` | Strength of consensus Y2 haplotype at SNPs where Y2 founders vary | `y_haplotype_id.ipynb` |
| Figure 1G (alignment) | `plots/Y2_align.svg` | Consensus Y2 haplotype at SNPs where Y2 founders vary | `y_haplotype_id.ipynb` |
| Figure 1H | `plots/Y2_mismatches.svg` | Number of Y2 rats deviating from their consensus at each SNP | `y_haplotype_id.ipynb` |
| Figure 2A | `plots/MT_founders_tree.svg` | NJ tree of HS founder MT | `nj_trees.ipynb` |
| Figure 2B | `plots/MT_over_time.svg` | Percentage of births with the MT1 haplotype over time | `group_use.ipynb` |
| Figure 2C | `plots/MT_ref_vs_alt.svg` | Count of reference vs. alternate MT alleles in modern and founder rats | `group_making.ipynb` |
| Figure 2D | `plots/all_MT_msa.svg` | Pseudo-alignment of modern MT haplotypes along with HS founder MT | `mt_haplotype_id.ipynb` |
| Figure 3 | `plots/Y_association.svg` | Results of Y haplotype association tests | `association_figures.R` |
| Figure 4 | `plots/MT_association.svg` | Results of MT haplotype association tests | `association_figures.R` |
| Table 1 | `associations/gene_expression_Y_sig_tests.csv` | Details of significant Y gene expression associations | `gene_expression_misc.R` |
| Table 2 | `associations/gene_expression_MT_sig_tests.csv` | Details of significant MT gene expression associations | `gene_expression_misc.R` |
| Figure S1 | `plots/heterozygosity.svg` | Distribution of heterozygosity across Y/MT | `group_use.ipynb` |
| Figure S2 | `plots/filtration_steps.svg` | Distribution of statistics used to filter genotypes when haplotyping Y and MT | `group_making.ipynb` |
| Figure S3 | `plots/pca.svg` | Biplot of Y/MT PCA | `group_use.ipynb` |
| Figure S4 | `plots/MT_missing.svg` | Distribution of missingness and restriction sites (for ddGBS) on MT |
| Figure S5 | `plots/metabolite.svg` | Effect plot of significant Y haplotype association with a metabolite | `phenotypes_misc.R` |
| Figure S6 | `plots/sig_MT_effect_plots.svg` | Effect plots of significant MT gene expression associations | `gene_expression_effect_plots.R` |
| Table S1 | `associations/phenotype_summary.csv` | Summary statistics of ex-GWAS phenotypes used | `phenotypes_misc.R` |
| Table S2 | `associations/gene_expression_summary.csv` | Summary statistics of gene expression phenotypes used | `gene_expression_misc.R` |
| Table S3 | `phenotypes/mt_vs_kidneys.csv` | Contingency table of MT haplotype and kidney count | `phenotypes_misc.R` |
