[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11493119.svg)](https://doi.org/10.5281/zenodo.11493119)

This repository contains code to reproduce the paper "Y and Mitochondrial
Chromosomes in the Heterogeneous Stock Rat Population".

A prior version of this paper, the preprint Okamoto *et al.* 2023, is in bioRxiv
doi:[10.1101/2023.11.29.566473](https://doi.org/10.1101/2023.11.29.566473)

Input data files ("Data") for this code are in UCSD Library Digital Collections:
Okamoto, Faith; Chitre, Apurva S.; NIDA Center for GWAS in Outbred Rats; Palmer,
Abraham A. (2024). Data from: Y and MT Chromosomes in the Heterogeneous Stock
Rat Population. In The Center for GWAS in Outbred Rats Database (C-GORD). UC San
Diego Library Digital Collections. https://doi.org/10.6075/J0VX0GQQ

Software requirements are in `requirements.csv`.

# Reproducing this paper's analysis #

Do not run `genotype_helpers.py`, `io_helpers.py`, `gene_expression_helpers.R`.

- Core analysis (do #1 first, other two may occur in any order)
    
    1. Making haplogroups: `group_making.ipynb` creates haplogroup files.
    2. GWAS phenotypes: `gwas_phenotype_association.R` calculates p-values. 
    Note that this file has a part that is run pre-GCTA, then a part that is run
    post-GCTA. Notes are provided on how to run GCTA in between.
    3. Gene expression: `gene_expression_association.R` calculates p-values.

- SNP annotation (must be done in order, after making haplogroups)

    1. Finding SNPs: each of `y`/`mt_haplotype_id.ipynb` determines consensus
    sequences of modern haplotypes, then makes a VCF of SNPs between haplotypes.
    2. Build mRatBN7.2 in snpEff: 
        a. Make a subdirectory named `mRatBN7.2` within `snpEff/data`
        b. From [mRatBN7.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_015227675.2/)
        in GenBank, download Genome Sequences (FASTA), Protein (FASTA), and
        Annotation Features (GTF). Put these within the `mRatBN7.2` directory,
        named `sequences.fa`, `protein.fa`, and `genes.gtf`, respectively. Make
        sure these have consistent chromosome names, refSeq-style (for MT,
        `NC_001665.2`, for Y, `NC_051357.1`), and that they include the entire
        genome (as of May 14, the FTP site's GTF only has the MT Chromosome).
        c. Edit `snpEff.config` to append these lines:
        ```
        # Rat genome, version mRatBN7.2
        mRatBN7.2.genome : Rat
        mRatBN7.2.reference : data/mRatBN7.2
        mRatBN7.2.NC_001665.2.codonTable : Invertebrate_Mitochondrial
        ```
        d. Run `java -jar snpEff.jar build -gtf22 -v mRatBN7.2 -nocheckcds`. If
        this worked,`snpEff/data/mRatBN7.2/` should have many new `.bin` files
    3. Annotate SNPs using snpEff:
        a. Copy `Y_SNPs.vcf` and `MT_SNPs.vcf` into `snpEff/`.
        b. Run, within `snpEff/` and once each for Y/MT, the command
        `java -Xmx8g -jar snpEff.jar mRatBN7.2 Y_SNPs.vcf > Y_SNPs.ann.vcf`
        c. Copy annotated VCFs back to `results/groups/`

- Create intermediate files (must be done in order, after making haplogroups)

    1. `gene_expression_misc.R` saves subsets of normalized RNA-seq data
    2. `dkc1_x_eqtls.R` runs X-WAS eQTLs, saving p-values and the top SNP

- Misc. figures and tables: now, all other scripts may be run in any order.

# Directory structure #

The code expects this directory structure within your working directory:

- code
- data
    * depth
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

Some directories should start populated with files:

- `code/` has the code included herein.

- Directly within `data/` is `sample_info.csv` (basic information about each
modern HS rat), which is directly in Data.

- `data/depth/` has average read depth along given regions. These include:
`founders_near_Med14Y.csv`, `modern_deep_near_Med14Y.csv`,
`modern_shallow_in_Med14Y.csv`, `MT_depth.csv`, and `Y_depth.csv`. These files
are from "Sequencing depths" in Data.

- `data/expression/` has `log2` BED files, i.e. `<abbr>.rn7.expr.log2.bed.gz`,
for each tissue abbreviation (`<abbr>`). There are also "Brain" `tpm` and `iqn`.
Finally, lookup tables for gene/tissue names (`gene`/`tissue_name_map.csv`) and
cis contigs for each of Y and MT (`Y_MT_cis_contigs.csv`) are also included.
These files are from "Gene expression" in Data.

- `data/genotypes/` has all input genotypes. These include:
    * VCFs (documented in `genotype_helpers.py`, from "Raw genotypes" in Data):
    `HS_founders.vcf.gz`, `modern_HS_deep_sequenced.vcf.gz`,
    `modern_HS_shallow_sequenced.vcf.gz`, `STRs.vcf.gz`
    * PLINK-binary genotypes (also from "Raw genotypes") for X SNPs in a subset
    of modern HS rats: `chrX.bed`, `chrX.bim`, `chrX.fam`
    * `mRatBN_7_2_MT.fasta`: reference MT (also from "Raw genotypes")
    * GRM files (from "Genetic relationship matrix" in Data): 
    `autosomes.grm.bin`, `autosomes.grm.id`, `autosomes.grm.N.bin`

- `data/phenotypes/` has phenotypes, including (from "GWAS phenotypes" in Data):
    * `GWAS_phenotype_table.csv` (RFIDs x phenotypes, trait names anonymized)
    * `trait_dictionary.csv` (ties trait names to their original projects)
    * `kidneys.csv` (kidney count at birth for all rats tested)

Most scripts write to `results/` subdirectories. Exceptions:
- `gene_expression_misc.R` writes to `data/expression/`
- `dkc1_x_eqtls.R` writes `data/genotypes/X_GWAS_Dkc1_male_top_snp.csv`
- `gwas_phenotype_association.R` writes to `data/phenotypes/`

# Core output files #

Many output files are produced in the course of running the analysis scripts.
Some of them are likely to be of more immediate interest for those wishing to
directly use the results here presented. Here is a short list of files which
are not represented by a table within the text or supplements.
The given filenames are within `results/`.

| File | Description | Format |
|-|-|-|
| `groups/Y_groups.csv` | All rats with a Y haplogroup | CSV file with `rfid` and `Y_group` columns |
| `groups/MT_groups.csv` | All rats with a MT haplotype | CSV file with `rfid` and `MT_group` columns |
| `associations/Y_gene_expression_tests.csv` | Results of Y haplogroup to gene expression associations, sorted by p-value | CSV file with columns `ensembl_id`, `chr` [of gene],  `tissue` [using short names], `n1` [number of Y1 rats], `n2` [number of Y2 rats], `p`, and `adj_p` [from BH] | 
| `associations/MT_gene_expression_tests.csv` | Results of MT haplotype to gene expression associations, sorted by p-value | CSV file with columns `ensembl_id`, `chr` [of gene], `tissue` [using short names], `n1` [number of MT1 rats], `n2` [number of MT2 rats], `p`, and `adj_p` [from BH] | 
| `associations/Y_GWAS_phenotype_tests.csv` | Results of Y haplogroup to GWAS phenotype associations, sorted by p-value | CSV file with columns `name` [of phenotype], `Freq`, `p` [both output by GCTA], and `adj_p` [from BH] |
| `associations/MT_GWAS_phenotype_tests.csv` | Results of MT haplotype to GWAS phenotype associations, sorted by p-value | CSV file with columns `name` [of phenotype], `Freq`, `p` [both output by GCTA], and `adj_p` [from BH] |

Notably, the latter four files are from "Association test results" in Data.

# Sources of tables/figures #

Locations of output files are relative to/within `results/`. 
All script names are relative to/within `code/`.

Figures and tables are in order of appearance, with supplements at the end.
Some things not mentioned in the table include:
- The count of haplotyped and phenotyped rats is from `group_use.ipynb`
- The mention of "Ranking loses raw abundance information by introducing ties"
  in "Methods -> Gene expression association" is supported by a section of 
  `gene_expression_misc.R`

| Item | Output filename | Description | Script |
|-|-|-|-|
| Figure 1A | `plots/Y_founders_tree.svg` | NJ tree of HS founder Y | `nj_trees.ipynb` |
| Figure 1B | `plots/Y_ref_vs_alt.svg` | Count of reference vs. alternate Y alleles in modern and founder rats | `group_making.ipynb` |
| Figure 1C | `plots/Y_over_time.svg` | Percentage of births with each Y haplogroup over time | `group_use.ipynb` |
| Figure 1D | `plots/Y<1/2>_founders_Y<1/2>_vary_pos.svg` | Founder haplotypes at SNPs with founders vary | `y_haplotype_id.ipynb` |
| Figure 2A | `plots/MT_founders_tree.svg` | NJ tree of HS founder MT | `nj_trees.ipynb` |
| Figure 2B | `plots/MT_ref_vs_alt.svg` | Count of reference vs. alternate MT alleles in modern and founder rats | `group_making.ipynb` |
| Figure 2C | `plots/MT_over_time.svg` | Percentage of births with each MT haplotype over time | `group_use.ipynb` |
| Figure 2D | `plots/all_MT_MSA.svg` | Pseudo-alignment of HS founder MT haplotypes at SNPs genotyped in modern HS rats | `mt_haplotype_id.ipynb` |
| Figure 3A | `plots/Y_GWAS_phenotype_QQplot.svg` | QQ plot for p values of Y haplogroup vs. GWAS phenotype associations | `qqplots.ipynb` |
| Figure 3B | `plots/Y_gene_expression_QQplot.svg` | QQ plot for p values of Y haplogroup vs. gene expression associations | `qqplots.ipynb` |
| Figure 3C | `plots/ENSRNOG00000057231_Brain_effect_plot.svg` | Effect plot of *Ddx3y* expression in Brain by Y haplogroup | `effect_plots.ipynb` |
| Figure 3D | `plots/ENSRNOG00000055562_Brain_effect_plot.svg` | Effect plot of *Dkc1* expression in Brain by Y haplogroup | `effect_plots.ipynb` |
| Figure 4A | `plots/MT_GWAS_phenotype_QQplot.svg` | QQ plot for p values of MT haplotype vs. GWAS phenotype associations | `qqplots.ipynb` |
| Figure 4B | `plots/MT_gene_expression_QQplot.svg` | QQ plot for p values of MT haplotype vs. gene expression associations | `qqplots.ipynb` |
| Figure 4C | `plots/ENSRNOG00000033615_Brain_effect_plot.svg` | Effect plot of *Mt-nd3* expression in Brain by MT haplotype | `effect_plots.ipynb` |
| Figure 4D | `plots/ENSRNOG00000043866_Brain_effect_plot.svg` | Effect plot of *16S RNA* expression in Brain by MT haplotype | `effect_plots.ipynb` |
| Figure 4E | `plots/ENSRNOG00000029971_Brain_effect_plot.svg` | Effect plot of *Mt-nd5* expression in Brain by MT haplotype | `effect_plots.ipynb` |
| Figure 4F | `plots/ENSRNOG00000030478_Brain_effect_plot.svg` | Effect plot of *12S RNA* expression in Brain by MT haplotype | `effect_plots.ipynb` |
| Figure S1 | `plots/heterozygosity.png` | Distribution of heterozygosity across Y/MT | `group_use.ipynb` |
| Figure S2 | `plots/filtration_steps.png` | Distribution of statistics used to filter genotypes when haplotyping Y and MT | `group_making.ipynb` |
| Figure S3 | `plots/PheWAS_hist.png` | Distribution of sample sizes used in PheWAS | `group_use.ipynb` |
| Figure S4 | `plots/Y_depth.png` | Sequencing depth along Y | `snp_breakdown.ipynb` |
| Figure S5 | `plots/Y_STR_tree.svg` | NJ tree of founder and modern rats, based on STRs | `nj_trees.ipynb` |
| Figure S6 | `plots/SNP_effects.png` | Annotations for SNPs between haplotypes | `snp_breakdown.ipynb` |
| Figure S7 | `plots/Y_deviants.png` | Number of rats deviating from their Y consensus at each SNP | `y_haplotype_id.ipynb` |
| Figure S8 | `plots/PCA.png` | Biplot of Y/MT PCA | `group_use.ipynb` |
| Figure S9 | `plots/MT_missing.svg` | Distribution of missingness and restriction sites (for ddGBS) on MT | `group_use.ipynb` |
| Figure S10 | `plots/MT_depth.png` | Sequencing depth along MT | `snp_breakdown.ipynb` |
| Figure S11 | `plots/metabolite_effect_plot.png` | Effect plot of significant Y haplogroup association with a metabolite | `effect_plots.ipynb` |
| Figure S12A | `plots/X_WAS_Dkc1_male.svg` | Manhattan plot of X-WAS against Dkc1 expression in male brain | `manhattan_plots.ipynb` |
| Figure S12B | `plots/X_WAS_Dkc1_female.svg` | Manhattan plot of X-WAS against Dkc1 expression in female brain | `manhattan_plots.ipynb` |
| Figure S12C | `plots/Dkc1_top_SNP_effect_plot.svg` | Figure 3D colored by genotype at most-associated X SNP | `effect_plots.ipynb` |
| Figure S13 | many `plots/*_effect_plot.svg` | Effect plots of significant MT gene expression associations | `effect_plots.ipynb` |
| Table S1 | `associations/GWAS_phenotype_summary.csv` | Summary statistics of GWAS phenotypes used | `group_use.ipynb` |
| Table S2 | `associations/gene_expression_summary.csv` | Summary statistics of gene expression phenotypes used | `gene_expression_misc.R` |
| Table S3 | `associations/Y_gene_expression_tests_sig.csv` | Details of significant Y gene expression associations | `gene_expression_misc.R` |
| Table S4 | `associations/MT_gene_expression_tests_sig.csv` | Details of significant MT gene expression associations | `gene_expression_misc.R` |
| Table S5 | `associations/MT_vs_kidneys.csv` | Contingency table of MT haplotype and kidney count | `gwas_phenotype_association.R` |
| File S1 Figure A | `plots/Med14Y_expr_effect_plot.png` | Effect plot of *Med14Y* expression by tissue | `effect_plots.ipynb` |
| File S1 Figure B | `plots/Dkc1_expr_by_Med14Y_effect_plot.png` | Effect plot of *Dkc1* expression by tissue, colored by *Med14Y* expression group | `effect_plots.ipynb` |
| File S1 Figure C | `plots/founder_Med14Y_depth.png` | Sequencing depth in founders around *Med14Y* | `med14y_depth.ipynb` |
| File S1 Figure D | `plots/deep_modern_Med14Y_depth.png` | Sequencing depth in high-coverage modern HS rats around *Med14Y* | `med14y_depth.ipynb` |
| File S1 Figure E | `plots/modern_Med14Y_avg_depth.png` | Distribution of average read depth within *Med14Y* in low-coverage modern HS rats | `med14y_depth.ipynb` |
