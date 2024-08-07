"""Helper for genotype processing notebooks.

Stores genotype-related functions/constants used in multiple notebooks.

< CONSTANTS >

- `Genotypes`: a type alias for `pd.DataFrame` which additionally
    implies the data are formatted to represent genotypes

- VCF files with all Y and MT genotypes used (see "Genotype datasets")
    * `SHALLOW_MODERN_FILE`: raw shallow-sequenced Y and MT genotypes
        from all our modern HS rats; only biallelic SNPs
    * `DEEP_MODERN_FILE`: processed deep-sequenced Y and MT genotypes
        from a subset of modern HS rats; SNPs and indels
    * `FOUNDERS_FILE`: processed deep-sequenced Y and MT genotypes 
        from a male rat of each HS founder strain; SNPs and indels
    * `STR_FILE`: processed deep-sequenced Y STR genotypes
        from the deep-sequenced modern HS rat subset and male founders

- `FOUNDER_NAMES`: short names of all HS founder strains
    as they appear in `FOUNDERS_FILE` and `STR_FILE`
- `CHR_NAMES`: lookup table with common chromosome names (Y/MT) to
    their RefSeq identifiers (NC_051357.1/NC_001665.2)

< FUNCTIONS >

All public functions have some basic imput validation when necessary.
Many also have default or optional arguments; see docstrings.

- `get_genotypes()`: read (and filter) a VCF file into a `Genotypes`
- `plot_align()`: plot a pseudo-alignment
- `filter_maf()`: filter genotypes to remove variants with MAF=0
- `are_variants_maf_0()`: check if variants in a `Genotypes` have MAF=0
"""

#region ---- Imports ----

import errno # error handling
import os # check for file existance

from typing import Tuple
from numpy.typing import ArrayLike

import allel # read VCF files
import matplotlib.pyplot as plt # make simple visualizations
import numpy as np # data management
import pandas as pd # data management
import pymsaviz # pseudo-alignment

import io_helpers # file I/O

print(f'''
This genotype helper file adds:
NumPy version {np.__version__},
pyMSAviz version {pymsaviz.__version__},
scikit-allel version {allel.__version__}
''')

#endregion
#region ---- Constants -----

Genotypes = pd.DataFrame
"""
`pd.DataFrame` with a shape of (samples x variants).

The index is labeled with sample names:
- Modern HS rats use RFIDs
- HS founders use strain short names (see `FOUNDER_NAMES`)

The columns are labeled with base-pair position.

Values may be in one of two paradigms:

- Numeric, i.e. `format="allele_num"` in `get_genotypes()`

    Alleles are numbered starting from 0, which is the reference.
    Missing genotypes will have the value `np.nan`. 
    The VCF files are diploid, i.e. two alleles per sample per locus.
    The alleles are combined into a single number as so:
    if the genotype call is `a/b` (unphased), they become `a * 10 + b`.

- String, i.e. `format="nucleotide"` in `get_genotypes()`

    Alleles are one-letter nucleic acid codes (A, C, G, T).
    Missing genotypes will have the value `None`.
    The VCF files are diploid, but the representation used is haploid.
    Consequently, this format cannot handle heterozygosity.

This orientation was chosen for consistency with SciPy functions:
- individual observations (here, samples) are rows
- cross-observation dimensions (here, variants) as columns

See `pd.DataFrame` documentation.
-----
"""

# plots are stored with their letter tags
_Plot = Tuple[str, plt.Axes]

# names of the items represented by labels on `Genotypes` axes
_GENOTYPES_AXIS_ITEM = {'index': 'sample', 'columns': 'variant'}

SHALLOW_MODERN_FILE = io_helpers.file_path('data', 'genotypes', 
                                          'modern_HS_shallow_sequenced.vcf.gz')
DEEP_MODERN_FILE = io_helpers.file_path('data', 'genotypes', 
                                       'modern_HS_deep_sequenced.vcf.gz')
FOUNDERS_FILE = io_helpers.file_path('data', 'genotypes', 'HS_founders.vcf.gz')
STR_FILE = io_helpers.file_path('data', 'genotypes', 'STRs.vcf.gz')

FOUNDER_NAMES = ['ACI', 'BUF', 'BN', 'F344', 'M520', 'MR', 'WKY', 'WN']
CHR_NAMES = {'Y': 'NC_051357.1', 'MT': 'NC_001665.2'}

sexes = io_helpers.read_col('data', 'sample_info.csv', colname='sex', 
                           index_col='rfid')
_MALES = sexes.index[sexes == 'M']
# all founder samples are male
_MALES = _MALES.to_list() + FOUNDER_NAMES

_VARIANT_MIN_NONMISS_PERCENT = 75
_SAMPLE_MIN_NONMISS_PERCENT = 50
_NUM_ALT_SEQUENCES_READ = 8
# allows all SNPs through
_MT_MIN_INFO_SCORE = 0.85
# more leniant than normal (0.9) to accommodate females
_Y_MIN_INFO_SCORE = 0.5

#endregion
#region ---- Public functions ----

def get_genotypes(vcf_file: str, chrom: str, format: str, 
                  variant_subset: ArrayLike = None, allow_het: bool = False,  
                  use_maf_filter: bool = True, use_miss_filter: bool = True, 
                  filter_plots: list[_Plot] = None) -> Genotypes:
    """Read a VCF file into a `Genotypes`, with optional filters.

    This function handles everything necessary to go from a VCF file
    to a fully-formatted, fully-filtered `Genotypes` data frame.

    If this is the Y chromosome, samples are subset to only males.
    If these are shallow-sequenced genotypes, an INFO score filter
    is applied; see `_MT_`/`_Y_MIN_INFO_SCORE` for thresholds.

    Defaults (all of these have options):
    - Variant positions are not subset before filtration.
    - Heterozygosity is treated as missing.
    - MAF and missingness filters are applied.
    - Diagnostic plots are not made for filters.
    
    Parameters
    ----------
    vcf_file: str
        VCF file to read data from. One of the filename constants
        exported by this module.
    chrom: str
        Chromosome to read data for; either "Y" or "MT".
    format: str
        Format to return data in; either "nucleotide" or "allele_num".
        See the `Genotypes` docstring for details.
    variant_subset: ArrayLike, optional
        If specified, a list of variant positions to subset data to.
        Removing other variants occurs before statistic-based filters.
    allow_het: bool, default=False
        Whether to treat heterozygous calls as valid or missing.
    use_maf_filter: float, default=True
        Whether to remove variants with MAF=0.
    use_miss_filter: float, default=True
        Whether to remove variants with too much missingness,
        and then samples with too much missingness.
        See, respectively, `_VARIANT_MIN_NONMISS_PERCENT` and 
        `_SAMPLE_MIN_NONMISS_PERCENT` for the minimum
        genotyping rate thresholds used.
    filter_plots: list[_Plot], optional
        If specified, a list of (tag, Axes) pairs to plot filter 
        diagnostics (i.e. distribution of test statistics) on.
        This option assumes that all filters will be used and
        thus expects to receive plots for all four filters.
        This option further assumes that the genotype data
        consist of purely SNPs, for the purpose of plot labels.
    
    Returns
    -------
    geno: Genotypes
        Genotypes as a (samples x variants) shape data frame.
    """

    # ---- input validation -----

    if not os.path.isfile(vcf_file):
        raise FileNotFoundError(
            errno.ENOENT, os.strerror(errno.ENOENT), vcf_file
            )
    if chrom != 'Y' and chrom != 'MT':
        raise ValueError('`chrom` must be either "Y" or "MT":'
                         f' `{chrom}` is not a valid value')
    if format != 'nucleotide' and format != 'allele_num':
        raise ValueError('`format` must be either "nucleotide" or "allele_num":'
                         f' `{format}` is not a valid value')
    if format == 'nucleotide' and allow_het:
        raise ValueError('`format="nucleotide"` is incompatible with'
                         ' `allow_het=True`')
    if filter_plots is not None:
        if len(filter_plots) != 4:
            raise ValueError('`filter_plots`, if set, must be of length 4:'
                             f' it is length {len(filter_plots)}')
        if not (vcf_file == SHALLOW_MODERN_FILE 
                and use_maf_filter and use_miss_filter):
            raise ValueError('If `filter_plots` is set, all filters (INFO score'
                             ', MAF, variant/sample missingness) must be used; '
                             'thus shallow-sequenced genotypes are required.')
        
    # ---- read raw data ----
    
    has_info_scores = vcf_file == SHALLOW_MODERN_FILE
    use_plots = filter_plots is not None

    geno, info_scores = _vcf_to_df(vcf_file, chrom, has_info_scores, format)

    if chrom == 'Y':
        geno = _subset_axis_by_label(geno, _MALES, axis='index', why='male')
    if variant_subset is not None:
        geno = _subset_axis_by_label(geno, variant_subset, axis='columns',
                                     why='in the given subset')

    # only for allele_num since nucleotide heterozygosity was already removed
    if format == 'allele_num': geno = _set_allele_num_missing(geno, allow_het)

    # ---- apply filters -----

    if has_info_scores:
        # Y's peak of good INFO scores is lower, thus, different threshold
        if chrom == 'MT': min_info_score = _MT_MIN_INFO_SCORE
        else: min_info_score = _Y_MIN_INFO_SCORE

        if use_plots:
            _plot_hist(info_scores, cutoff_val=min_info_score, 
                       cutoff_label=min_info_score, 
                       title=f'INFO score, by SNP in {chrom}', 
                       x_label='INFO score', y_label='# SNPs', 
                       plot=filter_plots[0])
        geno = _filter_axis_by_cutoff(geno, vals=info_scores, 
                                      min_val=min_info_score, axis='columns', 
                                      statistic='INFO score')
                     
    if use_maf_filter:
        if use_plots: _plot_maf(geno, chrom, plot=filter_plots[1])
        geno = filter_maf(geno)
    
    if use_miss_filter:
        geno = _filter_miss_by_axis(
            geno, min_nonmiss_percent=_VARIANT_MIN_NONMISS_PERCENT, 
            by_axis='columns', chrom=chrom, 
            plot = filter_plots[2] if use_plots else None
            )
        geno = _filter_miss_by_axis(
            geno, min_nonmiss_percent=_SAMPLE_MIN_NONMISS_PERCENT, 
            by_axis='index', chrom=chrom,
            plot = filter_plots[3] if use_plots else None
            )
        
    print(f'Genotypes for {geno.shape[0]} rats across {geno.shape[1]} variants')
    
    return geno

def write_snp_vcf(snps: pd.DataFrame, chrom: str) -> None:
    """Write a set of SNPs into a sample-less VCF.

    Takes SNPs presumed to be between haplotypes for some
    chromosome and writes them to a file in the VCF format.
    
    Parameters
    ----------
    snps: pd.DataFrame
        SNPs to write to the file. Must have columns `POS`
        (1-indexed locus), `REF`, and `ALT` (ref/alt alleles).
    chrom: str
        Chromosome these data are for; either "Y" or "MT".
    """ 

    if chrom != 'Y' and chrom != 'MT':
        raise ValueError('`chrom` must be either "Y" or "MT":'
                         f' `{chrom}` is not a valid value')

    for col in ['POS', 'REF', 'ALT']:
        if col not in snps.columns:
            raise ValueError(f'snps must include a {col} column')

    snps['#CHROM'] = CHR_NAMES[chrom]
    snps['ID'] = chrom + ':' + snps['POS'].astype(str)
    snps['QUAL'] = 50
    snps['FILTER'] = 'PASS'
    snps['INFO'] = '.'

    # VCF format: https://samtools.github.io/hts-specs/VCFv4.2.pdf
    snps = snps[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']]

    file = io_helpers.file_path('results', 'groups', f'{chrom}_SNPs.vcf')
    with open(file, 'w') as vcf:
        vcf.write(f'##fileformat=VCFv4.2\n##contig=<ID={CHR_NAMES[chrom]}>\n')
        snps.to_csv(vcf, sep='\t', index=False, lineterminator='\n')

def plot_align(geno: Genotypes, basename: str, 
               wrap_length: int = None) -> plt.Figure:
    """Plot a pseduo-alignment of SNPs.

    See "Visualizing SNP pseudo-alignments".
    Given SNP genotypes, write a .fasta file, then 
    use that to plot a corresponding `pyMSAviz` object.
    Suppresses the default position-number ticks.

    Parameters
    ----------
    geno: Genotypes
        Genotypes as a (samples x variants) shape data frame.
        Must be in string/nucleotide format. Must be of SNPs.
    basename: str
        Basename of file to save FASTA to; no directory or extension.
    wrap_length: int, optional
        How many nucleotides to wrap the alignment after.
        If not passed, then no wrapping will occur.

    Returns
    -------
    fig: plt.Figure
        The completed pseudo-alignment.
    """

    file = io_helpers.file_path('results', 'alignments', f'{basename}.fasta')
    
    # FASTA format: https://www.ncbi.nlm.nih.gov/genbank/fastaformat/
    with open(file, mode='w') as fasta:
        for sample in geno.index:
            # sequence ID is the sample's name
            fasta.write(f'>{sample}\n')
            fasta.write(geno.loc[sample].str.cat())
            fasta.write('\n')
    
    mv = pymsaviz.MsaViz(file, wrap_length=wrap_length)
    mv.set_plot_params(ticks_interval=None)

    return mv.plotfig()

def filter_maf(geno: Genotypes) -> Genotypes:
    """Remove variants from a `Genotypes` if they have MAF=0."""
    maf_0 = are_variants_maf_0(geno)
    print(f'Filtering out {maf_0.sum()} variants for MAF=0')
    return geno.loc[:, ~maf_0]

def are_variants_maf_0(geno: Genotypes) -> pd.Series: 
    """Determine whether each variant in a `Genotypes` has MAF=0."""
    return geno.nunique(axis='index', dropna=True) <= 1 

#endregion
#region ---- Private functions -----

def _vcf_to_df(vcf_file: str, chrom: str, has_info_scores: bool, 
               format: str) -> Tuple[Genotypes, np.ndarray]:
    """Read raw data from a VCF into a `Genotypes` data frame.

    No filters are applied.

    Parameters
    ----------
    vcf_file: str
        VCF file to read data from. One of the filename constants
        exported by this module.
    chrom: str
        Chromosome to read data for; either "Y" or "MT".
    has_info_scores: bool
        Whether this file has INFO score information to read.
    format: str
        Format to return data in; either "nucleotide" or "allele_num".
        See the `Genotypes` docstring for details.
    
    Returns
    -------
    geno: Genotypes
        Unfiltered genotypes as a (samples x variants) shape data frame.
    info_score: np.ndarray
        INFO scores corresponding to each variant in `geno`.
        Or, if `has_info_scores=False`, then `None`.
    """

    # ---- read into dictionary of NumPy arrays ----
    
    fields = ['calldata/GT', 'samples', 'variants/POS']
    if has_info_scores: fields.append('variants/INFO_SCORE')
    # allele sequences are only necessary for string format
    if format == 'nucleotide': fields += ['variants/REF', 'variants/ALT']
        
    # `read_vcf()` errors if tabix is not set to None, since none is available
    vcf = allel.read_vcf(vcf_file, alt_number=_NUM_ALT_SEQUENCES_READ, 
                         region=chrom, tabix=None, fields=fields)
    print(f'Read {vcf["variants/POS"].size} variants'
          f' across {vcf["samples"].size} samples')
    
    # ----- process into `Genotypes` -----
    
    # turn variants x samples x 2 array into variants x samples
    if format == 'nucleotide':
        alleles = np.concatenate((vcf['variants/REF'][:, np.newaxis], 
                                  vcf['variants/ALT']), axis=1).astype(str)
        gt = vcf['calldata/GT']

        if np.max(gt) >= alleles.shape[1]:
            raise ValueError(f'An allele numbered {np.max(gt)} appears in the'
                             f' genotypes; only {alleles.shape[1]} (0-indexed)'
                             f' alleles were read from {vcf_file}. Try'
                             ' increasing `_NUM_ALT_SEQUENCES_READ`.')
        gt = _gt_to_nucleotide_arr(gt, alleles)
    elif format == 'allele_num': gt = _gt_to_allele_num_arr(vcf['calldata/GT'])

    geno = pd.DataFrame(gt, index=vcf['variants/POS'], columns=vcf['samples'])
    info_scores = vcf['variants/INFO_SCORE'] if has_info_scores else None

    return geno.transpose(), info_scores

def _gt_to_nucleotide_arr(gt: np.ndarray, alleles: np.ndarray) -> np.ndarray:
    """Convert a 3D genotype array to a 2D nucleotide array.

    Each homozygous call (e.g. `0/0` or `1/1`), located in 
    `[row, col, :]`, is converted to its haploid nucleotide sequence.
    For example, if the reference allele is A and the first alternate
    is T, `0/0` becomes A and `1/1` becomes T.

    Heterozygous calls (e.g. `0/1`) are treated as missing. 
    Missing calls (`-1/-1`) in general are set to None.

    Parameters
    ----------
    gt: np.ndarray
        Allele numbers, shape (variants x samples x 2).
    alleles: np.ndarray
        Allele sequences, shape (variants x allele count).
    
    Returns
    -------
    gt: np.ndarray
        Haploid nucleotide genotypes, shape (variants x samples).
    """

    # missing genotypes are encoded as -1/-1
    is_missing = np.all(gt == np.array([-1, -1]), axis=-1)

    a1 = np.take_along_axis(alleles, gt[:, :, 0], 1)
    a2 = np.take_along_axis(alleles, gt[:, :, 1], 1)

    # both missing AND heterozygous positions are set to None/missing
    return np.where(np.logical_or(is_missing, a1 != a2), None, a1)

def _gt_to_allele_num_arr(gt: np.ndarray) -> np.ndarray:
    """Convert a 3D genotype array to a 2D allele number array.

    Each homozygous call (e.g. `0/0` or `1/1`), located in 
    `[row, col, :]`, is combined by the equation `a/b -> a * 10 + b`.
    For example, `0/1` becomes 11:
    - `[row, col, 0] = 0`
    - `[row, col, 1] = 1`
    - `0 * 10 + 1 = 11`

    Missing calls (`-1/-1`) end up as -11.

    Parameters
    ----------
    gt: np.ndarray
        Allele numbers, shape (variants x samples x 2).
    
    Returns
    -------
    gt: np.ndarray
        Numeric genotypes, shape (variants x samples).
    """

    return np.sum(gt * [10, 1], axis=-1)

def _subset_axis_by_label(geno: Genotypes, subset: ArrayLike, axis: str, 
                          why: str) -> Genotypes:
    """Subset a `Genotypes` by axis labels.

    Parameters
    ----------
    geno: Genotypes
        Genotypes as a (samples x variants) shape data frame.
    subset: ArrayLike
        Labels to subset the axis to.
    axis: str
        Axis to subset labels of; either "index" or "columns".
    why: str
        Reason/description of this subset. 
        Fills in "removed for not being {why}".
    
    Returns
    -------
    geno: Genotypes
        Input Genotypes with some rows or columns removed.
    """

    axis_num = 0 if axis == 'index' else 1

    old_size = geno.shape[axis_num]
    geno = geno.filter(subset, axis=axis)

    print(f'{old_size - geno.shape[axis_num]} {_GENOTYPES_AXIS_ITEM[axis]}s'
          f' removed for not being {why}')
    
    return geno

def _set_allele_num_missing(geno: Genotypes, allow_het: bool) -> Genotypes:
    """Set missing positions to `np.nan` in a numeric `Genotypes`.

    Parameters
    ----------
    geno: Genotypes
        Genotypes as a (samples x variants) shape data frame.
        Must be in numeric/allele_num format.
    allow_het: bool
        Whether to treat heterozygosity as missing.

    Returns
    -------
    geno: Genotypes
        Input Genotypes with some values set to `np.nan`.
    """
    
    if not allow_het:
        print('Setting heterozygous calls to missing')
        # homozygous `a/b` have `a = b`: `a * 10 + b = a * 11`, `% 11 == 0`
        geno.loc[:] = np.where(geno.to_numpy() % 11 != 0, -11, geno.to_numpy())
    
    print('Setting missing calls to np.nan')
    return geno.replace(-11, np.nan)

def _plot_hist(vals: ArrayLike, cutoff_val: float, cutoff_label: str, 
               title: str, x_label: str, y_label: str, plot: _Plot) -> None:
    """Plot a histogram with a cutoff line.

    Parameters
    ----------
    vals: ArrayLike
        Values to plot the distribution of.
    cutoff_val: float
        X-value to draw a vertical cutoff line at.
    cutoff_label: str
        What to label the cutoff with.
    title: str
        Title of the histogram.
    x_label: str
        X-axis label of the histogram.
    y_label: str
        Y-axis label of the histogram.
    plot: _Plot
        Subplot to put the histogram on.
    """
    
    tag, ax = plot
    
    ax.hist(vals)
    ax.axvline(x=cutoff_val, color='red', label=cutoff_label)
    ax.legend(framealpha=1, frameon=True)
    
    ax.set_title(title)
    ax.set_xlabel(x_label, fontsize=14)
    ax.set_ylabel(y_label, fontsize=14)
    io_helpers.write_panel_tag(tag, x=-0.2, y=1.1, ax=ax)
    
def _plot_maf(geno: Genotypes, chrom: str, plot: _Plot) -> None:
    """Plot of pie chart of MAF=0 vs. MAF>0.
    
    Parameters
    ----------
    geno: Genotypes
        Genotypes as a (samples x variants) shape data frame.
    chrom: str
        Chromosome these data are for; for plot title.
    plot: _Plot
        Subplot to put the pie chart on.    
    """

    tag, ax = plot
    maf_0 = are_variants_maf_0(geno)

    ax.pie([maf_0.sum(), (~maf_0).sum()], labels=['MAF = 0', 'MAF > 0'])

    # title & tag have to be shifted up to align with histograms
    ax.text(0.5, 1.08, f'MAF filter in {chrom}', transform=ax.transAxes, 
            fontsize=13, va='center', ha='center')
    io_helpers.write_panel_tag(tag, x=-0.2, y=1.1, ax=ax)
    
def _filter_miss_by_axis(geno: Genotypes, min_nonmiss_percent: float, 
                         by_axis: str, chrom: str, plot: _Plot) -> Genotypes:
    """Filter a `Genotypes` axis by missingness.

    The filter applied is a greater-than-or-equal-to cutoff.

    Parameters
    ----------
    geno: Genotypes
        Genotypes as a (samples x variants) shape data frame.
    min_nonmiss_percent: float
        Minimum genotyping rate/nonmissing call percent to allow.
    by_axis: str
        Axis to filter by; either "index" or "columns".
    chrom: str
        Chromosome these data are for; for plot title.
    plot: _Plot
        Subplot to put the histogram on.

    Returns
    -------
    geno: Genotypes
        Input Genotypes with some rows or columns removed.
    """

    # missingness is by/for `by_axis`, and thus across the other axis
    across_axis = 'index' if by_axis == 'columns' else 'columns'
    nonmiss_percent = 100 * (1 - geno.isnull().mean(axis=across_axis))

    if plot is not None:
        item = _GENOTYPES_AXIS_ITEM[by_axis]
        # plots assume they have SNPs, not generic variants
        if item == 'variant': item = 'SNP'

        _plot_hist(nonmiss_percent, cutoff_val=min_nonmiss_percent, 
                   cutoff_label=f'{100 - min_nonmiss_percent}% missingness',
                   title=f'Nonmissing SNPs, by {item} in {chrom}', 
                   x_label='% of nonmissing calls', y_label=f'# {item}s', 
                   plot=plot)
    return _filter_axis_by_cutoff(geno, vals=nonmiss_percent, 
                                  min_val=min_nonmiss_percent, axis=by_axis, 
                                  statistic='% nonmissing calls')

def _filter_axis_by_cutoff(geno: Genotypes, vals: np.ndarray, min_val: float, 
                           axis: str, statistic: str) -> Genotypes:
    """Filter a `Genotypes` axis by a statistic and threshold.

    The filter applied is a greater-than-or-equal-to cutoff.

    Parameters
    ----------
    geno: Genotypes
        Genotypes as a (samples x variants) shape data frame.
    vals: np.ndarray
        Values for a statistic corresponding to each row/column.
    min_val: float
        Minimum value which passes the cutoff.
    axis: str
        Axis to filter by; either "index" or "columns".
    statistic: str
        Name of the statistic used for filtration.

    Returns
    -------
    geno: Genotypes
        Input Genotypes with some rows or columns removed.
    """

    passed = vals >= min_val
    print(f'Filtering out {sum(~passed)} {_GENOTYPES_AXIS_ITEM[axis]}s'
          f' for {statistic} < {min_val}')
    return geno.loc[passed] if axis == 'index' else geno.loc[:, passed] 

#endregion