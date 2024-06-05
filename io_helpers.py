"""Helper for file I/O.

Handles file reading and plot saving from a single place.

< CONSTANTS >

- `COLORBLIND_PALETTE`: greys and then 4 colorblind-friendly colors

< FUNCTIONS >

All public functions have some basic imput validation when necessary.
Many also have default or optional arguments; see docstrings.

- `load_sample_info_with_haplogroups()`: load sample info+haplogroups.
- `load_haplogroups()`: load pre-existing haplotype group assignments.
- `read_csv()`: read a CSV file; `file_path()` syntax for the path.
- `read_col()`: read one column from CSV to make a `pd.Series`.
- `save_plot()`: save a plot as an image file.
- `file_path()`: build a filepath within the working directory.
- `write_panel_tag()`: add a letter tag to a figure panel.
- `roman_format()`: format a word as non-italic for Matplotlib.

All but the last of these are essentially wrappers:
`load_sample_info_with_haplogroups()` and `load_haplogroups()` for
`read_csv()` to handle the common case of the sample information and
Y/MT haplogroup files, `read_csv()` for the corresponding pandas
function, `save_plot()` runs Matplotlib's `figure.savefig()` to the
`plots` directory, `file_path()` is `os.path.join()`, and
`write_panel_tag()` is Matplotlib's `ax.text()`.
"""

#region ---- Imports ----

import os # use robust file paths
import sys # access Python version number

import matplotlib # access version number
import matplotlib.pyplot as plt # make simple visualizations
import pandas as pd # read CSVs

print(f'''
This I/O helper file uses:
Python version {sys.version},
Matplotlib version {matplotlib.__version__},
pandas version {pd.__version__}
''')

#endregion
#region ---- Constants -----

# move working directory up from the 'code' subfolder
_WD = os.path.dirname(os.getcwd())

# https://www.nature.com/articles/nmeth.1618
COLORBLIND_PALETTE = ['#444444', '#888888', '#E69F00', 
                      '#56B4E9', '#009E73', '#F0E442']

#endregion
#region ---- Public functions ----

def load_sample_info_and_haplogroups(add_female_to_y: bool) -> pd.DataFrame:
    """Load sample information table.
    
    Adds columns for pre-existing haplogroups.

    Parameters
    ----------
    add_female_to_y: bool
        Whether to assign females a `Y_group` of "female".
    
    Returns
    -------
    sample_info: pd.DataFrame
        Sample table, with RFID, sex, library prep method, and DOB.
    """

    sample_info = read_csv('data', 'sample_info.csv', index_col='rfid')

    sample_info = sample_info.join(load_haplogroups('Y'))
    if add_female_to_y:
        sample_info.loc[sample_info['sex'] == 'F', 'Y_group'] = 'female'
    sample_info = sample_info.join(load_haplogroups('MT'))

    return sample_info

def load_haplogroups(chrom: str, keep_founders: bool = False) -> pd.Series:
    """Load Y/MT haplotype groups.
    
    Parameters
    ----------
    chrom: str
        Chromosome for file to read; either "Y" or "MT".
    keep_founders: bool, default=False
        Whether to retain founder annotations.
        If False, then e.g. Y1 (ACI) -> Y1 conversion is done.
    
    Returns
    -------
    groups: pd.Series
        Haplotype groups (index is RFID).
    """

    if chrom != 'Y' and chrom != 'MT':
        raise ValueError('`chrom` must be either "Y" or "MT":'
                         f' `{chrom}` is not a valid value')
    
    groups = read_col('results', 'groups', f'{chrom}_groups.csv', 
                      colname=f'{chrom}_group', index_col='rfid')
    if not keep_founders: groups = groups.str.slice(0, 2 if chrom == 'Y' else 3)
    
    return groups

def read_csv(*args, **kwargs) -> pd.DataFrame:
    """Read a CSV file using `file_path()` syntax."""
    return pd.read_csv(file_path(*args), **kwargs)

def read_col(*args, colname: str, index_col: str, **kwargs) -> pd.Series:
    """Read a CSV column."""
    return read_csv(*args, index_col=index_col, **kwargs)[colname]

def save_plot(filename: str, extension: str = 'svg', 
              figure: plt.Figure = None) -> None:
    """Save a figure (defaults to current) into `results/plots/`."""
    if figure is None: figure = plt.gcf()
    file = file_path('results', 'plots', f'{filename}.{extension}')
    figure.savefig(file, transparent=True)
 
def file_path(*args) -> str: return os.path.join(_WD, *args)

def write_panel_tag(tag: str, x: float, y: float, fontsize: int = 16,
                    ax: plt.Axes = None) -> None:
    """Write a letter tag on a (sub)plot as a panel label.

    Parameters
    ----------
    tag: str
        Letter to use as a tag.
    x: float
        X-coordinate within the Axes to write the tag at.
    y: float
        Y-coordinate within the Axes to write the tag at.
    fontsize: int, default=16
        Font size to write the tag with.
    ax: plt.Axes, optional
        Plot to write the tag on. Defaults to the current axis.
    """

    if len(tag) != 1: 
        raise ValueError('`tag` must be a single letter'
                         f': `{tag}` is not a valid value')
    if ax is None: ax = plt.gca()

    ax.text(x, y, tag, transform=ax.transAxes, fontsize=fontsize, 
            fontweight='bold', va='top', ha='left')

def roman_format(text: str) -> str:
    """Set roman format (blocks italics)."""
    return ' '.join(r'$\rm{' + word + r'}$' for word in text.split())

#endregion