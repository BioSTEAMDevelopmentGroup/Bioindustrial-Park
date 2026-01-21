# -*- coding: utf-8 -*-
"""
PreFerS Plot Utilities
======================

Output management utilities for consistent file saving and directory structure.

@author: Dr. Ouwen Peng
@institute: PreFerS - Centre for Precision Fermentation and Sustainability
"""

import os
from datetime import datetime
import matplotlib.pyplot as plt


# =============================================================================
# Output Directory Management
# =============================================================================

def get_analysis_dirs(script_path=None, config="", timestamp=None, create=True):
    """
    Create and return standardized analysis output directories.
    
    Creates the standard structure:
        analyses/
        ├── data/      # For Excel, CSV, pickle files
        └── figure/    # For PNG, PDF, SVG plots
    
    Parameters
    ----------
    script_path : str, optional
        Path to the calling script. If None, uses current working directory.
    config : str
        Configuration name to include in folder name (e.g., 'config1').
    timestamp : str, optional
        Custom timestamp. If None, generates current timestamp (YYYYMMDD_HHMM).
    create : bool
        Whether to create the directories if they don't exist.
    
    Returns
    -------
    dict
        Dictionary with 'base', 'data', and 'figure' directory paths.
    
    Examples
    --------
    >>> dirs = get_analysis_dirs(__file__, config='config1')
    >>> print(dirs['data'])   # .../analyses/results_config1_20260121_0030/data/
    >>> print(dirs['figure']) # .../analyses/results_config1_20260121_0030/figure/
    """
    if script_path is None:
        base_dir = os.getcwd()
    else:
        base_dir = os.path.dirname(os.path.abspath(script_path))
    
    # If script is already in an 'analyses' folder, use that instead of creating nested one
    if os.path.basename(base_dir).lower() == 'analyses':
        analyses_parent = base_dir
    else:
        analyses_parent = os.path.join(base_dir, "analyses")
    
    if timestamp is None:
        timestamp = datetime.now().strftime('%Y%m%d_%H%M')
    
    # Build folder name
    if config:
        folder_name = f"results_{config}_{timestamp}"
    else:
        folder_name = f"results_{timestamp}"
    
    # Create analyses base directory
    analyses_base = os.path.join(analyses_parent, folder_name)
    data_dir = os.path.join(analyses_base, "data")
    figure_dir = os.path.join(analyses_base, "figure")
    
    if create:
        os.makedirs(data_dir, exist_ok=True)
        os.makedirs(figure_dir, exist_ok=True)
    
    return {
        'base': analyses_base,
        'data': data_dir,
        'figure': figure_dir,
        'timestamp': timestamp,
        'config': config
    }


def get_output_dir(script_path=None, config="", timestamp=None, 
                   create=True, subfolder=None):
    """
    Create and return a timestamped output directory.
    
    Follows the convention: results_{config}_{YYYYMMDD_HHMM}/
    
    Parameters
    ----------
    script_path : str, optional
        Path to the calling script. If None, uses current working directory.
    config : str
        Configuration name to include in folder name.
    timestamp : str, optional
        Custom timestamp. If None, generates current timestamp.
    create : bool
        Whether to create the directory if it doesn't exist.
    subfolder : str, optional
        Additional subfolder (e.g., 'figures', 'tables').
    
    Returns
    -------
    str
        Absolute path to output directory.
    
    Examples
    --------
    >>> output_dir = get_output_dir(__file__, config='config1')
    >>> # Returns: /path/to/script/results_config1_20260121_0030/
    """
    if script_path is None:
        base_dir = os.getcwd()
    else:
        base_dir = os.path.dirname(os.path.abspath(script_path))
    
    if timestamp is None:
        timestamp = datetime.now().strftime('%Y%m%d_%H%M')
    
    if config:
        folder_name = f"results_{config}_{timestamp}"
    else:
        folder_name = f"results_{timestamp}"
    
    output_dir = os.path.join(base_dir, folder_name)
    
    if subfolder:
        output_dir = os.path.join(output_dir, subfolder)
    
    if create:
        os.makedirs(output_dir, exist_ok=True)
    
    return output_dir


def get_figures_dir(script_path=None, config="", timestamp=None):
    """
    Get the figures subdirectory within the output directory.
    
    Parameters
    ----------
    script_path : str, optional
        Path to the calling script.
    config : str
        Configuration name.
    timestamp : str, optional
        Custom timestamp.
    
    Returns
    -------
    str
        Path to figures directory.
    """
    return get_output_dir(script_path, config, timestamp, 
                          create=True, subfolder='figures')


def get_tables_dir(script_path=None, config="", timestamp=None):
    """
    Get the tables subdirectory within the output directory.
    
    Parameters
    ----------
    script_path : str, optional
        Path to the calling script.
    config : str
        Configuration name.
    timestamp : str, optional
        Custom timestamp.
    
    Returns
    -------
    str
        Path to tables directory.
    """
    return get_output_dir(script_path, config, timestamp, 
                          create=True, subfolder='tables')


# =============================================================================
# Figure Saving
# =============================================================================

def save_figure(fig, name, output_dir, formats=('png',), dpi=300, 
                close=True, verbose=True):
    """
    Save a figure with consistent settings.
    
    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure to save.
    name : str
        Base filename (without extension).
    output_dir : str
        Directory to save to.
    formats : tuple of str
        File formats to save (e.g., ('png', 'pdf')).
    dpi : int
        Resolution for raster formats.
    close : bool
        Whether to close the figure after saving.
    verbose : bool
        Whether to print save confirmation.
    
    Returns
    -------
    list of str
        Paths to saved files.
    
    Examples
    --------
    >>> fig, ax = plt.subplots()
    >>> ax.plot([1, 2, 3], [1, 4, 9])
    >>> save_figure(fig, 'my_plot', output_dir, formats=('png', 'pdf'))
    """
    os.makedirs(output_dir, exist_ok=True)
    
    saved_paths = []
    for fmt in formats:
        filepath = os.path.join(output_dir, f"{name}.{fmt}")
        fig.savefig(filepath, dpi=dpi, bbox_inches='tight', 
                    facecolor='white', edgecolor='none')
        saved_paths.append(filepath)
        
        if verbose:
            print(f"  ✓ Saved: {filepath}")
    
    if close:
        plt.close(fig)
    
    return saved_paths


def save_table(df, name, output_dir, formats=('xlsx',), verbose=True):
    """
    Save a DataFrame as table file.
    
    Parameters
    ----------
    df : pd.DataFrame
        Data to save.
    name : str
        Base filename (without extension).
    output_dir : str
        Directory to save to.
    formats : tuple of str
        File formats ('xlsx', 'csv').
    verbose : bool
        Whether to print save confirmation.
    
    Returns
    -------
    list of str
        Paths to saved files.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    saved_paths = []
    for fmt in formats:
        filepath = os.path.join(output_dir, f"{name}.{fmt}")
        
        if fmt == 'xlsx':
            df.to_excel(filepath)
        elif fmt == 'csv':
            df.to_csv(filepath)
        else:
            df.to_csv(filepath.replace(f'.{fmt}', '.csv'))
        
        saved_paths.append(filepath)
        
        if verbose:
            print(f"  ✓ Saved: {filepath}")
    
    return saved_paths


# =============================================================================
# Timestamp Utilities
# =============================================================================

def get_timestamp(format_str='%Y.%m.%d-%H.%M'):
    """
    Get current timestamp string for file naming.
    
    Parameters
    ----------
    format_str : str
        strftime format string.
    
    Returns
    -------
    str
        Formatted timestamp.
    """
    return datetime.now().strftime(format_str)


def get_timestamp_compact():
    """Get compact timestamp: YYYYMMDD_HHMM"""
    return datetime.now().strftime('%Y%m%d_%H%M')
