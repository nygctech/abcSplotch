import os
import glob
import h5py
import pickle
import numpy as np
import matplotlib.pyplot as plt
import rpy2.robjects as ro

from pathlib import Path
from rpy2.robjects.packages import importr

plt.rcParams['svg.fonttype'] = 'none'


def make_label_colors(labels, cmap_name="gist_ncar"):
    """
    Generate visually distinct colors for categorical labels.

    Parameters
    ----------
    labels : list
        List of categorical labels.
    cmap_name : str
        Matplotlib colormap name.

    Returns
    -------
    dict
        Mapping: label -> RGBA color
    """

    unique_labels = list(dict.fromkeys(labels))  # preserve order
    n = len(unique_labels)

    cmap = plt.get_cmap(cmap_name)

    # Evenly sample the colormap
    colors = [cmap(i / n) for i in range(n)]

    return dict(zip(unique_labels, colors))

# Scatterplot of cell-level annotations for a given tissue
def plot_annotations(tissue_name, splotch_inputs, figsize=(10,10)):
    '''
    Parameters:
    ----------
    tissue_name: str
        name of tissue to plot
    splotch_inputs: path
        path to directory containing abcSplotch input files (output by generate_splotch_inputs)

    Returns:
    -------
    ax: axis
    '''
    if not os.path.exists(os.path.join(splotch_inputs, 'information.p')):
        raise ValueError('No "information.p" file found in specified splotch inputs directory')
    sinfo = pickle.load(open(os.path.join(splotch_inputs, 'information.p'), 'rb'))

    # Example input file from which to get annotation data
    infiles = glob.glob(os.path.join(splotch_inputs, "*", 'data_*.R'))
    if len(infiles) == 0:
        raise ValueError('No abcSplotch input files found in "%s/*/data_*.R"' % splotch_inputs)
    
    rstan = importr("rstan")
    dat = rstan.read_rdump(infiles[0])
    annots = np.array(dat.rx2('D'))

    # Find spots under specified tissue
    tissues = np.array([Path(x[0]).stem.split('.')[0] for x in sinfo['filenames_and_coordinates']])
    coords = np.array([tuple(map(float, x[1].split('_'))) for x in sinfo['filenames_and_coordinates']])
    if tissue_name not in tissues:
        raise ValueError('Tissue "%s" not found among:' % tissue_name, np.unique(tissues))
    t_idx = tissues == tissue_name

    # Distinct colors for all labels in dataset
    cdict = make_label_colors(sinfo['annotation_mapping'])

    fig, ax = plt.subplots(1, figsize=figsize)
    ax.axis('off')
    ax.invert_xaxis()

    for a, lbl in enumerate(sinfo['annotation_mapping']):
        a_idx = np.logical_and(t_idx, annots==a+1)
        if np.sum(a_idx) > 0:
            ax.scatter(coords[a_idx,0], coords[a_idx,1], s=2, label=lbl, c=cdict[lbl])
    
    # Make legend with larger color points for visual clarity
    ax.legend(
        markerscale=5,
        scatterpoints=1
    )
    
    return ax
    
# Scatterplot of raw counts for a given gene on a given tissue
def plot_rawcounts(gene_name, tissue_name, splotch_inputs, figsize=(10,10), vmin=None, vmax=None):
    '''
    Parameters:
    ----------
    gene_name: str
        name of gene to plot
    tissue_name: str
        name of tissue to plot
    splotch_inputs: path
        path to directory containing abcSplotch input files (output by generate_splotch_inputs)

    Returns:
    -------
    ax: axis
    '''
    if not os.path.exists(os.path.join(splotch_inputs, 'information.p')):
        raise ValueError('No "information.p" file found in specified splotch inputs directory')
    sinfo = pickle.load(open(os.path.join(splotch_inputs, 'information.p'), 'rb'))

    if gene_name not in sinfo['genes']:
        raise ValueError('Gene "%s" not found among splotch-indexed genes' % gene_name)
    gid = list(sinfo['genes']).index(gene_name)+1

    # Read input file for specified gene
    infile = os.path.join(splotch_inputs, str(gid//100), 'data_%d.R' % gid)
    if not os.path.exists(infile):
        raise ValueError('%s (gid=%d) not found among inputs (%s)' % (gene_name, gid, splotch_inputs))

    rstan = importr("rstan")
    dat = rstan.read_rdump(infile)
    counts = np.array(dat.rx2('counts'))

    if vmin is None:
        vmin = counts.min()
    if vmax is None:
        vmax = np.percentile(counts, 99)

    # Find spots under specified tissue
    tissues = np.array([Path(x[0]).stem.split('.')[0] for x in sinfo['filenames_and_coordinates']])
    coords = np.array([tuple(map(float, x[1].split('_'))) for x in sinfo['filenames_and_coordinates']])
    if tissue_name not in tissues:
        raise ValueError('Tissue "%s" not found among:' % tissue_name, np.unique(tissues))
    t_idx = tissues == tissue_name

    fig, ax = plt.subplots(1, figsize=figsize)
    ax.axis('off')
    ax.invert_xaxis()

    pts = ax.scatter(coords[t_idx,0], coords[t_idx,1], c=counts[t_idx], s=2, vmin=vmin, vmax=vmax)
    plt.colorbar(pts, label='Counts')
    
    return ax

# Scatterplot of cell-level lambdas for a given gene on a given tissue
def plot_lambdas(gene_name, tissue_name, splotch_inputs, splotch_outputs, figsize=(10,10), vmin=None, vmax=None):
    '''
    Parameters:
    ----------
    gene_name: str
        name of gene to plot
    tissue_name: str
        name of tissue to plot
    splotch_inputs: path
        path to directory containing abcSplotch input files (output by generate_splotch_inputs)
    splotch_outputs: path
        path to direcotry containing abcSplotch output files (output by abcsplotch.r)

    Returns:
    -------
    ax: axis
    '''
    if not os.path.exists(os.path.join(splotch_inputs, 'information.p')):
        raise ValueError('No "information.p" file found in specified splotch inputs directory')
    sinfo = pickle.load(open(os.path.join(splotch_inputs, 'information.p'), 'rb'))

    if gene_name not in sinfo['genes']:
        raise ValueError('Gene "%s" not found among splotch-indexed genes' % gene_name)
    gid = list(sinfo['genes']).index(gene_name)+1

    # Read output file for specified gene
    outfile = os.path.join(splotch_outputs, str(gid//100), 'results_%d.R' % gid)
    if not os.path.exists(outfile):
        raise ValueError('%s (gid=%d) not found among outputs (%s)' % (gene_name, gid, splotch_outputs))
    readRDS = ro.r["readRDS"]
    dat = readRDS(outfile)
    
    lambdas = dat.rx2("lambda").rx2("mean")
    lambdas = np.array(lambdas)
    
    if vmin is None:
        vmin = lambdas.min()
    if vmax is None:
        vmax = np.percentile(lambdas, 99)

    # Find spots under specified tissue
    tissues = np.array([Path(x[0]).stem.split('.')[0] for x in sinfo['filenames_and_coordinates']])
    coords = np.array([tuple(map(float, x[1].split('_'))) for x in sinfo['filenames_and_coordinates']])
    if tissue_name not in tissues:
        raise ValueError('Tissue "%s" not found among:' % tissue_name, np.unique(tissues))
    t_idx = tissues == tissue_name

    fig, ax = plt.subplots(1, figsize=figsize)
    ax.axis('off')
    ax.invert_xaxis()

    pts = ax.scatter(coords[t_idx,0], coords[t_idx,1], c=lambdas[t_idx], s=2, vmin=vmin, vmax=vmax)
    plt.colorbar(pts, label='lambda (inferred counts)')
    
    return ax