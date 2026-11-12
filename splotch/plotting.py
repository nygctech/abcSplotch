import os
import glob
import h5py
import pickle
import numbers
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

def plot_tissue_comparison(
    gene_name,
    tissue_name,
    splotch_inputs,
    splotch_outputs,
    figsize=(18, 6),
    marker_size=2,
    counts_vmin=None,
    counts_vmax=None,
    lambda_vmin=None,
    lambda_vmax=None,
    cmap="viridis",
    legend=True,
):
    """
    Plot tissue annotations, raw counts, and inferred lambdas
    side by side with identical spatial aspect ratios.

    Returns
    -------
    fig : matplotlib.figure.Figure
    axes : ndarray of matplotlib.axes.Axes
    """
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    # Load metadata
    info_path = os.path.join(splotch_inputs, "information.p")
    with open(info_path, "rb") as f:
        sinfo = pickle.load(f)

    genes = list(sinfo["genes"])
    if gene_name not in genes:
        raise ValueError(f"Gene '{gene_name}' not found")

    gid = genes.index(gene_name) + 1

    # Load gene-specific counts and annotations
    infile = os.path.join(
        splotch_inputs, str(gid // 100), f"data_{gid}.R"
    )
    assert os.path.exists(infile), f"Input file {infile} does not exist"
    
    rstan = importr("rstan")
    dat = rstan.read_rdump(infile)

    counts = np.asarray(dat.rx2("counts"))
    annots = np.asarray(dat.rx2("D"))

    # Load inferred lambdas
    outfile = os.path.join(
        splotch_outputs, str(gid // 100), f"results_{gid}.R"
    )
    assert os.path.exists(outfile), f"Output file {outfile} does not exist"
    
    res = ro.r["readRDS"](outfile)
    lambdas = np.asarray(res.rx2("lambda").rx2("mean"))

    # Extract tissue coordinates
    entries = sinfo["filenames_and_coordinates"]
    tissues = np.array([
        Path(x[0]).stem.split(".")[0] for x in entries
    ])
    coords = np.array([
        tuple(map(float, x[1].split("_"))) for x in entries
    ])

    mask = tissues == tissue_name
    if not np.any(mask):
        raise ValueError(f"Tissue '{tissue_name}' not found")

    xy = coords[mask]
    annots = annots[mask]
    counts = counts[mask]
    lambdas = lambdas[mask]

    # Common spatial limits, with small padding
    xmin, ymin = xy.min(axis=0)
    xmax, ymax = xy.max(axis=0)
    dx = max(xmax - xmin, 1)
    dy = max(ymax - ymin, 1)
    xlim = (xmin - 0.02 * dx, xmax + 0.02 * dx)
    ylim = (ymin - 0.02 * dy, ymax + 0.02 * dy)

    # Color scaling (matching original percentile convention)
    if counts_vmin is None:
        counts_vmin = counts.min()
    if counts_vmax is None:
        counts_vmax = np.percentile(counts, 99)
    if lambda_vmin is None:
        lambda_vmin = lambdas.min()
    if lambda_vmax is None:
        lambda_vmax = np.percentile(lambdas, 99)

    fig, axes = plt.subplots(
        1, 3, figsize=figsize,
        gridspec_kw={"wspace": 0.25}
    )

    # Panel 1: Annotations
    labels = sinfo["annotation_mapping"]
    colors = make_label_colors(labels)

    for i, label in enumerate(labels):
        idx = annots == i + 1
        if np.any(idx):
            axes[0].scatter(
                xy[idx, 0], xy[idx, 1],
                s=marker_size,
                c=[colors[label]],
                label=label,
                rasterized=True,
            )

    if legend:
        axes[0].legend(
            loc="upper left",
            bbox_to_anchor=(0, -0.05),
            ncol=2,
            markerscale=4,
            fontsize=8,
            frameon=False,
        )

    axes[0].set_title("Annotations")

    # Panels 2 and 3: Expression
    for ax, values, vmin, vmax, title, cbar_label in [
        (axes[1], counts, counts_vmin, counts_vmax,
         f"{gene_name}: Raw counts", "Counts"),
        (axes[2], lambdas, lambda_vmin, lambda_vmax,
         f"{gene_name}: Inferred lambda", "Lambda"),
    ]:
        pts = ax.scatter(
            xy[:, 0], xy[:, 1],
            c=values,
            s=marker_size,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            rasterized=True,
        )
        ax.set_title(title)

        # Separate colorbar axis preserves plot aspect ratio
        divider = make_axes_locatable(ax)
        cax = divider.append_axes(
            "right", size="4%", pad=0.08
        )
        fig.colorbar(pts, cax=cax, label=cbar_label)

    # Enforce identical spatial geometry
    for ax in axes:
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.invert_xaxis()
        ax.set_aspect("equal", adjustable="box")
        ax.axis("off")

    return fig, axes

def resolve_splotch_labels(values, mapping, kind):
    """
    Convert human-readable labels or 1-based integer indices
    to (index, display_label) pairs.

    mapping: ordered sequence of names, or dictionary
             mapping integer indices to names.
    """
    if isinstance(mapping, dict):
        if all(isinstance(k, numbers.Integral) for k in mapping):
            lookup = {int(k): str(v) for k, v in mapping.items()}
        else:
            raise ValueError(
                f"{kind} mapping must be an ordered sequence "
                "or an integer-to-name dictionary"
            )
    else:
        lookup = {
            i: str(name)
            for i, name in enumerate(mapping, start=1)
        }

    reverse = {}
    for idx, name in lookup.items():
        if name in reverse:
            raise ValueError(
                f"Duplicate {kind} name '{name}' in mapping"
            )
        reverse[name] = idx

    resolved = []

    for value in values:
        if isinstance(value, numbers.Integral):
            idx = int(value)
            if idx not in lookup:
                raise ValueError(
                    f"Unknown {kind} index {idx}. "
                    f"Valid indices: {list(lookup)}"
                )
        elif isinstance(value, str):
            if value not in reverse:
                raise ValueError(
                    f"Unknown {kind} name '{value}'. "
                    f"Available names: {list(reverse)}"
                )
            idx = reverse[value]
        else:
            raise TypeError(
                f"{kind} must be specified by name or integer index"
            )

        resolved.append((idx, lookup[idx]))

    return resolved

def plot_beta_kde(
    gene_name,
    splotch_inputs,
    splotch_outputs,
    level=1,
    conditions=None,
    annotations=None,
    comparison="conditions",
    scale="log",
    figsize=None,
    ax=None,
    cmap="tab10",
    alpha=0.9,
    linewidth=2,
    fill=False,
    show_zero=True,
    sharex=True,
    sharey=True,
):
    """
    Plot marginal posterior densities of effective beta coefficients.

    Conditions, annotations, and cell types may be specified using
    either human-readable names or 1-based integer indices.

    Parameters
    ----------
    gene_name : str
        Gene to visualize.
    splotch_inputs, splotch_outputs : str
        abcSplotch input and output directories.
    level : {1, 2, 3}
        Hierarchical level of effective beta coefficients.
    conditions : list of str or int
        Conditions at the requested hierarchical level.
    annotations : list
        Non-compositional:
            Annotation names or integer IDs, e.g.
            ["Muscle", "Connective Tissue"] or [1, 2].

        Compositional:
            (annotation, celltype) pairs, with either names
            or integer IDs, e.g.
            [("Muscle", "Myonuclei"), (2, 3)].
    comparison : {"conditions", "annotations"}
        "conditions": one subplot per annotation, curves by condition.
        "annotations": one subplot per condition, curves by annotation.
    scale : {"log", "expression"}
        "log": plot effective beta on its original log scale.
        "expression": plot exp(beta), transforming posterior
        densities using the appropriate Jacobian.
    figsize : tuple or None
        Figure size. Defaults to (6*n_panels, 4).
    ax : matplotlib Axes or None
        Optional existing axis for single-panel plots.
    cmap : str
        Matplotlib colormap.
    alpha : float
        Curve opacity.
    linewidth : float
        Curve line width.
    fill : bool
        Fill area beneath posterior density curves.
    show_zero : bool
        Show reference line at beta=0 (or exp(beta)=1).
    sharex, sharey : bool
        Share axes between subplots.

    Returns
    -------
    fig : matplotlib.figure.Figure
    axes : ndarray of matplotlib.axes.Axes
    """
    # --------------------------------------------------
    # Validate arguments
    # --------------------------------------------------
    if level not in (1, 2, 3):
        raise ValueError("level must be 1, 2, or 3")
    
    if comparison not in ("conditions", "annotations"):
        raise ValueError(
            "comparison must be 'conditions' or 'annotations'"
        )
    
    if scale not in ("log", "expression"):
        raise ValueError(
            "scale must be 'log' or 'expression'"
        )

    # --------------------------------------------------
    # Load gene index and results
    # --------------------------------------------------
    with open(
        os.path.join(splotch_inputs, "information.p"), "rb"
    ) as f:
        sinfo = pickle.load(f)

    genes = list(sinfo["genes"])

    if gene_name not in genes:
        raise ValueError(
            f"Gene '{gene_name}' not found in information.p"
        )

    gid = genes.index(gene_name) + 1

    outfile = os.path.join(
        splotch_outputs,
        str(gid // 100),
        f"results_{gid}.R",
    )

    if not os.path.isfile(outfile):
        raise FileNotFoundError(outfile)

    res = ro.r["readRDS"](outfile)
    beta = res.rx2("beta_effective")
    level_key = f"beta_l{level}"

    available_levels = list(
        ro.r["names"](beta.rx2("marginals"))
    )

    if level_key not in available_levels:
        raise ValueError(
            f"{level_key} not found. "
            f"Available levels: {available_levels}"
        )

    marginals = beta.rx2("marginals").rx2(level_key)
    metadata = beta.rx2("metadata").rx2(level_key)

    # --------------------------------------------------
    # Extract coefficient metadata
    # --------------------------------------------------
    columns = list(ro.r["names"](metadata))

    def get_column(name):
        return list(metadata.rx2(name))

    coef_names = [
        str(x) for x in get_column("coefficient")
    ]

    condition_col = (
        f"level_{level}_condition"
        if level > 1 else "condition"
    )

    condition_ids = [
        int(x) for x in get_column(condition_col)
    ]

    mroi_ids = [
        int(x) for x in get_column("mroi")
    ]

    compositional = "cell_type" in columns

    if compositional:
        cell_ids = [
            int(x) for x in get_column("cell_type")
        ]
    else:
        cell_ids = [None] * len(coef_names)

    # Map (condition, MROI, celltype) -> coefficient name
    coef_lookup = {
        (c, m, k): name
        for c, m, k, name in zip(
            condition_ids,
            mroi_ids,
            cell_ids,
            coef_names,
        )
    }

    # Map coefficient name -> posterior marginal
    marginal_names = list(ro.r["names"](marginals))

    marginal_lookup = {
        str(name): marginals.rx2(str(name))
        for name in marginal_names
    }

    # --------------------------------------------------
    # Resolve condition names for the requested level
    # --------------------------------------------------
    condition_mapping = sinfo["beta_mapping"][
        f"beta_level_{level}"
    ]
    
    if conditions is None:
        conditions = list(range(1, len(condition_mapping) + 1))
    else:
        conditions = list(conditions)
    
    resolved_conditions = resolve_splotch_labels(
        conditions,
        condition_mapping,
        f"level-{level} condition",
    )
    
    # --------------------------------------------------
    # Resolve annotation and cell-type names
    # --------------------------------------------------
    annotation_mapping = sinfo["annotation_mapping"]
    
    if annotations is None:
        if compositional:
            # All annotation-celltype combinations
            celltype_mapping = sinfo["celltype_mapping"]
    
            annotations = [
                (m, k)
                for m in range(1, len(annotation_mapping) + 1)
                for k in range(1, len(celltype_mapping) + 1)
            ]
        else:
            annotations = list(
                range(1, len(annotation_mapping) + 1)
            )
    else:
        annotations = list(annotations)
    
    if compositional:
        if not all(
            isinstance(a, (tuple, list)) and len(a) == 2
            for a in annotations
        ):
            raise ValueError(
                "Compositional annotations must be "
                "(annotation, celltype) pairs"
            )
    
        mroi_values = [a[0] for a in annotations]
        cell_values = [a[1] for a in annotations]
    
    else:
        if any(
            isinstance(a, (tuple, list))
            for a in annotations
        ):
            raise ValueError(
                "Non-compositional annotations must be "
                "individual annotation names or IDs"
            )
    
        mroi_values = annotations
    
    resolved_mrois = resolve_splotch_labels(
        mroi_values,
        annotation_mapping,
        "annotation",
    )
    
    if compositional:
        resolved_cells = resolve_splotch_labels(
            cell_values,
            sinfo["celltype_mapping"],
            "cell type",
        )
    
        parsed_annotations = [
            (
                (m_idx, k_idx),
                f"{m_name}, {k_name}",
            )
            for (m_idx, m_name), (k_idx, k_name)
            in zip(resolved_mrois, resolved_cells)
        ]
    
    else:
        parsed_annotations = [
            ((m_idx, None), m_name)
            for m_idx, m_name in resolved_mrois
        ]

    # --------------------------------------------------
    # Organize panels and comparison curves
    # --------------------------------------------------
    if comparison == "conditions":
        panels = parsed_annotations
        curves = resolved_conditions

        def make_key(panel, curve):
            (m, k), _ = panel
            c, _ = curve
            return (c, m, k)

        def panel_label(panel):
            return panel[1]

        def curve_label(curve):
            return curve[1]

    else:
        panels = resolved_conditions
        curves = parsed_annotations

        def make_key(panel, curve):
            c, _ = panel
            (m, k), _ = curve
            return (c, m, k)

        def panel_label(panel):
            return panel[1]

        def curve_label(curve):
            return curve[1]

    # --------------------------------------------------
    # Create figure and axes
    # --------------------------------------------------
    n_panels = len(panels)

    if ax is not None:
        if n_panels != 1:
            raise ValueError(
                "ax can only be supplied for single-panel plots"
            )

        fig = ax.figure
        axes = np.array([ax])

    else:
        if figsize is None:
            figsize = (6 * n_panels, 4)

        fig, axes = plt.subplots(
            1,
            n_panels,
            figsize=figsize,
            sharex=sharex,
            sharey=sharey,
            squeeze=False,
        )

        axes = axes.ravel()

    colors = plt.get_cmap(cmap)(
        np.linspace(0, 1, max(len(curves), 2))
    )

    # --------------------------------------------------
    # Plot posterior marginal densities
    # --------------------------------------------------
    for ax_i, panel in zip(axes, panels):
        n_plotted = 0

        for j, curve in enumerate(curves):
            key = make_key(panel, curve)
            coef = coef_lookup.get(key)

            if coef is None:
                print(f"Missing coefficient: {key}")
                continue

            if coef not in marginal_lookup:
                print(f"Missing marginal: {coef}")
                continue

            marginal = np.asarray(
                marginal_lookup[coef],
                dtype=float,
            )

            if (
                marginal.ndim != 2
                or marginal.shape[1] != 2
            ):
                raise ValueError(
                    f"Unexpected marginal shape for {coef}: "
                    f"{marginal.shape}"
                )

            x = marginal[:, 0]
            density = marginal[:, 1]

            # Transform beta -> exp(beta)
            if scale == "expression":
                x = np.exp(x)
                density = density / x

            order = np.argsort(x)
            x = x[order]
            density = density[order]

            color = colors[j]
            label = curve_label(curve)

            ax_i.plot(
                x,
                density,
                color=color,
                alpha=alpha,
                linewidth=linewidth,
                label=label,
            )

            if fill:
                ax_i.fill_between(
                    x,
                    density,
                    color=color,
                    alpha=0.12,
                )

            n_plotted += 1

        # Reference line
        if show_zero:
            ax_i.axvline(
                0 if scale == "log" else 1,
                color="gray",
                linestyle="--",
                linewidth=1,
                alpha=0.6,
                zorder=0,
            )

        ax_i.set_title(panel_label(panel))

        ax_i.set_xlabel(
            "Effective beta (log expression)"
            if scale == "log"
            else "exp(effective beta)"
        )

        ax_i.set_ylabel("Posterior density")

        if n_plotted:
            ax_i.legend(frameon=False)

        ax_i.spines["top"].set_visible(False)
        ax_i.spines["right"].set_visible(False)

    fig.suptitle(
        f"{gene_name} — Effective beta, hierarchy level {level}",
        fontsize=13,
    )

    if ax is None:
        fig.tight_layout()

    return fig, axes