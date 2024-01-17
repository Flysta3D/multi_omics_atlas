from typing import Optional, Union
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.sparse import issparse
from anndata import AnnData

try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal


def basic_stats(
    adata: AnnData,
    mito_label: Optional[str] = None,
    copy: bool = False,
):
    adata = adata.copy() if copy else adata
    adata.var_names_make_unique()

    # The number of genes expressed in the count matrix.
    # adata.obs["nGenes"]: deteced genes per cell.
    # adata.obs["nCounts"]: total gene counts per cell.
    adata.obs["nGenes"], adata.obs["nCounts"] = np.array(
        (adata.X > 0).sum(1)
    ), np.array((adata.X).sum(1))

    # The number of cells counts in the count matrix.
    adata.var["nCells"], adata.var["nCounts"] = np.array(
        (adata.X > 0).sum(0).T
    ), np.array((adata.X).sum(0).T)

    if adata.var_names.inferred_type == "bytes":
        adata.var_names = adata.var_names.astype("str")
    mito_label = "MT-" if mito_label is None else mito_label.upper()
    mito_genes = adata.var_names.str.upper().str.startswith(mito_label)

    # The percentage of counts in mitochondrial genes.
    if sum(mito_genes) > 0:
        # The percentage of counts in mitochondrial genes.
        adata.obs["pMito"] = np.array(
            adata.X[:, mito_genes].sum(1).reshape((-1, 1))
            / adata.obs["nCounts"].values.reshape((-1, 1))
        )
        # The number of mitochondrial genes.
        if issparse(adata.X):
            adata.obs["nMito"] = adata.X[:, mito_genes].getnnz(axis=1)
        else:
            adata.obs["nMito"] = np.count_nonzero(adata.X[:, mito_genes], axis=1)

    else:
        print("No mitochondria genes detected.")

    return adata if copy else None


def save_fig(
    path=None,
    prefix=None,
    dpi=None,
    ext="pdf",
    transparent=True,
    close=True,
    verbose=True,
):
    """Save a figure from pyplot.
    code adapated from http://www.jesshamrick.com/2012/09/03/saving-figures-from-pyplot/
    Parameters
    ----------
         path: `string`
            The path (and filename, without the extension) to save_fig the
            figure to.
        prefix: `str` or `None`
            The prefix added to the figure name. This will be automatically set
            accordingly to the plotting function used.
        dpi: [ None | scalar > 0 | 'figure' ]
            The resolution in dots per inch. If None, defaults to rcParams["savefig.dpi"].
            If 'figure', uses the figure's dpi value.
        ext: `string` (default='pdf')
            The file extension. This must be supported by the active
            matplotlib backend (see matplotlib.backends module).  Most
            backends support 'png', 'pdf', 'ps', 'eps', and 'svg'.
        close: `boolean` (default=True)
            Whether to close the figure after saving.  If you want to save_fig
            the figure multiple times (e.g., to multiple formats), you
            should NOT close it in between saves or you will have to
            re-plot it.
        verbose: boolean (default=True)
            Whether to print information about when and where the image
            has been saved.
    """
    import matplotlib.pyplot as plt

    prefix = os.path.normpath(prefix)
    if path is None:
        path = os.getcwd() + "/"

    # Extract the directory and filename from the given path
    directory = os.path.split(path)[0]
    filename = os.path.split(path)[1]
    if directory == "":
        directory = "."
    if filename == "":
        filename = "dyn_savefig"

    # If the directory does not exist, create it
    if not os.path.exists(directory):
        os.makedirs(directory)

    # The final path to save_fig to
    savepath = (
        os.path.join(directory, filename + "." + ext)
        if prefix is None
        else os.path.join(directory, prefix + "_" + filename + "." + ext)
    )

    if verbose:
        print(f"Saving figure to {savepath}...")

    # Actually save the figure
    plt.savefig(
        savepath,
        dpi=300 if dpi is None else dpi,
        transparent=transparent,
        format=ext,
        bbox_inches="tight",
    )

    # Close it
    if close:
        plt.close()

    if verbose:
        print("Done")


def basic_stats_multi(
    adatas: Union[AnnData, list],
    slice_key: Optional[str] = None,
    inner: Literal["box", "quartiles", "point", None] = None,
    palette: str = "rainbow",
    x_rotation: int = -30,
    orient: Literal["row", "col"] = "row",
    save_show_or_return: Literal["show", "save", "return"] = "show",
    save_kwargs: Optional[dict] = None,
):
    """
    Plot the basic statics (area, nGenes, nCounts and pMito) of each category of adata.

    Args:
        adatas: An Anndata object or a list of Anndata object.
        slice_key: Which slice to facet the data into subplots.
        inner: Representation of the data points in the violin interior.
               * `inner` == `'box'`: draw a miniature boxplot;
               * `inner` == `'quartiles'`: draw the quartiles of the distribution;
               * `inner` == `'point'` or `'stick'`: show each underlying datapoint;
               * `inner` == None: draw unadorned violins.
        palette: Colors to use for the different levels of the hue variable.
        x_rotation: The rotation angle of the x-axis labels.
        orient: The direction of the drawing.
                * `'row'`: means to draw multiple rows;
                * `'col'`: means to draw multiple columns.
        save_show_or_return: Whether to save, show or return the figure.
        save_kwargs: A dictionary that will be pass to the save_fig function.
                     By default, it is an empty dictionary and the save_fig function will use the {"path": None,
                     "prefix": 'basic_stats', "dpi": 300, "ext": 'pdf', "transparent": True, "close": True, "verbose":
                     True} as its parameters. Otherwise, you can provide a dictionary that properly modify those keys
                     according to your needs.

    Returns:
        A violin plot that shows the fraction of each category, produced by seaborn.
    """
    # The final data used for plotting.
    data_col = ["area", "nGenes", "nCounts", "pMito"]
    plot_data = pd.DataFrame()

    adatas = adatas if isinstance(adatas, list) else [adatas]
    for i, adata in enumerate(adatas):
        if len(adata.obs.columns.intersection(["area", "nGenes", "nCounts", "pMito"])) != 4:
            raise ValueError(
                "Basic statics (area, nGenes, nCounts and pMito) do not exist in adata.obs."
            )
        if slice_key is not None and slice_key not in adata.obs.columns:
            raise ValueError("`slice_key` do not exist in adata.obs.")
        sub_data = adata.obs[data_col]
        sub_data["slices"] = (
            np.array([f"slice{i}"] * len(sub_data.index))
            if slice_key is None
            else adata.obs[slice_key]
        )
        plot_data = pd.concat([plot_data, sub_data], axis=0)

    res = plot_data.melt(value_vars=["area", "nGenes", "nCounts", "pMito"], id_vars=["slices"])
    # Draw a combination of boxplot, scatterplot and kernel density estimate.
    figsize_y = 4
    figsize_x = len(plot_data["slices"].unique().tolist()) * 1.2
    FacetGrid_kws = dict(
        sharex=False,
        sharey=False,
        hue="variable",
        height=figsize_y,
        aspect=figsize_x / figsize_y,
    )

    if orient == "row":
        g = sns.FacetGrid(res, row="variable", **FacetGrid_kws)
    else:
        g = sns.FacetGrid(res, col="variable", **FacetGrid_kws)

    violin_kws = dict(
        order=np.sort(res["slices"].unique()), width=0.8, palette=palette, saturation=0.8
    )
    strip_kws = dict(
        order=np.sort(res["slices"].unique()),
        color="whitesmoke",
        edgecolor="gray",
        size=1,
        jitter=True,
    )
    if inner == "point":
        g.map_dataframe(sns.violinplot, x="slices", y="value", inner=None, **violin_kws)
        g.map_dataframe(sns.stripplot, x="slices", y="value", **strip_kws)
    else:
        g.map_dataframe(sns.violinplot, x="slices", y="value", inner=inner, **violin_kws)

    # Set Labels
    for ax, ylabel in zip(g.axes.flat, data_col):
        plt.setp(ax.texts, text="")
        ax.set_ylabel(ylabel)

    g.set_titles("")
    g.set_xlabels("")
    g.set_xticklabels(rotation=x_rotation)
    g.set(ylim=(0, None))

    if save_show_or_return == "save":
        s_kwargs = {
            "path": None,
            "prefix": "basic_stats",
            "dpi": 300,
            "ext": "pdf",
            "transparent": True,
            "close": True,
            "verbose": True,
        }
        s_kwargs.update(
            (k, save_kwargs[k]) for k in s_kwargs.keys() & save_kwargs.keys()
        )
        save_fig(**s_kwargs)
    elif save_show_or_return == "show":
        plt.tight_layout()
        plt.show()
    elif save_show_or_return == "return":
        return g
