import math
import operator
import warnings
from typing import Optional, Union

import matplotlib as mpl
import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.sparse import isspmatrix
from sklearn.preprocessing import minmax_scale

pal = [
    "#DC143C",
    "#0000FF",
    "#20B2AA",
    "#FFA500",
    "#9DBEB9",
    "#7CFC00",
    "#FFFF00",
    "#666666",
    "#DBF6E9",
    "#C9CBFF",
    "#00416D",
    "#40A8C4",
    "#E5D549",
    "#808000",
    "#A03C78",
    "#1B9E77",
    "#D95F02",
    "#7570B3",
    "#E7298A",
    "#66A61E",
    "#E6AB02",
    "#A6761D",
    "#D35D6E",
    "#EFB08C",
    "#F5F1DA",
    "#FF5722",
    "#ADE498",
    "#B1B493",
    "#35D0BA",
    "#864000",
    "#E4BAD4",
    "#00E8FF",
    "#FF00DD",
    "#B5FF00",
    "#FF7700",
    "#006DFF",
    "#00FF91",
    "#1AAB27",
    "#7F90F0",
    "#B774C5",
    "#A28706",
    "#5A2686",
    "#003FFF",
    "#6787E8",
    "#488D6B",
    "#FF00A0",
    "#E8DB53",
    "#75AA4A",
    "#1F3EEE",
    "#E368E7",
    "#F65145",
    "#0087FF",
    "#FFDD00",
    "#F8A76F",
]

# Convert sparse matrix to dense matrix.
to_dense_matrix = lambda X: np.array(X.todense()) if isspmatrix(X) else np.asarray(X)


def groups_cmap(
    groups: list,
    colormap: Union[str, list, dict] = "rainbow",
    mask_color: str = "gainsboro",
) -> dict:

    # Create a dictionary that stores all the groups and their colors.
    raw_groups = groups.copy()
    raw_groups.sort()

    gcdict = {}

    # Processing mask.
    if "mask" in groups:
        groups.remove("mask")
        mask_hex = mpl.colors.to_hex(mask_color, keep_alpha=True)
        gcdict["mask"] = mask_hex

    # Processing the rest of the groups
    if colormap == "auto":
        colormap = pal[: len(groups)]

    if (
        colormap != "auto"
        and isinstance(colormap, str)
        and not (colormap in list(mpl.colormaps))
    ):
        colormap = [colormap] * len(groups)

    if (
        colormap != "auto"
        and isinstance(colormap, str)
        and colormap in list(mpl.colormaps)
    ):
        lscmap = mpl.cm.get_cmap(colormap)
        colormap = [
            mpl.colors.to_hex(lscmap(i)) for i in np.linspace(0, 1, len(groups))
        ]

    if isinstance(colormap, list):
        for group, color in zip(groups, colormap):
            gcdict[group] = mpl.colors.to_hex(color, keep_alpha=True)
    elif isinstance(colormap, dict):
        for group in groups:
            gcdict[group] = mpl.colors.to_hex(colormap[group], keep_alpha=True)

    # Check if the gcdict contains the correct number of groups
    gcdict_keys = [key for key in gcdict.keys()]
    gcdict_keys.sort()
    if operator.eq(gcdict_keys, raw_groups):
        return gcdict
    else:
        raise ValueError("Wrong color and transparency settings for groups.")


def space_multi(
    adata: AnnData,
    groupby: Union[str, tuple] = None,
    spatial_key: str = "spatial",
    slice_key: Optional[str] = "slices",
    colormap: Union[str, list, dict] = "auto",
    alphamap: Union[float, str] = "auto",
    mask: Union[str, int, float, list] = None,
    mask_color: str = "gainsboro",
    mask_alpha: float = 1.0,
    point_size: float = 0.5,
    ncol: int = 6,
    filename: Optional[str] = "space_2d.png",
    show: bool = False,
):
    """
    Args:
        adata: AnnData object.
        groupby: The key that stores clustering or annotation information in adata.obs,
                 a gene's name or a list of genes' name in adata.var.
        spatial_key: The key in `.obsm` that corresponds to the spatial coordinate of each bucket.
        slice_key: The key in `.obs` that corresponds to the slice labels.
        colormap: Colors to use for plotting data.
        alphamap: The opacity of the color to use for plotting data.
        mask: The part that you don't want to be displayed.
        mask_color: Color to use for plotting mask information.
        mask_alpha: The opacity of the color to use for plotting mask information.
        point_size: The size of the plotting points.
        ncol: The maximum number of subplots that can be drawn per row in the figure.
        filename: Filename of output file. Writer type is inferred from the extension of the filename.
        show: Whether to open the visualization window to display the image.
    """

    from plotnine import (
        aes,
        element_blank,
        element_rect,
        element_text,
        facet_wrap,
        geom_point,
        ggplot,
        ggsave,
        labs,
        scale_color_cmap,
        scale_color_gradient,
        scale_color_manual,
        theme,
        theme_classic,
    )

    warnings.filterwarnings("ignore")

    # The final data used for plotting.
    plot_data = pd.DataFrame()

    # The `spatial_key` array in original adata.obsm.
    plot_data["x"] = adata.obsm[spatial_key][:, 0].astype(np.float64)
    plot_data["y"] = adata.obsm[spatial_key][:, 1].astype(np.float64)

    # Slice number and quantity.
    plot_data["slice"] = (
        np.array(["slice"] * adata.obs.shape[0])
        if slice_key is None
        else adata.obs[slice_key].values
    )
    slices = plot_data["slice"].drop_duplicates().values.tolist()
    n_slices = int(len(slices))
    ncol = n_slices if n_slices <= ncol else ncol
    nrow = math.ceil(n_slices / ncol)

    # The`groupby` array in original adata.obs or adata.X.
    mask_list = mask if isinstance(mask, list) else [mask]
    obs_names = set(adata.obs_keys())
    gene_names = set(adata.var_names.tolist())

    if groupby is None:
        title = "Coordinates"

        groups_type = "Groups"
        plot_data[groups_type] = np.array(["None"] * adata.obs.shape[0])
        plot_data["Alpha"] = 1.0 if alphamap == "auto" else alphamap
        _cmap = scale_color_manual(values=mask_color)

    elif groupby in obs_names:
        title = "Clustering"

        groups_type = "Groups"
        plot_data[groups_type] = (
            adata.obs[groupby]
            .map(lambda x: "mask" if x in mask_list else str(x))
            .values
        )
        alphamap = 1.0 if alphamap == "auto" else alphamap
        plot_data["Alpha"] = plot_data[groups_type].map(
            lambda x: mask_alpha if x == "mask" else alphamap
        )

        gc_dict = groups_cmap(
            groups=plot_data[groups_type].unique().tolist(),
            colormap=colormap,
            mask_color=mask_color,
        )
        _cmap = scale_color_manual(values=gc_dict)

    elif groupby in gene_names or set(groupby) <= gene_names:
        gn = "+".join(groupby) if isinstance(groupby, tuple) else groupby
        title = f"Gene(s) name: {gn}"
        groups_type = "Gene_exp"
        groupby = list(groupby) if isinstance(groupby, tuple) else [groupby]
        plot_data[groups_type] = adata[:, groupby].X.sum(axis=1).round(2)
        p5 = plot_data[groups_type].quantile(0.05)
        p95 = plot_data[groups_type].quantile(0.95)
        plot_data[groups_type][plot_data[groups_type]>p95] = p95
        plot_data[groups_type][plot_data[groups_type]<p5] = p5
        plot_data["Alpha"] = (
            minmax_scale(plot_data["Gene_exp"]) * 0.5 + 0.5
            if alphamap == "auto"
            else alphamap
        )
        if colormap == "auto":
            _cmap = scale_color_gradient(low="gainsboro", high="darkblue")
        else:
            _cmap = scale_color_cmap(cmap_name=colormap)

    else:
        raise ValueError(
            "\n`groupby` value is wrong."
            "\n`groupby` can be a string and one of adata.obs_names or adata.var_names. "
            "\n`groupby` can also be a list and is a subset of adata.var_names"
        )

    # Set ggplot object.
    per_x = plot_data["x"].max() - plot_data["x"].min()
    per_y = plot_data["y"].max() - plot_data["y"].min()

    q = (
        ggplot(aes(x="x", y="y"), plot_data)
        + geom_point(aes(color=groups_type, alpha="Alpha"), size=point_size)
        + _cmap
        + facet_wrap("slice", ncol=ncol)
        + labs(title=title)
        + theme_classic()
        + theme(
            text=element_text(color="black", weight="bold", family="sans-serif"),
            axis_text_x=element_blank(),
            axis_text_y=element_blank(),
            figure_size=(math.ceil(per_x / per_y * 3) * ncol, 3 * nrow),
            panel_border=element_rect(),
            dpi=300,
        )
    )

    if not (filename is None):
        ggsave(q, filename, limitsize=False)

    if show is True:
        print(q)
