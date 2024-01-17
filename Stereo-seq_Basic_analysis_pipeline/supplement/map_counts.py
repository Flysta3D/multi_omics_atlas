import os
from typing import Optional

import numpy as np
import pandas as pd
import spateo as st
from anndata import AnnData
from scipy.sparse import isspmatrix


def map_counts(
    norm_adata: AnnData,
    counts_adata: AnnData,
    key_added: Optional[str] = None,
    norm_X: Optional[str] = None,
    counts_X: Optional[str] = None,
) -> AnnData:
    """
    Mapping the gene count matrix to the normalized gene expression matrix.
    Args:
        norm_adata: An adata containing normalized gene expression matrix.
        counts_adata: An adata containing gene count matrix.
        key_added: The key under which to add the gene count matrix.
                   If key_added is None, it defaults to norm_adata.X.
        norm_X: The key under which to contain the normalized gene expression matrix.
                If norm_X is None, it defaults to norm_adata.X.
        counts_X: The key under which to contain gene count matrix.
                  If counts_X is None, it defaults to counts_adata.X.
    Returns:
        norm_adata: An anndata object.
    """
    norm_adata = norm_adata.copy()
    counts_adata = counts_adata.copy()

    # Check obs_names
    norm_obs_index = norm_adata.obs_names.values.flatten()
    counts_obs_index = counts_adata.obs_names.values.flatten()
    obs_diff = np.setdiff1d(
        ar1=norm_obs_index, ar2=counts_obs_index, assume_unique=False
    )
    if len(obs_diff) != 0:
        raise ValueError(
            "There is an error in the `obs_names`. "
            "Some values in `norm_adata.obs_names` are not present in `counts_adata.obs_names`."
        )

    # Chekc var_names
    norm_var_index = norm_adata.var_names.values.flatten()
    counts_var_index = counts_adata.var_names.values.flatten()
    var_diff = np.setdiff1d(
        ar1=norm_var_index, ar2=counts_var_index, assume_unique=False
    )
    if len(var_diff) != 0:
        raise ValueError(
            "There is an error in the `var_names`. "
            "Some values in `norm_adata.var_names` are not present in `counts_adata.var_names`."
        )

    # Check X
    def check_X(adata, X_label, X_label_name):
        if X_label is None:
            X = adata.X
        elif X_label in adata.obsm_keys():
            X = adata.obsm[norm_X]
        elif X_label in adata.layers.keys():
            X = adata.layers[norm_X]
        else:
            raise ValueError(f"`{X_label_name}` value is wrong.")

        to_dense_matrix = (
            lambda x: np.array(x.todense()) if isspmatrix(x) else np.asarray(x)
        )
        return to_dense_matrix(X)

    n_X = check_X(adata=norm_adata, X_label=norm_X, X_label_name="norm_X")
    c_X = check_X(adata=counts_adata, X_label=counts_X, X_label_name="counts_X")

    # map counts matrix to norm matrix
    norm_data = pd.DataFrame(
        n_X, columns=norm_adata.var_names.values, index=norm_adata.obs_names.values
    )
    counts_data = pd.DataFrame(
        c_X, columns=counts_adata.var_names.values, index=counts_adata.obs_names.values
    )
    map_data = counts_data.loc[norm_data.index, norm_data.columns]

    if key_added is None:
        norm_adata.X = map_data.values
    else:
        norm_adata.layers[key_added] = map_data.values
    return norm_adata