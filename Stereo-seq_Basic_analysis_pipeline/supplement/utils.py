import os
from typing import Literal

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.sparse import csr_matrix, issparse


def check_nan_gene(path: str, nan_gene: str = "CG5842") -> pd.DataFrame:
    """
    Since some gene names may be `nan`, python cannot recognize them as normal strings, so we need to replace the gene names with other names.

    Args:
        path: Path to read file.
        nan_gene: The replaced gene name.
    """
    data = pd.read_csv(
        path,
        sep="\t",
        dtype={
            "geneID": str,  # geneID
            "x": np.uint32,  # x
            "y": np.uint32,  # y
            "MIDCounts": np.uint16,  # total
        },
        comment="#",
    )

    data["geneID"].fillna(value=nan_gene, inplace=True)
    return data


def convert_matrix(matrix, direct: Literal["to_sparse", "to_dense"] = "to_sparse"):
    """Convert between sparse and dense matrices."""

    # Convert sparse matrix to dense matrix.
    to_dense_matrix = lambda X: np.array(X.todense()) if issparse(X) else np.asarray(X)
    # Convert dense matrix to sparse matrix.
    to_sparse_matrix = lambda X: csr_matrix(X) if not issparse(X) else np.asarray(X)

    if direct == "to_sparse":
        return to_sparse_matrix(matrix)
    elif direct == "to_dense":
        return to_dense_matrix(matrix)


def convert_X_matrix(
    adata: AnnData,
    direct: Literal["to_sparse", "to_dense"] = "to_sparse",
    inplace: bool = True,
):
    """Converts the format of the Anndata.X matrix between sparse and dense matrices."""

    adata = adata if inplace else adata.copy()
    adata.X = convert_matrix(matrix=adata.X, direct=direct)
    return None if inplace else adata


def convert_layer_matrix(
    adata: AnnData,
    layer: str,
    direct: Literal["to_sparse", "to_dense"] = "to_sparse",
    inplace: bool = True,
):
    """Converts the format of the Anndata.layers[layer] matrix between sparse and dense matrices."""

    adata = adata if inplace else adata.copy()
    adata.layers[layer] = convert_matrix(matrix=adata.layers[layer], direct=direct)
    return None if inplace else adata
