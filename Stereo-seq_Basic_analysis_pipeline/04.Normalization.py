#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os,re,sys
import pandas as pd
import numpy as np
import anndata as ad
import spateo as st
import dynamo as dyn


# In[2]:


from typing import Optional, Union
import os
import matplotlib.pyplot as plt
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


# In[3]:


from typing import Any, List, Optional, Tuple, Union
from scipy.sparse import csr_matrix, isspmatrix, spmatrix
from anndata import AnnData
def integrate(
    adatas: List[AnnData],
    batch_key: str = "slices",
) -> AnnData:
    """
    Concatenating all anndata objects.
    Args:
        adatas: AnnData matrices to concatenate with.
        batch_key: Add the batch annotation to :attr:`obs` using this key.
    Returns:
        integrated_adata: The concatenated AnnData, where adata.obs[batch_key] stores a categorical variable labeling the batch.
    """
    to_dense_matrix = lambda X: np.array(X.todense()) if isspmatrix(X) else np.asarray(X)
    
    # Merge the obsm data and varm data of all anndata objcets separately.
    obsm_dict, varm_dict = {}, {}
    obsm_keys, varm_keys = adatas[0].obsm.keys(), adatas[0].varm.keys()
    n_obsm_keys, n_varm_keys = len(obsm_keys), len(varm_keys)
    if n_obsm_keys > 0:
        for key in obsm_keys:
            if key != 'X_pca':
                if key != 'X_umap':
                    obsm_dict[key] = np.concatenate([to_dense_matrix(adata.obsm[key]) for adata in adatas], axis=0)
    if n_varm_keys > 0:
        for key in varm_keys:
            varm_dict[key] = np.concatenate([to_dense_matrix(adata.varm[key]) for adata in adatas], axis=0)
    
    # Delete obsm and varm data.
    for adata in adatas:
        del adata.obsm, adata.varm
    
    # Concatenating obs and var data.
    batch_ca = [adata.obs[batch_key][0] for adata in adatas]
    integrated_adata = adatas[0].concatenate(
        *adatas[1:], 
        #batch_key=batch_key, 
        #batch_categories=batch_ca, 
        join="outer", fill_value=0, uns_merge="unique"
    )
    
    # Add Concatenated obsm data and varm data to integrated anndata object.
    if n_obsm_keys > 0:
        for key, value in obsm_dict.items():
            integrated_adata.obsm[key] = value
    if n_varm_keys > 0:
        for key, value in varm_dict.items():
            integrated_adata.varm[key] = value
    
    return integrated_adata


# In[8]:


sample = 'PUPA-24h'

align_cellbin_dir = f'../Cellseg/{sample}/02_align/'


# In[9]:


k = os.listdir(align_cellbin_dir)
k.sort()
k,len(k)


# In[10]:


# preprocessing
raw_slices = []
norm_slices = []

for i in k:
    adata_path = os.path.join(align_cellbin_dir,i)
    adata = ad.read(adata_path)
    z = adata.obs["slices"].map(lambda x: int(x[-2:]) * 14)
    adata.obsm["align_spatial"] = np.c_[adata.obsm["align_spatial"], z]
    basic_stats(adata=adata, mito_label="mt:")
    raw_slices.append(adata)
    
    adata = adata[adata.obs.pMito < 0.1, :]
    st.pp.filter.filter_cells(adata=adata, min_area=20, min_expr_genes=20, inplace=True)
    st.pp.filter.filter_genes(adata=adata, min_cells=3, min_counts=1, inplace=True)
    basic_stats(adata=adata, mito_label="mt:")
    st.tl.pearson_residuals(adata=adata, n_top_genes=None)
    #adata.write_h5ad(filename=os.path.join(save_h5ad_folder, bin_file), compression="gzip")
    adata.layers["counts_X"] = adata.X.copy()
    
    norm_adata = adata.copy()
    norm_adata.X = norm_adata.obsm["pearson_residuals"]
    del norm_adata.layers, norm_adata.obsm["pearson_residuals"]
    norm_slices.append(norm_adata)


# In[11]:


# integration
raw_adata = integrate(adatas=raw_slices, batch_key="slices")
norm_adata = integrate(adatas=norm_slices, batch_key="slices")


# In[17]:


norm_adata.obs['slices'] = ['Pupa-24h'+i[8:] for i in norm_adata.obs['slices'].tolist()]


# In[18]:


raw_adata.obs['slices'] = ['Pupa-24h'+i[8:] for i in raw_adata.obs['slices'].tolist()]


# In[19]:


raw_adata


# In[20]:


norm_adata


# In[21]:


norm_adata.layers['counts_X'] = raw_adata[norm_adata.obs.index,norm_adata.var_names].X

del norm_adata.var


# In[22]:


norm_adata


# In[23]:


#save
raw_adata.write_h5ad(f'../Cellseg/{sample}/03.norm/{sample}_raw.h5ad',compression='gzip')


# In[24]:


norm_adata.write_h5ad(f'../Cellseg/{sample}/03.norm/{sample}_norm.h5ad',compression='gzip')

