#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os,re,sys
from pathlib import Path
import spateo as st
import pandas as pd
import numpy as np
import anndata as ad
from typing import Any, List, Optional, Tuple, Union
from scipy.sparse import csr_matrix, isspmatrix, spmatrix
from anndata import AnnData


# In[2]:


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


# In[2]:


project = 'L1-early_a'


# In[3]:


cellbin_dir = os.path.join('../Cellseg/',project,'cellbin_data')


# In[22]:


k = os.listdir(cellbin_dir)
k.sort()
k


# In[23]:


slices = []
for i in k:
    adata_path = os.path.join(cellbin_dir,i)
    adata = ad.read(adata_path)
    adata.obs['slices'] = [i[11:] for i in adata.obs['slices'].tolist()]
    slices.append(adata)
slices, len(slices)


# In[25]:


st.pl.multi_models(
    slices,
    spatial_key="spatial",
    id_key="slices",
    mode="single",
    cpo="xy",
    jupyter="trame",
    shape=(5,5),
    filename=os.path.join(out_image_path, f"{sample_id}_raw_slices_align_space.png")
)


# In[26]:


align_slices, pis, _ = st.align.morpho_align(
    models=slices,
    layer="X",
    spatial_key="spatial",
    key_added="align_spatial",
    # beta2=10,
    max_iter=400,
    #device="0",
    SVI_mode=False,
    verbose=False
)

'''
#use GPU
align_slices, pis, _ = st.align.morpho_align(
    models=slices,
    layer="X",
    spatial_key="spatial",
    key_added="align_spatial",
    # beta2=10,
    max_iter=400,
    device="0",
    SVI_mode=True,
    verbose=False
)

'''


# In[29]:


st.pl.multi_models(
    align_slices,
    spatial_key="align_spatial",
    id_key="slices",
    mode="single",
    cpo="xy",
    jupyter="trame",
    shape=(5,5),
    filename=os.path.join(out_image_path, f"{sample_id}_raw_slices_align_space.png")
)


# In[67]:


st.pl.multi_models(
    align_slices,
    spatial_key="Rigid_3d_align_spatial",
    id_key="slices",
    mode="single",
    cpo="xy",
    jupyter="trame",
    shape=(5,5),
    filename=os.path.join(out_image_path, f"{sample_id}_raw_slices_align_space.png")
)


# In[68]:


align_slices


# In[ ]:


outdir = '../Cellseg/L1-early_a/02.align/align_cellbin/'
# 
for i in align_slices:
    if 'iter_spatial' in i.uns:
        del i.uns['iter_spatial']
    if 'VecFld_morpho' in i.uns:
        del i.uns['VecFld_morpho']
    adata_path = os.path.join(outdir,i.obs['slices'][0]+'_cellbin.h5ad')
    i.write_h5ad(adata_path,compression='gzip')


all_adata = integrate(adatas=align_slices, batch_key="slices")

all_adata.write_h5ad(os.path.join('../Cellseg/L1-early_a/02.align/','L1-early_a_cellbin_align.h5ad'),compression='gzip')





