#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os,sys,getopt,re
try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal

import warnings
import pandas as pd
import spateo as st
warnings.filterwarnings("ignore")
import math
import operator
from typing import Optional, Union
import matplotlib as mpl
import numpy as np
from anndata import AnnData
from scipy.sparse import isspmatrix
from sklearn.preprocessing import minmax_scale
import anndata as ad

import matplotlib.pyplot as plt
import seaborn as sns

import dynamo as dyn
import scanpy as sc
from scipy import sparse
import logging as logg
import time

try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal

from plotnine import *


# In[4]:


def cluster_small_multiples(adata, clust_key, size=60, frameon=False, legend_loc=None, **kwargs):
    tmp = adata.copy()
    
    for i,clust in enumerate(adata.obs[clust_key].cat.categories):
        tmp.obs[clust] = adata.obs[clust_key].isin([clust]).astype('category')
        if type(adata.uns[clust_key+'_colors'])==dict:
            tmp.uns[clust+'_colors'] = ['#d3d3d3', list(adata.uns[clust_key+'_colors'].values())[i]]
        else:
            tmp.uns[clust+'_colors'] = ['#d3d3d3', adata.uns[clust_key+'_colors'][i]]

    sc.pl.umap(tmp, 
               groups=tmp.obs[clust].cat.categories[1:].values, 
               color=adata.obs[clust_key].cat.categories.tolist(), 
               size=size, 
               frameon=frameon, 
               legend_loc=legend_loc, 
               show = False,
               **kwargs)
    
    del tmp


def rank_genes_groups_df(adata, key='rank_genes_groups',subgroup = None):
        dd = []
        groupby = adata.uns[key]['params']['groupby']
        if subgroup != None:
            cate = subgroup
        else:
            cate = adata.obs[groupby].cat.categories
        #print(cate)
        for group in cate:
            cols = []
            # inner loop to make data frame by concatenating the columns per group
            for col in adata.uns[key].keys():
                if col != 'params':
                       cols.append(pd.DataFrame(adata.uns[key][col][str(group)], columns=[col]))
            df = pd.concat(cols,axis=1)
            df['group'] = group
            dd.append(df)
        # concatenate the individual group data frames into one long data frame
        rgg = pd.concat(dd)
        rgg['group'] = rgg['group'].astype('category')
        return rgg

def get_top(df,num):
    top = []
    for i in df['group'].unique().tolist():
        ddf = df[df['group']==i].head(int(num))
        top.append(ddf)
    gene_df = pd.concat(top)
    return gene_df

def cluster(norm_adata,save_folder,resolution,sample_name = "S"):
    key = "scc_"+str(resolution)
    save_image_folder = os.path.join(save_folder,key)
    if not os.path.exists(save_image_folder):
        os.makedirs(save_image_folder)
    # clu
    st.tl.scc(
        adata=norm_adata,
        spatial_key="align_spatial", key_added=key, pca_key="X_pca",
        e_neigh=30, s_neigh=6,cluster_method="louvain",
        resolution=resolution,
        copy=False
    )
    # fig
    os.chdir(save_image_folder)
    sc.pl.umap(norm_adata,color = ['slices',key],save = '_'+key+'.png')
    
    sc.tl.rank_genes_groups(norm_adata, key, method="t-test")
    df = rank_genes_groups_df(norm_adata)
    gene_df = get_top(df,50)
    gene_df.to_csv(key+'_markers_top50.csv')
    
    top10 = get_top(df,10)['names'].unique().tolist()
    sc.pl.dotplot(norm_adata,top10, groupby= key,save='_'+key+'_top10_markers.pdf')
    
    cluster_small_multiples(norm_adata,key,size = 2, save = '_'+key+'_eachclu_on_umap.png')
    
    pal = ["#DC143C","#0000FF","#20B2AA","#FFA500","#9DBEB9","#7CFC00","#FFFF00","#666666","#DBF6E9","#C9CBFF","#00416D","#40A8C4","#E5D549",
          "#808000","#A03C78","#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#D35D6E","#EFB08C","#F5F1DA","#FF5722",
          "#ADE498","#B1B493","#35D0BA","#864000","#E4BAD4","#00e8ff","#ff00dd","#b5ff00","#ff7700","#006dff","#00ff91","#1aab27","#7f90f0","#b774c5","#a28706","#5a2686","#003fff","#6787e8","#488d6b","#ff00a0"]
    
    pic = {
        'E0-1h_e':[0.3,(25,8)],
    'E0-1h_c':[0.1,(12,11)],
       'E1-2h_d':[0.5,(25,6)],
    'E2.5-3h_c':[0.1,(25,6)],
    'E2.5-3h_d':[0.1,(24,6)],
    'E2.5-3h_f':[0.1,(24,4)],
    'E3-3.5h_c':[0.1,(24,12)],
    'E2-2.5h_b':[0.1,(12,12)],
    'E2-2.5h_a':[0.2,(20,9)],
       'E3.5-4h_d':[0.8,(14,20)],
       'E3.5-4h_g':[0.6,(20,20)],
       'E4-8_e':[0.3,(25,15)],
       'E4-8_h':[0.5,(20,10)],
    'E4-8h_f':[0.5,(16,8)],
       'E4-8_g':[0.2,(28,7)],
       'E8-10_c':[0.5,(20,10)],
       'E8-10_b':[0.2,(25,10)],
     'E10-12h_d':[0.1,(15,18)],
    'E10-12h_e':[0.1,(16,18)],
       'E10-12h_f':[0.1,(15,18)],
    'E12-14h_a':[0.1,(22,12)],
    'E12-14h_b':[0.1,(24,16)],
       'E12-14h_c':[0.3,(25,20)],
       'E14-16h_d':[0.2,(16,28)],
       'E14-16h_c':[0.3,(16,18)],
       'E14-16h_f':[0.2,(14,19)],
    'E16-18h_a':[0.1,(16,12)],
       'E16-18h_e':[0.3,(28,14)],
       'E16-18h_b':[0.1,(18,14)],
    'E18-20h_d':[0.1,(14,18)],
    'E20-22h_a':[0.1,(30,8)],
    'E20-22h_c':[0.2,(24,4)],
    'E20-22h_e':[0.2,(20,10)],

    'L1-early_a':[0.2,(25,28)],
    'L1-early_b':[0.2,(25,28)],
    'L1-early_e':[0.2,(25,28)],
       'L1-late_c':[0.1,(25,32)],
       'L1-late_d':[0.1,(38,15)],
    'L2-early_a':[0.1,(16,35)],
       'L2-early_b':[0.1,(20,35)],
       'L2-late_c':[0.1,(20,40)],
       'L2-late_d':[0.1,(36,16)],

    'L3-early':[0.2,(18,32)],
    'L3-late':[0.1,(16,56)],
    'Pupa-12h_a':[0.1,(25,25)],
    'Pupa-24h_a':[0.1,(35,40)],
    'Pupa-36h_a':[0.1,(20,30)],
    'Pupa-48h_a':[0.1,(40,40)],
    'Pupa-60h_a':[0.1,(30,50)],
    'Pupa-72h_a':[0.1,(40,40)],
      }
    
    norm_adata.obs['align_x'] = np.split(norm_adata.obsm['align_spatial'],3,axis=1)[0]
    norm_adata.obs['align_y'] = np.split(norm_adata.obsm['align_spatial'],3,axis=1)[1]
    norm_adata.obs['align_z'] = np.split(norm_adata.obsm['align_spatial'],3,axis=1)[2]
    
    project = sample_name
    dot_size = pic[project][0]
    pic_size = pic[project][1]
    p = (ggplot(aes(x = 'align_x',y = 'align_y'),norm_adata.obs)
             +geom_point(aes(color = key),size = dot_size)
             +facet_wrap('slices',ncol = 6)
             +labs(title = project)
             +scale_color_manual(values=pal)
             #+theme_grey()
             +theme(
                 figure_size = pic_size,
                 plot_background=element_rect(fill='black'),
                 panel_background=element_rect(fill='black'),
                 axis_ticks=element_line(),
                 line=element_line(color='black'),
                 title = element_text(color = 'white'),
                 panel_border=element_rect(color='grey', size=1),
                 legend_title = element_text(color = 'white'),
                 legend_key=element_rect(fill='black', alpha=1),
                 legend_background=element_rect(color='white', size=1, fill='black'),
                 legend_text=element_text(weight='bold',color = 'white'),
                 legend_key_size=30,
                                )
            )
    
    ggsave(p,project+'_'+key+'_sp.pdf',limitsize = False)
    
    each_clu_SP = os.path.join(save_image_folder,'SP_each_clu')
    if not os.path.exists(each_clu_SP):
        os.mkdir(each_clu_SP)
    for c in norm_adata.obs[key].unique():
        p = (ggplot(aes(x = 'align_x',y = 'align_y'),norm_adata.obs)
                 +geom_point(color = '#9C9C9C',size = dot_size)
                 +geom_point(aes('align_x','align_y'),norm_adata.obs[norm_adata.obs[key]==c],color = 'red',size = (dot_size+0.1))
                 #+geom_point(aes(color = 'Annotation_2_tissue'),size = dot_size)
                 +facet_wrap('slices',ncol = 6)
                 +labs(title = project+' cluster : '+c)
                 #+scale_color_manual(values=anno_color)
                 #+theme_grey()
                 +theme(
                     figure_size = pic_size,
                     title = element_text(color = 'white'),
                     plot_background=element_rect(fill='black'),
                     panel_background=element_rect(fill='black'),
                     axis_ticks=element_line(),
                     line=element_line(color='black'),
                     panel_border=element_rect(color='grey', size=1),
                     legend_title = element_text(color = 'white'),
                     legend_key=element_rect(fill='black', alpha=1),
                     legend_background=element_rect(color='white', size=1, fill='black'),
                     legend_text=element_text(weight='bold',color = 'white'),
                     legend_key_size=30,
                                    )
                )
        ggsave(p,os.path.join(each_clu_SP,project+'_'+c+'.png'),limitsize = False)
    
    return(norm_adata)


project = 'L1-early_a'

save_folder = f'G:/tzc/4T_tzc_20230408/00.test/Cellseg/{project}/04.clu/'

adata = ad.read(f'G:/tzc/4T_tzc_20230408/00.test/Cellseg/{project}/03.norm/{project}_norm.h5ad')


sc.pl.violin(adata,keys=['nGenes','nCounts'],groupby='slices',rotation = 90,)


dyn.tl.neighbors(
    adata=adata,
    X_data=adata.obsm['align_spatial'],
    n_neighbors=6,
    n_pca_components=50,
    result_prefix='spatial_'
)
adata.X = adata.X.astype(np.float64)
print(adata.X[np.isnan(adata.X)])

import gc
gc.collect()


st.tl.pca_spateo(adata=adata, n_pca_components=50, pca_key="X_pca")

dyn.tl.reduceDimension(adata=adata, basis="umap", n_components=2)

sc.pl.umap(adata,color = ['slices'])

res = [1.0,1.2,1.4,1.6,1.8,2.0]


## clustering
for resolution in res:
    adata = cluster(adata,save_folder,resolution = resolution,sample_name = project)
    print('===== resolution =',resolution,' finished ! =====')

adata.write_h5ad(os.path.join(save_folder,project+'_clu.h5ad'),compression = 'gzip')


