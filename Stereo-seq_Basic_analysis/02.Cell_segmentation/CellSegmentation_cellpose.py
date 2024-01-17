#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
from pickle import TRUE
import cv2
from typing import Literal, Optional, Union,Tuple
import numpy as np
import spateo as st
from anndata import AnnData
from cellpose import models
from skimage.filters import threshold_local, threshold_multiotsu, unsharp_mask
from skimage.morphology import dilation, disk

import matplotlib.pyplot as plt
import os,re,sys

from cellpose import utils,transforms
from cellpose import models, io
from cellpose import plot
import pandas as pd
from shapely import wkb


# In[2]:


def equalhist_transfer(
    img, method="global", cliplimit=20, tileGridSize=(7, 7)
) -> np.ndarray:
    """
    Histogram equalization for image enhancement, including:
        global histogram equalization,
        adaptive local histogram equalization.

    Args:
        img: Single-channel uint8 type grayscale image[0, 255].
        method: Histogram equalization methods. Available `method` are:
            * `'global'`: global histogram equalization;
            * `'local'`: adaptive local histogram equalization.
        cliplimit: Threshold to limit contrast when method = 'local'.
        tileGridSize: Size of grid for histogram equalization when method = 'local'.

    Returns:
        Image matrix after image enhancement.
    """
    if method == "global":
        # Global histogram equalization.
        return cv2.equalizeHist(img)
    elif method == "local":
        # Adaptive local histogram equalization.
        clahe = cv2.createCLAHE(clipLimit=cliplimit, tileGridSize=tileGridSize)
        return clahe.apply(img)
    else:
        raise ValueError(
            f"Available histogram equalization methods are `global` and `local`."
        )
def generate_bin1_image(
    gem_path: str,
    EH_method: Literal["global", "local"] = "global",
    EH_cliplimit: int = 20,
    color_flip: bool = False,
    gray_factor: int = 1,
    output_path: Optional[str] = None,
) -> Tuple[pd.DataFrame, np.ndarray]:
    """
    Preprocessing of original lasso data.

    Args:
        gem_path: Path to read file.
        EH_method: Histogram equalization methods. Available `method` are:
            * `'global'`: global histogram equalization;
            * `'local'`: adaptive local histogram equalization.
        EH_cliplimit: Threshold to limit contrast when method = 'local' .
        color_flip: If color_flip is True, flip the color of the image.
        gray_factor:  Increasing the value in the grayscale image.
        output_path: If save is not None, save is the path to save the image.

    Return:
        The bin1 image.
    """

    # Read bgi file and generate the original image.
    print(f"Generate bin1 image ... {gem_path}")
    adata = st.io.read_bgi_agg(path=gem_path)
    image = adata.X.todense() if isspmatrix(adata.X) else np.asarray(adata.X)

    #  Increasing the value in the grayscale image.
    image = image * gray_factor
    image[image > 255] = 255

    # Image enhancement based on Histogram equalization.
    image = equalhist_transfer(image.astype("uint8"), method=EH_method, cliplimit=EH_cliplimit)

    # Flip the color of the image.
    image = cv2.bitwise_not(src=image) if color_flip else image

    # Save the image.
    cv2.imwrite(output_path, image)
    print(f"Successfully ... {output_path}")
    return image


# In[3]:


def cellpose(
    img: np.ndarray,
    model_type: Literal["cyto", "nuclei"] = "nuclei",
    **kwargs,
) -> np.ndarray:
    """Run Cellpose on the provided image.
    Args:
        img: Image as a Numpy array.
        model_type: Cellpose model to use. Can be one of the two pretrained
        models:
            * cyto: Labeled cytoplasm
            * nuclei: Labeled nuclei
            Or any generic CellposeModel model.
    **kwargs: Additional keyword arguments to :func:`Cellpose.eval`
    function.
        Returns:
        Numpy array containing cell labels.
    """
    model = models.Cellpose(model_type=model_type)
    res, _, _, _ = model.eval(img, channels=[0, 0], **kwargs)
    return res

def custom_segmentation(
    path_bin1: str,
    path_ssdna: str,
    radius: Optional[int] = None,
    amount: int = 10,
    block_size: int = 11,
    offset: int = 0,
    otsu_class: int = 3,
    otsu_index: int = 0,
    dilation_distance: float = 2.0,
    diameter: Optional[int] = None,
    flow_threshold: float = 0.8,
    min_size: int = -1,
    **kwargs,
) -> AnnData:
    """Run custom_segmentation
    Args:
        path_bin1: Bin1 lasso file start with coordinate 0
        path_ssdna: Cropped ssDNA
        radius: Default:None
        amount: Recommend 5-30,The higher the value, the brighter the
        image(overexposure).Default:10
        block_size: Image block processing with uneven brightness. Default:9
        offset: Try increase this value if too few are reserved. Default:0
        otsu_class: Classify image histograms into several categories. Bigger
        is slower.Default:3
        otsu_index: Pix less than this threshold will be filtered out. If the
        final image retains too much noise cell, increase it. range:(0,otsu_class - 2).
        Default:1
        dilation_distance: If the cell segmentation retains too many noise
        cells, decrease it. Default:2
        diameter: (float (optional, default 1)), if set to None, then
        diameter is automatically estimated if size model is loaded
        flow_threshold: (float (optional, default 0.8)), flow error
        threshold (all cells with errors below threshold are kept) (not used for 3D)
        min_size:(int (optional, default -1)), minimum number of pixels per
        mask, can turn off with -1
        **kwargs: Additional keyword arguments to :func:`Cellpose.eval`
        function.https://cellpose.readthedocs.io/en/latest/api.html#cellposemodel
    """
    adata = st.io.read_bgi_agg(path_bin1, path_ssdna)
    img = cv2.imread(path_ssdna, 0)
    
    if radius != None:
        image = unsharp_mask(img, radius=radius, amount=amount)
    else:
        image = img
    
    local_thresh = threshold_local(image, block_size, offset=offset)
    image[image < local_thresh] = 0
    
    thresholds = threshold_multiotsu(image, classes=otsu_class)
    
    image[image < thresholds[otsu_index]] = 0
    
    if dilation_distance != 0:
        footprint = disk(dilation_distance)
        image = dilation(image, footprint)
    
    segmented_cells = cellpose(
        img=image,
        diameter=diameter,
        flow_threshold=flow_threshold,
        min_size=min_size,
        #masks = masks, 
        **kwargs,
    ) #
    # add
    adata.layers["stain_mask"] = image
    adata.layers["cellpose_labels"] = segmented_cells
    return adata

def image_to_rgb(img0, channels=[0,0]):
    """ image is 2 x Ly x Lx or Ly x Lx x 2 - change to RGB Ly x Lx x 3 """
    img = img0.copy()
    img = img.astype(np.float32)
    if img.ndim<3:
        img = img[:,:,np.newaxis]
    if img.shape[0]<5:
        img = np.transpose(img, (1,2,0))
    if channels[0]==0:
        img = img.mean(axis=-1)[:,:,np.newaxis]
    for i in range(img.shape[-1]):
        if np.ptp(img[:,:,i])>0:
            img[:,:,i] = np.clip(transforms.normalize99(img[:,:,i]), 0, 1)
            img[:,:,i] = np.clip(img[:,:,i], 0, 1)
    img *= 255
    img = np.uint8(img)
    RGB = np.zeros((img.shape[0], img.shape[1], 3), np.uint8)
    if img.shape[-1]==1:
        RGB = np.tile(img,(1,1,3))
    else:
        RGB[:,:,channels[0]-1] = img[:,:,0]
        if channels[1] > 0:
            RGB[:,:,channels[1]-1] = img[:,:,1]
    return RGB

def get_seg_labels(seg_adata,cell_layer):
    # cell labels
    import pandas as pd
    cells = pd.DataFrame(seg_adata.layers[cell_layer])
    cells["x"] = cells.index
    cells = pd.melt(cells, id_vars=["x"])
    cells = cells[cells["value"] != 0]
    cells.columns = ["x", "y", cell_layer]
    cells.sort_values(by=[cell_layer, "x", "y"], inplace=True)
    return cells

def plot_seg(long_df,img,cell_layer):
        long_df = long_df.reset_index(drop=True)
        check_labels = np.zeros([long_df['x'].max()+1, long_df['y'].max()+1], dtype=np.int32)
        check_labels[long_df['x'].astype(int), long_df['y'].astype(int)] = long_df[cell_layer]
        img0 = img.copy()
        if img0.shape[0] < 4:
            img0 = np.transpose(img0, (1,2,0))
        if img0.shape[-1] < 3 or img0.ndim < 3:
            img0 = image_to_rgb(img0, channels=[[0,0]])
        else:
            if img0.max()<=50.0:
                img0 = np.uint8(np.clip(img0*255, 0, 1))

        outlines_test = utils.masks_to_outlines(check_labels)
        outX, outY = np.nonzero(outlines_test)

        imgout= img0.copy()
        imgout[outX, outY] = np.array([255,0,0])
        return imgout




out = '../Cellseg/'


project = 'L1-early_b'

outdir = os.path.join(out,project)
if not os.path.exists(outdir):
    os.mkdir(outdir)


slices = 'S01'


# In[4]:


path_bin1 =project+"/"+project+"_"+slices+"_bin1_final.gem.gz"
path_ssdna =project+"_"+slices+"_cropped_img.tif"

print('path_bin1:',path_bin1)
print('path_ssdna:',path_ssdna)



dilation_distance = 1   #2
amount = 40   #10
block_size = 199  # 11
radius = 2  # None
otsu_class = 5  # 3
otsu_index = 0   # 0 
diameter = None  # None
flow_threshold = 0.8   # 0.8
min_size = -1   # -1

seg_adata = custom_segmentation(
    
            path_bin1 = path_bin1,
            path_ssdna = path_ssdna,
    
            dilation_distance = dilation_distance,
            amount = amount,
            block_size = block_size,
            radius = radius,  
            otsu_class = otsu_class,
            otsu_index = otsu_index,       
            diameter = diameter,
            flow_threshold = flow_threshold,
            min_size = min_size,
       )

# cellseg plot
long_df = get_seg_labels(seg_adata,'cellpose_labels')
img = io.imread(path_ssdna)
imgout = plot_seg(long_df,img,'cellpose_labels')
plt.figure(figsize=(8,8))
plt.imshow(imgout)
plt.show()


plt.close()


# expand
exp_distance = 2
st.cs.expand_labels(seg_adata, "cellpose_labels", distance = exp_distance, max_area=np.inf, out_layer = "cell_labels_expanded")

# expand plot
long_df = get_seg_labels(seg_adata,'cell_labels_expanded')
img = io.imread(path_ssdna)

imgout = plot_seg(long_df,img,'cell_labels_expanded')

#plt.figure(figsize=(6,2))
plt.figure(figsize=(8,8))
plt.imshow(imgout)
plt.show()
#plt.title(slice,fontsize = 30)
#plt.axis('off')
#plt.savefig(img_save_path)


k = "_"+str(dilation_distance)+"_"+str(amount)+"_"+str(block_size)+"_"+str(radius)+"_"+str(otsu_class)+"_"+str(otsu_index)+'_'+str(diameter)+'_'+str(flow_threshold)+'_'+str(min_size)+'_'+str(exp_distance)
k

seg_adata_save = os.path.join(outdir,project+'_'+slices+k+'_seg.h5ad')
seg_fig_save = os.path.join(outdir,project+'_'+slices+k+'_seg.png')
seg_all_fig = os.path.join(outdir,project+'_'+slices+k+'_seg_all.pdf')
adata_save = os.path.join(outdir,project+'_'+slices+k+'_cellbin.h5ad')

# seg_adata save
seg_adata.write(seg_adata_save, compression="gzip")

# seg qc save
#plt.figure(figsize=(20,10))
long_df = get_seg_labels(seg_adata,'cell_labels_expanded')
img = io.imread(path_ssdna)
imgout = plot_seg(long_df,img,'cell_labels_expanded')

fig, axes = plt.subplots(ncols=2, figsize=(30,15), tight_layout=True)
axes[0].imshow(img)
axes[1].imshow(imgout)
plt.savefig(seg_fig_save)

# all seg plot
n_layers = len(seg_adata.layers.keys())
fig, axes = plt.subplots(ncols=n_layers+2, figsize=(30,6), tight_layout=True)
axes[0].imshow(img)
for i, layer in enumerate(seg_adata.layers.keys()):
    st.pl.imshow(seg_adata, layer, ax=axes[i+1], labels=True, save_show_or_return="return")
axes[-1].imshow(imgout)
plt.savefig(seg_all_fig)

# adata save
adata = st.io.read_bgi(
            path=path_bin1,
            segmentation_adata=seg_adata,
            labels_layer='cell_labels_expanded',
            seg_binsize=1,
        )
adata.obs["slices"] = str(project)+'_'+slices
adata.write(adata_save, compression="gzip")

plt.close()


