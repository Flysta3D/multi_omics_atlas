from typing import  Optional, Union,Tuple
try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal
import os,sys,re
from scipy.sparse import issparse
import warnings
import spateo as st
from anndata import AnnData
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")
import pandas as pd
import anndata as ad
import numpy as np
import cv2
import torch

'''
sys.path.append("/hwfssz5/ST_SUPERCELLS/P21Z10200N0206/xiangrong/software/miniconda3/envs/spateo0513/lib/python3.8/site-packages/")

try:
    from cellpose.models import Cellpose, CellposeModel
except ModuleNotFoundError:
    Cellpose = None
    CellposeModel = None



def _cellpose(
    img: np.ndarray,
    model: Union[Literal["cyto", "nuclei"], "CellposeModel"] = "nuclei",
    device: Optional[torch.device] = torch.device('cuda:0'),
    **kwargs,
) -> np.ndarray:
    """
    Run Cellpose on the provided image.
    See Also: https://cellpose.readthedocs.io/en/latest/api.html
    Args:
        img: Image as a Numpy array.
        model: Cellpose model to use. Can be one of the two pretrained models:
            * cyto: Labeled cytoplasm
            * nuclei: Labeled nuclei
            Or any generic CellposeModel model.
        device: Device used for model running / training (torch.device(‘cuda’) or torch.device(‘cpu’)), overrides gpu
                input, recommended if you want to use a specific GPU (e.g. torch.device(‘cuda:1’))
        **kwargs: Additional keyword arguments to :func:`Cellpose.eval`
            function.
    Returns:
        Numpy array containing cell labels.
    """
    if isinstance(model, str):
        model = Cellpose(model_type=model, gpu=True, device=device)  # Use GPU if available
        # model = CellposeModel(pretrained_model='/home/yao/.cellpose/models/size_nucleitorch_0.npy', model_type=model, gpu=True, device=device)

    masks, flows, styles, diams = model.eval(img, **kwargs)
    return masks


def segmentation_spateo(
    path_bin1: str,
    path_ssdna: str,
    seg_method: Literal["watershed", "stardist", "w+s", "cellpose"] = "cellpose",
    seg_kwargs: Optional[dict] = None,
    expand_method: Optional[Literal["certain_distance", "affinity"]] = None,
    expand_kwargs: Optional[dict] = None,
    expand_distance: Union[int, float] = 0,
) -> AnnData:
    """
    Args:
        path_bin1: The file path of Stereo-seq data.
        path_ssdna: The file path of ssDNA (nuclei) staining image.
        seg_method: Methods for segmenting and labeling nuclei from stained images.
                * 'watershed': Spateo includes a custom Watershed-based approach to segment and label nuclei from the
                               staining image. At a high level, it uses a combination of global and local thresholding
                               to first obtain a mask of nuclei (recall that this is called segmentation), and then
                               uses Watershed to assign labels (recall that this is called labeling).
                * 'stardist': Spateo includes a variety of existing deep-learning approaches for fluorescent nuclei
                              segmentation, such as StarDist, Cellpose and DeepCell. We’ve found StarDist to perform
                              most consistently the best among these methods, so that is what we will be using here.
                              Again, note that deep learning-based methods for segmenting cell staining images often
                              combine the segmentation and labeling (i.e. the model outputs the final cell labels).
                * 'w+s': Though both the Watershed and StarDist methods perform well in segmenting nuclei, they each
                          have their limitations. The Watershed approach tends to result in rough cell boundaries due
                          to the nature of thresholding, while StarDist sometimes has difficulty identifying nuclei in
                          dense regions (resulting in “holes”). Additionally, because we apply histogram equalization
                          with the equalize parameter (which can be turned off by setting equalize=-1), sometimes noise
                          in empty regions get amplified and is falsely identified as a cell. We can mitigate these by
                          augmenting the StarDist labels with Watershed labels, copying over Watershed labels that do
                          not overlap with any Stardist labels, and removing Stardist labels that do not overlap with
                          any Watershed labels.
        seg_kwargs:
                When `seg_method` == 'watershed':
                    * 'otsu_index': Which threshold index should be used for classifying cells.
                    * 'min_distance': Minimum distance, in pixels, between peaks.
                                      Set according to the actual cell distance, within the minimum cell distance.
                                      The unit is the number of spots.
                    * 'k': Size of the kernel to use for Gaussian blur.
                           Set to be similar to the diameter of the cell.
                           The unit is the number of spots.
                When `seg_method` == 'stardist':
                    * 'tilesize': Run prediction separately on tiles of size `tilesize` x `tilesize`
                                  and merge them afterwards. Useful to avoid out-of-memory errors. Can be
                                  set to <= 0 to disable tiling. When `min_overlap` is also provided, this
                                  becomes the `block_size` parameter to :func:`StarDist2D.predict_instances_big`.
                                  Using tilesize=-1 to turn off tiling, since this is a small dataset.
                                  Usually, we can leave this to its default value.
                    * 'equalize': Controls the `clip_limit` argument to the :func:`clahe` function.
                                  Set this value to a non-positive value to turn off equalization.
                                  Using equalize=-1 to turn off histogram equalization.
                                  Usually, we can leave this to its default value.
                when `seg_method` == 'cellpose':
                    * 'model': Cellpose model to use. Can be one of the two pretrained models: cyto(Labeled cytoplasm)
                               and nuclei(Labeled nuclei).Or any generic CellposeModel model. Defaults to `"nuclei"`.
                    * 'diameter':  Expected diameter of each segmentation (cells for `model="cyto"`, nuclei for
                                   model="nuclei"`). Can be `None` to run automatic detection.
                    * 'normalize': Whether to percentile-normalize the image.
                    * 'min_size': minimum number of pixels per mask, can turn off with -1.
                    * 'flow_threshold': flow error threshold (all cells with errors below threshold are kept).
        expand_method: Methods for expanding nuclei labels to cell labels.
                * 'certain_distance': Simply expand each label by some distance.
                * 'affinity': we take advantage of the fact that ssDNA staining weakly stains the cytoplasm, and thus
                              we can identify cytoplasmic regions by thresholding the image with a lenient threshold.
        expand_kwargs:
                When `expand_method` == 'certain_distance':
                    * 'distance': Distance to expand.
                                  Internally, this is used as the number of iterations of distance 1 dilation.
                    * 'max_area': Maximum area of each label.
                When `seg_method` == 'affinity':
                    * 'otsu_index': Which threshold index should be used for classifying cells.
        expand_distance: Expand all the cell labels by some amount to mitigate the effects of RNA diffusion.
                         The distance that should be expanded will depend on the level of RNA diffusion in the data.
    Examples:
        When `seg_method` == 'watershed':
        >>> seg_adata = segmentation_spateo(
        ...     path_bin1=path_bin1, path_ssdna=path_ssdna, seg_method="watershed",
        ...     seg_kwargs=dict(otsu_index=0, min_distance=5, k=10)
        ...)
        When `seg_method` == 'stardist':
        >>> seg_adata = segmentation_spateo(
        ...     path_bin1=path_bin1, path_ssdna=path_ssdna, seg_method="stardist",
        ...     seg_kwargs=dict(tilesize=2000, equalize=2.0)
        ...)
        When `seg_method` == 'w+s':
        >>> seg_adata = segmentation_spateo(
        ...     path_bin1=path_bin1, path_ssdna=path_ssdna, seg_method="w+s",
        ...     seg_kwargs=dict(otsu_index=0, min_distance=5, k=10, tilesize=2000, equalize=2.0)
        ...)
        When `seg_method` == 'cellpose':
        >>> seg_adata = segmentation_spateo(
        ...     path_bin1=path_bin1, path_ssdna=path_ssdna, seg_method="cellpose",
        ...     seg_kwargs=dict(model="nuclei", diameter=None, normalize=True, flow_threshold=0.8, min_size=-1, device=torch.device('cuda:0'))
        ...)
    """

    adata = st.io.read_bgi_agg(path_bin1, path_ssdna)

    # Segment and label nuclei from the staining image.
    expand_markers_layer = None
    if seg_method in ["watershed", "w+s"]:
        expand_markers_layer = "watershed_labels"
        default_seg_kwargs = {"otsu_index": 0, "min_distance": 5, "k": 10}
        if not (seg_kwargs is None):
            default_seg_kwargs.update((k, seg_kwargs[k]) for k in default_seg_kwargs.keys() & seg_kwargs.keys())
        st.cs.mask_nuclei_from_stain(adata, otsu_index=default_seg_kwargs["otsu_index"])
        st.cs.find_peaks_from_mask(adata, layer='stain', min_distance=default_seg_kwargs["min_distance"])
        st.cs.watershed(adata, layer="stain", k=default_seg_kwargs["k"], out_layer=expand_markers_layer)
    if seg_method in ["stardist", "w+s"]:
        expand_markers_layer = "stardist_labels"
        st.cs.stardist(adata, layer="stain", out_layer=expand_markers_layer, **seg_kwargs)
    if seg_method == "w+s":
        expand_markers_layer = "augmented_labels"
        st.cs.augment_labels(adata, "watershed_labels", "stardist_labels", out_layer=expand_markers_layer)
    if seg_method == "cellpose":
        expand_markers_layer = "cellpose_labels"
        labels = _cellpose(adata.layers["stain"], channels=[0, 0], **seg_kwargs)
        adata.layers[expand_markers_layer] = labels

    # Check nuclear labels.
    check_labels = adata.layers[expand_markers_layer].copy()
    check_label_types = np.unique(check_labels)[1:]

    delete_labels = []
    for label_type in check_label_types:
        check_label = check_labels.copy()
        check_label[check_label != label_type] = 0
        check_label[check_label == label_type] = 1
        contours = cv2.findContours(check_label.astype(np.uint8), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)[0]
        if len(contours) != 1:
            delete_labels.append(label_type)
    check_labels[np.isin(check_labels, delete_labels)] = 0
    adata.layers[expand_markers_layer] = check_labels

    # Expand nuclei labels to cytoplasm.
    if expand_method == "certain_distance":
        # Simply expand each label by some distance.
        st.cs.expand_labels(adata, expand_markers_layer, out_layer="cell_labels", **expand_kwargs)
    elif expand_method == "affinity":
        # we will take advantage of the fact that ssDNA staining weakly stains the cytoplasm,
        # and thus we can identify cytoplasmic regions by thresholding the image with a lenient threshold.
        st.cs.mask_cells_from_stain(adata, out_layer="stain_cell_mask", **expand_kwargs)
        st.cs.watershed(
            adata, "stain",
            mask_layer="stain_cell_mask",
            markers_layer=expand_markers_layer,
            out_layer="cell_labels",
        )
    else:
        pass

    # Expand all the cell labels by some amount to mitigate the effects of RNA diffusion.
    # The distance that should be expanded will depend on the level of RNA diffusion in the data.
    if expand_distance != 0:
        st.cs.expand_labels(adata, 'cell_labels', distance=expand_distance, out_layer='cell_labels_expanded')

    return adata


def image_enhancement(
    img: np.ndarray,
    radius=1.0,
    amount=1.0,
    block_size=3,
    offset=0,
    otsu_class=3,
    otsu_index=0,
    dilation_distance=2
):
    """Rong Xiang's image enhancemnet method.
    Args:
        img: Grayscale input image.
        radius: If a scalar is given, then its value is used for all dimensions. If sequence is given, then there must
                be exactly one radius for each dimension except the last dimension for multichannel images. Note that 0
                radius means no blurring, and negative values are not allowed.
        amount: The details will be amplified with this factor. The factor could be 0 or negative. Typically, it is a
                small positive number, e.g. 1.0.
        block_size: Image block processing with uneven brightness.
        offset: Try increase this value if too few are reserved.
        otsu_class: Classify image histograms into several categories. Bigger is slower.
        otsu_index: Pix less than this threshold will be filtered out. If the final image retains too much noise cell,
                    increase it. range:(0,otsu_class - 2).
        dilation_distance: If the cell segmentation retains too many noise cells, decrease it.
    """
    from skimage.morphology import dilation, disk
    from skimage.filters import threshold_multiotsu, threshold_local, unsharp_mask

    image = img.copy()
    # Unsharp masking is a linear image processing technique which sharpens the image.
    image = unsharp_mask(image, radius=radius, amount=amount)

    # Compute a threshold mask image based on local pixel neighborhood.
    local_thresh = threshold_local(image, block_size=block_size, method='gaussian', offset=offset)
    image[image < local_thresh] = 0

    # Generate classes-1 threshold values to divide gray levels in image, following Otsu’s method for multiple classes.
    thresholds = threshold_multiotsu(image, classes=otsu_class)
    image[image < thresholds[otsu_index]] = 0

    # Return grayscale morphological dilation of an image.
    image = dilation(image, disk(dilation_distance))

    return 

'''

def_path = '/hwfssz5/ST_SUPERCELLS/P21Z10200N0206/public/drosophila_pipe_v1/'
sys.path.append(def_path)
from supplement import crop_by_ssdna
from supplement import bin1tobinx


def qc_seg(
    project:str,
    seg_adata:AnnData,
    cropped_img:str,
    save=None,
):
    if not os.path.exists(save):
        os.makedirs(save)
    qc_save = os.path.join(save,'cellbin')
    if not os.path.exists(qc_save):
        os.makedirs(qc_save)
    qc_00_save = os.path.join(qc_save,project+'_qc_00_ssDNA.pdf')
    qc_01_save = os.path.join(qc_save,project+'_qc_01_segmentation.pdf')
    
    #qc_00
    st.pl.imshow(seg_adata,'stain',save_show_or_return="return")
    plt.savefig(qc_00_save, dpi=100)
    #qc_01
    n_layers = len(seg_adata.layers.keys())
    fig, axes = plt.subplots(ncols=n_layers+1, figsize=(18, 5), tight_layout=True)
    axes[0].imshow(plt.imread(cropped_img))
   # st.pl.imshow(seg_adata,'X', ax=axes[1],vmax = 10,vmin = 2,save_show_or_return="return")
    for i, layer in enumerate(seg_adata.layers.keys()):
        st.pl.imshow(seg_adata,'X', ax=axes[i+1],vmax = 10,vmin = 2, save_show_or_return="return")
        if i==0:
            label_=False
        else:
            label_=True
        st.pl.imshow(seg_adata, layer, ax=axes[i+1], labels=label_,alpha=0.6, save_show_or_return="return")
    plt.savefig(qc_01_save, dpi=100)
    print('===== ',project,'cellbin qc_01 saved! =====')

def qc_cellbin(
    project:str,
    adata:AnnData,
    save=None,
):
    if not os.path.exists(save):
        os.makedirs(save)
    qc_save = os.path.join(save,'cellbin')
    if not os.path.exists(qc_save):
        os.makedirs(qc_save)
    qc_02_save = os.path.join(qc_save,project+'_qc_02_violin.pdf')
    qc_03_save = os.path.join(qc_save,project+'_qc_03_scatter.pdf')
    qc_04_save = os.path.join(qc_save,project+'_qc_04_geo.pdf')
    qc_05_save = os.path.join(qc_save,project+'_qc_05_space.pdf')
    qc_txt_save = os.path.join(qc_save,'qc_cellbin.txt')
    
    #qc_02
    import scanpy as sc
    import numpy as np
    sc.pp.filter_cells(adata, min_genes=1)
    sc.pp.filter_genes(adata, min_cells=1)
    adata.obs['n_counts'] = adata.X.sum(axis=1)
    P1 = np.mean(adata.obs['n_genes'])
    Q1 = np.std(adata.obs['n_genes'],ddof = 1)
    P2 = np.mean(adata.obs['n_counts'])
    Q2 = np.std(adata.obs['n_counts'],ddof = 1)
    f = open(qc_txt_save, "a" ,newline='')
    n_cells = len(adata.obs)
    print('project:',project,' (cellbin)','\nn_cells:',n_cells,'\nn_genes:',P1,'±',Q1,'\nn_counts:',P2,'±',Q2,'\ndetails:','\n',adata,'\n',adata.obs[0:5],file = f)
    f.close()
    print('===== ',project,'cellbin qc_txt saved! =====')

    #qc_03
    #小提琴图
    sc.pl.violin(adata, ['n_genes', 'n_counts',],jitter=0.4, multi_panel=True)
    plt.savefig(qc_02_save, dpi=100)
    print('===== ',project,'cellbin qc_02 saved! =====')
    sc.pl.scatter(adata, x='n_counts', y='n_genes')
    plt.savefig(qc_03_save, dpi=100)
    print('===== ',project,'cellbin qc_03 saved! =====')

    #qc_04
    #面积图
    # cell bin visualization - geo
    adata.obs["bigger"] = "The cell smaller than bin20"
    adata.obs.loc[adata.obs["area"] >= 400, "bigger"] = "The cell larger than bin20"

    st.pl.geo(
        adata, genes=['area'], color=['bigger'], ncols=2, show_legend="upper left",
        figsize=(8, 8), dpi=300, save_show_or_return="save", boundary_width=0,
        save_kwargs={"path": qc_04_save,"prefix": "geo", "dpi": 300, "ext": "pdf"}
        )
    print('===== ',project,'cellbin qc_04 saved! =====')

    # cell bin visualization - point
    st.pl.space(
        adata=adata, genes=["area"], color=['bigger'], pointsize=0.1, dpi=300, 
        boundary_width=0,show_legend="upper left", ncols=2, save_show_or_return="save",
        save_kwargs={"path": qc_05_save,"prefix": "point", "dpi": 300, "ext": "pdf"}
        )

    print('===== ',project,'cellbin qc_05 saved! =====')


def qc_bin20(
    project:str,
    adata:AnnData,
    save=None,
):
    if not os.path.exists(save):
        os.makedirs(save)
    bin20_out_qc = os.path.join(save,'bin20')
    if not os.path.exists(bin20_out_qc):
        os.makedirs(bin20_out_qc)
    bin20_qc_txt_save = os.path.join(bin20_out_qc,'qc_bin20.txt')
    bin20_qc_02_save = os.path.join(bin20_out_qc,project+'_qc_02_violin_bin20.pdf')
    bin20_qc_03_save = os.path.join(bin20_out_qc,project+'_qc_03_scatter_bin20.pdf')
    #qc_02
    import scanpy as sc
    import numpy as np
    sc.pp.filter_cells(adata, min_genes=1)
    sc.pp.filter_genes(adata, min_cells=1)
    adata.obs['n_counts'] = adata.X.sum(axis=1)
    P1 = np.mean(adata.obs['n_genes'])
    Q1 = np.std(adata.obs['n_genes'],ddof = 1)
    P2 = np.mean(adata.obs['n_counts'])
    Q2 = np.std(adata.obs['n_counts'],ddof = 1)
    f = open(bin20_qc_txt_save, "a" ,newline='')
    n_cells = len(adata.obs)
    print('project:',project,' (bin20)','\nn_cells:',n_cells,'\nn_genes:',P1,'±',Q1,'\nn_counts:',P2,'±',Q2,'\ndetails:','\n',adata,'\n',adata.obs[0:5],file = f)
    f.close()
    print('===== ',project,'bin20 qc_txt saved! =====')

    #qc_03
    #小提琴图
    sc.pl.violin(adata, ['n_genes', 'n_counts',],jitter=0.4, multi_panel=True)
    plt.savefig(bin20_qc_02_save, dpi=100)
    print('===== ',project,'bin20 qc_02 saved! =====')
    sc.pl.scatter(adata, x='n_counts', y='n_genes')
    plt.savefig(bin20_qc_03_save, dpi=100)
    print('===== ',project,'bin20 qc_03 saved! =====')
