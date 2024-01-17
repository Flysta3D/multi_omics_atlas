#https://github.com/bgi-drosophila/spatial_pipeline/tree/main/cell_segmentation

from pickle import TRUE
import cv2
import sys
sys.path.append("/hwfssz5/ST_SUPERCELLS/P21Z10200N0206/xiangrong/software/miniconda3/envs/spateo0513/lib/python3.8/site-packages/")
from cellpose import models
import spateo as st
from skimage.morphology import dilation,disk 
from skimage.filters import  threshold_multiotsu, threshold_local,unsharp_mask
from anndata import AnnData
from typing import Literal, Optional, Union
import numpy as np
def cellpose(
	img: np.ndarray,
	model_type: Literal["cyto", "nuclei"] = "nuclei",
	**kwargs,
) -> np.ndarray:
	"""Run Cellpose on the provided image.
	Args:
		img: Image as a Numpy array.
		model_type: Cellpose model to use. Can be one of the two pretrained models:
			* cyto: Labeled cytoplasm
			* nuclei: Labeled nuclei
			Or any generic CellposeModel model.
		**kwargs: Additional keyword arguments to :func:`Cellpose.eval`
			function.
			
	Returns:
		Numpy array containing cell labels.
	"""
	model = models.Cellpose(model_type=model_type)
	res,_, _,_ = model.eval(img,channels=[0, 0],**kwargs)
	return res

def custom_segmentation(
	path_bin1: str,
	path_ssdna: str,
	radius: Optional[int] = None,
	amount: int = 10,
	block_size: int = 9,
	offset: int = 0,
    otsu_class: int = 3,
    otsu_index: int = 1,
	dilation_distance: int = 2,
	diameter:Optional[int] = None,
	flow_threshold: float = 0.8,
	min_size: int = -1,
	**kwargs,
)-> AnnData:
	""" Run custom_segmentation
	Args:
		path_bin1: Bin1 lasso file start with coordinate 0
		path_ssdna: Cropped ssDNA
		radius: Default:None
		amount: Recommend 5-30,The higher the value, the brighter the image(overexposure).Default:10
		block_size: Image block processing with uneven brightness. Default:9
		offset: Try increase this value if too few are reserved. Default:0 
		otsu_class: Classify image histograms into several categories. Bigger is slower.Default:3
		otsu_index: Pix less than this threshold will be filtered out. If the final image retains too much noise cell, increase it. range:(0,otsu_class - 2). Default:1
		dilation_distance: If the cell segmentation retains too many noise cells, decrease it. Default:2
		diameter: (float (optional, default 1)), if set to None, then diameter is automatically estimated if size model is loaded
		flow_threshold:  (float (optional, default 0.8)), flow error threshold (all cells with errors below threshold are kept) (not used for 3D)
		min_size:(int (optional, default -1)), minimum number of pixels per mask, can turn off with -1
		**kwargs: Additional keyword arguments to :func:`Cellpose.eval` 
			function.https://cellpose.readthedocs.io/en/latest/api.html#cellposemodel
	"""
	adata = st.io.read_bgi_agg(path_bin1, path_ssdna)
	img = cv2.imread(path_ssdna,0)
	## 增强,图像够亮不需要增强
	if radius != None:
		image = unsharp_mask(img, radius = radius, amount= amount)
	else:
		image = img
	## 局部过滤
	local_thresh = threshold_local(image, block_size, offset = offset)
	image[image < local_thresh] = 0
	##全局过滤
	thresholds = threshold_multiotsu(image, classes = otsu_class)
	image[image<thresholds[otsu_index]] = 0
	##膨胀
	if dilation_distance != 0:
		footprint = disk(dilation_distance) 
		image = dilation(image, footprint)
	##分割
	segmented_cells = cellpose(img = image,
							   diameter = diameter,
							   flow_threshold = flow_threshold,
							   min_size = min_size,
							   **kwargs
							   )
	## add
	adata.layers['stain_mask'] = image
	# Check nuclear labels.
	check_labels = segmented_cells
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
	adata.layers['cellpose_labels'] = check_labels
	return adata


###qc & calculation distance

def distance_to_center(seg_adata = None,bin1_df_path = None,save_img_path = None,layer = "cell_labels"):
	#seg_adata = st.read_h5ad("/hwfssz5/ST_SUPERCELLS/P21Z10200N0206/project/drosophila/03.ST/02.segmentation/E18-20h_b/01_segmentation/seg_adata.h5ad")
	cells = pd.DataFrame(seg_adata.layers[layer])
	cells["x"] = cells.index
	cells = pd.melt(cells, id_vars=["x"])
	cells.columns = ["x", "y", layer]
	cells = cells[cells[layer] != 0]
	bin1_df = pd.read_csv(
		bin1_df_path,
		sep="\t",
		dtype={
			"geneID": "category",  # geneID
				"x": np.uint32,  # x
				"y": np.uint32,  # y
				"MIDCounts": np.uint16,  # total
			},
			comment="#",
		)
	cell_bin1_df = pd.merge(cells, bin1_df, on=["x", "y"], how="left")
	cell_bin1_df.dropna(inplace=True) 
	#cell_i = cell_bin1_df[layer].unique()[0]
	cell_bin1_df = cell_bin1_df.sort_values(layer,ascending =True)
	distance_list = []
	for cell_i in cell_bin1_df[layer].unique():
		tdf = cell_bin1_df.loc[cell_bin1_df[layer] ==cell_i,].reset_index()
		min_x = tdf.x.min()
		max_x = tdf.x.max()
		min_y = tdf.y.min()
		max_y = tdf.y.max()
		center_x = (max_x - min_x)/2+min_x
		center_y = (max_y - min_y)/2+min_y
		for i in tdf.index:
			x,y = tdf['x'][i],tdf['y'][i]
			distance_to_center = math.sqrt(((x-center_x)**2)+((y-center_y)**2))
			distance_list.append(distance_to_center)
	cell_bin1_df["distance_to_cellcenter"] = distance_list
	# nuclear
	nuclear_gene = ["cas","grk","CG33941","Max","lncRNA:CR42646","lncRNA:CR31044","asRNA:CR42871","lncRNA:CR32690","lncRNA:CR42862","asRNA:CR43464","asRNA:CR43476",
	"lncRNA:CR30121","lncRNA:CR42651","asRNA:CR44052","asRNA:CR43872","lncRNA:CR40053","asRNA:CR44872","lncRNA:CR44931","lncRNA:flam","lncRNA:iab4",
	"lncRNA:bxd","lncRNA:roX1","lncRNA:Hsromega","lncRNA:CR33963","lncRNA:CR42549","CR31032"]
	cell_bin1_df1 = cell_bin1_df.loc[cell_bin1_df["geneID"].isin(nuclear_gene),].reset_index()
	#cell_bin1_df1 = cell_bin1_df.loc[cell_bin1_df["geneID"].str.contains("^snoRNA\:"),].reset_index()
	print(cell_bin1_df1["geneID"].value_counts())
	## consider counts
	distance_to_cellcenter = [[str(cell_bin1_df1["distance_to_cellcenter"][i])]*int(cell_bin1_df1["MIDCounts"][i]) for i in cell_bin1_df1.index]
	## flatten
	distance_to_cellcenter = [float('%.2f' % float(element)) for sublist in distance_to_cellcenter for element in sublist]
	df_nuclear = pd.DataFrame({"value":distance_to_cellcenter})
	df_nuclear["type"] = "Nuclear RNA"
	## mt
	cell_bin1_df1 = cell_bin1_df.loc[cell_bin1_df["geneID"].str.contains("^mt\:"),].reset_index()
   # print(cell_bin1_df1["geneID"].value_counts())
	## consider counts
	distance_to_cellcenter = [[str(cell_bin1_df1["distance_to_cellcenter"][i])]*int(cell_bin1_df1["MIDCounts"][i]) for i in cell_bin1_df1.index]
	## flatten
	distance_to_cellcenter = [float('%.2f' % float(element)) for sublist in distance_to_cellcenter for element in sublist]
	df_mt = pd.DataFrame({"value":distance_to_cellcenter})
	df_mt["type"] = "Mitochondria RNA"
	xlim = 30 if re.search("random",layer) == None else 60
	fig = plt.figure() 
	fig = plt.xlim(0, int(xlim))
	fig = sns.distplot(df_nuclear["value"], hist=True,label ="Nuclear RNA",color = "green",bins=50 )
	fig = sns.distplot(df_mt["value"], hist=True,label ="Mitochondria RNA",color = "red",bins=50)
	plt.legend()
	fig.figure.savefig(save_img_path)
	return None


def qc(seg_adata,folder_cellbin_labels,bin_folder, bin1_file,folder_seg_image,folder_seg_h5ad,method = "cell_labels"):
	cells = pd.DataFrame(seg_adata.layers[method])
	cells["x"] = cells.index
	cells = pd.melt(cells, id_vars=["x"])
	cells = cells[cells["value"] != 0]
	cells.columns = ["x", "y", method]
	cells.to_csv(os.path.join(folder_cellbin_labels, f"{bin1_file[:-7]}_{method}.gem.gz"), sep="\t", index=False)
	## plot nuclear and mt gene distribution
	save_img_nuclear_mt = os.path.join(folder_seg_image, f"{method}_nuclear_mt_gene_distribution.pdf")
	distance_to_center(seg_adata = seg_adata,bin1_df_path = os.path.join(bin_folder, bin1_file),layer = method,save_img_path = save_img_nuclear_mt)
	### 
	adata = st.io.read_bgi(
		path=os.path.join(bin_folder, bin1_file),
		segmentation_adata=seg_adata,
		labels_layer=method,
		seg_binsize=1,
	)
	adata.obs["slices"] = str(bin1_file[:-7])
	basic_stats(adata=adata, mito_label="mt:")
	adata.write(os.path.join(folder_seg_h5ad, f"{bin1_file[:-7]}_{method}.h5ad"), compression="gzip")
	# cell bin visualization - geo	
	## for mouse bin30
	adata.obs["bigger"] = "The cell smaller than bin30"
	adata.obs.loc[adata.obs["area"] >= 900, "bigger"] = "The cell larger than bin30"
	# cell bin visualization - point
	st.pl.space(adata=adata, genes=["area"], color=['bigger'], pointsize=0.1, dpi=300, boundary_width=0,
				show_legend="upper left", ncols=2, save_show_or_return="save",
				save_kwargs={"path":os.path.join(folder_seg_image, f"{bin1_file[:-7]}_{method}"),
							"prefix": "point", "dpi": 300, "ext": "pdf"})
	# qc plot
	## 
	# basic_stats before qc
	basic_stats_multi(adatas=adata, slice_key="slices", inner="box", save_show_or_return="save",
				  save_kwargs={"path": os.path.join(folder_seg_image, f"{bin1_file[:-7]}_{method}"),
				   "prefix": 'basic_stats', "dpi": 300, "ext": 'pdf'})
	print(method)
	print(np.mean(adata.obs['nGenes']))


# path_bin1 = "/hwfssz5/ST_SUPERCELLS/P21Z10200N0206/project/drosophila/03.ST/02.segmentation/E18-20h_b/01_segmentation/bin1_data/E18-20h_b_S9_bin1_final.gem.gz"
# path_ssdna = "/hwfssz5/ST_SUPERCELLS/P21Z10200N0206/project/drosophila/03.ST/02.segmentation/E18-20h_b/01_segmentation/cropped_img/E18-20h_b_S9_cropped_img.tif"
# id = "s9"

# seg_adata = custom_segmentation(path_bin1 = path_bin1,path_ssdna = path_ssdna,radius = 2)
# fig, axes = plt.subplots(ncols=4, figsize=(18, 6), tight_layout=True)
# #axes[0].imshow(plt.imread(path_ssdna))
# st.pl.imshow(seg_adata, 'stain', ax=axes[0],cmap='gray')
# st.pl.imshow(seg_adata, "stain_mask", ax=axes[1],cmap='gray',alpha=1)
# st.pl.imshow(seg_adata, 'cellpose_labels',labels=True, alpha=1, ax=axes[2])
# st.pl.imshow(seg_adata, 'X',vmax = 10,vmin = 2 alpha=1, ax=axes[3])
# plt.savefig(os.path.join(folder_seg_image, f"{id}_model_1.png"), dpi=300)