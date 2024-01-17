
import os,re,sys,getopt

#input
def getargv(argv):
    try:
        opts, args = getopt.getopt(argv,"-h-l:-t:-o:",["help","lasso_from_web_folder=","afterPS_tif_folder=","save_folder="])
    except getopt.GetoptError:
        print("02_drosophila_pipe_v1_after_PS.py -l <lasso_from_web_folder> -t <afterPS_tif_folder> -o <save_folder>")
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print("02_drosophila_pipe_v1_after_PS.py -l <lasso_from_web_folder> -t <afterPS_tif_folder> -o <save_folder>")
            sys.exit()
        elif opt in ("-l", "--lasso_from_web_folder"):
            lasso_path = arg
        elif opt in ("-t", "--afterPS_tif_folder"):
            tifdir_path = arg
        elif opt in ("-o", "--save_folder"):
            save_folder = arg 
    return lasso_path,tifdir_path,save_folder

lasso_path,tifdir_path,aout = getargv(sys.argv[1:])

#
def_path = '../'
sys.path.append(def_path)
from supplement import crop_by_ssdna,bin1tobinx,custom_segmentation,cellpose,distance_to_center
from supplement import qc_seg,qc_cellbin,qc_bin20,equalhist_transfer,output_img
from supplement import check_nan_gene
from supplement import deHaze, output_img, segmentation_spateo

#aout
name_list = ['cropped_img','bin1_data','seg_02','bin20']
for aaa in name_list:
    out = os.path.join(aout,aaa)
    if not os.path.exists(out):
        os.makedirs(out)

out_02 = os.path.join(aout,'seg_02')

name_l =['qc','seg_data','seg_gem']
for aaa in name_l:
    out = os.path.join(out_01,aaa)
    if not os.path.exists(out):
        os.makedirs(out)
    out = os.path.join(out_02,aaa)
    if not os.path.exists(out):
        os.makedirs(out)

import cv2
import numpy as np
import warnings
from pathlib import Path
from typing import Literal, Optional, Union,Tuple
import spateo as st
from anndata import AnnData
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")
import pandas as pd
import anndata as ad
import torch
from skimage.filters import unsharp_mask
from matplotlib.backends.backend_pdf import PdfPages

def seg(
    k,
    tifdir_path:str,
    aout:str,
    min_txt_save:str,
    tiflist:str,
    pdf_save:str
    ):
    nrow = len(tif_list)
    with PdfPages(pdf_save) as pdf:
        fig1 = plt.figure()
        fig1, axes1 = plt.subplots(nrows=nrow,ncols=3, figsize=(18, 5*nrow), tight_layout=True)
        fig2 = plt.figure()
        fig2, axes2 = plt.subplots(nrows=nrow,ncols=5, figsize=(18, 5*nrow), tight_layout=True)
        fig3 = plt.figure()
        fig3, axes3 = plt.subplots(nrows=nrow,ncols=5, figsize=(18, 5*nrow), tight_layout=True)
        row = 0
        for i in tiflist:
            # get slice name : E8-10h_a_S01
            project = '_'.join(re.split(r'[_]\s*',i[:-4])[:3])
            m = re.findall(r'(SS.+?).tif',i)[0]
            while (m in k.keys()) == False:
                m = '_'.join(re.split(r'[_]\s*',m)[:len(re.split(r'[_]\s*',m))-1])
                if len(re.split(r'[_]\s*',m))==1:
                    print('tif name error!!')
                    break

            # raw bin1 of each part
            raw_bin1_path = k[m]

            # tif after PS
            img_after_PS_path = os.path.join(tifdir_path,i)
            
            print('===== ',project,'start =====')
            print('===== raw_bin1_path:',raw_bin1_path,' =====')
            print('===== img_after_PS_path:',img_after_PS_path,' =====')
            
            # data_output
            cropped_img_save = os.path.join(aout,'cropped_img',project+'_cropped_img.tif')
            new_bin1_save = os.path.join(aout,'bin1_data',project+'_bin1_final.gem.gz')
            
            # 01_crop_by_ssdna
            cropped_data, _ = crop_by_ssdna(
                path=raw_bin1_path, 
                path_ssdna=img_after_PS_path, 
                save_lasso=None, 
                save_img=cropped_img_save, 
                show=False
            )
            print('===== ',project,'cellbin cropped data finished! =====')

            # 02_min data
            ## make x&y from 0 to np.inf
            min_data = cropped_data.copy()
            a = min_data["x"].min()
            b = min_data["y"].min()
            min_data["x"] = min_data["x"] - a
            min_data["x"] = min_data["x"].astype(int)
            min_data["y"] = min_data["y"] - b
            min_data["y"] = min_data["y"].astype(int)
            min_data.to_csv(new_bin1_save, sep="\t", index=False)
            ## check nan
            data = check_nan_gene(
                path=new_bin1_save, nan_gene="CG5842"
            )
            data.to_csv(new_bin1_save, sep="\t", index=False)
   
            f = open(min_txt_save, "a" ,newline='')
            print([project,a,b],file = f)
            f.close()
            print('===== ',project,'cellbin new bin1 saved! =====')

            #####################################################
            ######################  bin20  ######################
            #####################################################
            
            tout = os.path.join(aout,'bin20','bin20_gem')
            if not os.path.exists(tout):
                os.makedirs(tout) 
            
            tout = os.path.join(aout,'bin20','bin20_data')
            if not os.path.exists(tout):
                os.makedirs(tout)
            
            tout = os.path.join(aout,'bin20','bin20_qc')
            if not os.path.exists(tout):
                os.makedirs(tout)
            
            bin20_gem_save = os.path.join(aout,'bin20','bin20_gem',project+'_bin20.gem.gz')
            bin20_data_save = os.path.join(aout,'bin20','bin20_data',project+'_bin20.h5ad')
           
            seg_out_qc_bin20 = os.path.join(aout,'bin20','bin20_qc')
            
            bin1_data = pd.read_csv(new_bin1_save, sep="\t")
            bin20_data = bin1tobinx(
                            bin1_data=bin1_data, 
                            binx=int(20), 
                            save=bin20_gem_save
                        )
            adata = st.io.read_bgi(
                        path=bin20_gem_save, 
                        binsize=int(20)
                    )

            adata.obs["slices"] = str(project)
            adata.uns['raw_min'] = {project+'_x':a,project+'_y':b}

            adata.write(bin20_data_save, compression="gzip")
            print('===== ',project,' bin20 data saved! =====')
            
            # qc_bin20
            qc_bin20(
                project = project,
                adata = adata,
                save = seg_out_qc_bin20,
            )

            print('===== ',project,' bin20 finished! =====')

            # data_output
            out_02 = os.path.join(aout,'seg_02')
            seg_data_save_02 = os.path.join(out_02,'seg_data',project+'_seg.h5ad')
            enhance_ssdna_image2_folder = os.path.join(out_02,'img_enhance')
            enhance_ssdna_image2_save = os.path.join(enhance_ssdna_image2_folder,project+'_enhanced_img.tif')
            if not os.path.exists(enhance_ssdna_image2_folder):
                os.makedirs(enhance_ssdna_image2_folder)
            seg_out_qc_02 = os.path.join(out_02,'qc')
            #qc_txt_save_02 = os.path.join(out_02,'qc','all_slices_qc.txt')
            seg_gem_save_02 = os.path.join(out_02,'seg_gem',project+'_seg.gem.gz')
            
            #image_enhance
            raw_ssdna_image = cv2.imread(cropped_img_save)
            enhance_ssdna_image1 = deHaze(src=raw_ssdna_image.copy(), min_r=5, guide_r=20)
            enhance_ssdna_image1 = cv2.cvtColor(
                enhance_ssdna_image1.astype(np.uint8), cv2.COLOR_RGB2GRAY
            )
            enhance_ssdna_image2 = equalhist_transfer(
                img=enhance_ssdna_image1.astype(np.uint8),
                method="local",
                cliplimit=7,
                tileGridSize=(10, 10),
            )
            output_img(
                img=enhance_ssdna_image2,
                filename=enhance_ssdna_image2_save,
                show_img=False
            )

            # Segmentation
            seg_adata = segmentation_spateo(
                path_bin1=new_bin1_save,
                path_ssdna=enhance_ssdna_image2_save,
                seg_method="cellpose",
                seg_kwargs=dict(
                    model="nuclei",
                    diameter=None,
                    normalize=True,
                    flow_threshold=1,
                    min_size=-1,
                    device=torch.device("cpu"),
                ),
                expand_method="certain_distance",
                expand_kwargs=dict(distance=2, max_area=400),
            )
            seg_adata.layers["original_ssDNA"] = raw_ssdna_image
            seg_adata.layers["enhanced_ssDNA1"] = enhance_ssdna_image1
            seg_adata.layers["enhanced_ssDNA2"] = enhance_ssdna_image2
            
            #qc
            #st.pl.imshow(seg_adata,'X', ax=axes1[int(row),1],vmax = 10,vmin =2,save_show_or_return="return")
            st.pl.imshow(seg_adata,'original_ssDNA', ax=axes1[int(row),1],labels=False,save_show_or_return="return")
            st.pl.imshow(seg_adata,'X', ax=axes1[int(row),2],vmax = 10,vmin =2,save_show_or_return="return")
            st.pl.imshow(seg_adata,'cellpose_labels', ax=axes1[int(row),2],labels=True,alpha=0.6, save_show_or_return="return")
            ax = axes1[int(row),2]
            ax.set_title(project+'_seg_02')
            
            # qc_seg
            st.pl.imshow(seg_adata,"original_ssDNA",ax=axes3[int(row),0],labels=False,save_show_or_return="return",cmap="gray",)
            st.pl.imshow(seg_adata,"enhanced_ssDNA1",ax=axes3[int(row),1],labels=False,save_show_or_return="return",cmap="gray",)
            st.pl.imshow(seg_adata,"enhanced_ssDNA2",ax=axes3[int(row),2],labels=False,save_show_or_return="return",cmap="gray",)
            st.pl.imshow(seg_adata,'X', ax=axes3[int(row),3],vmax = 10,vmin =2,save_show_or_return="return")
            st.pl.imshow(seg_adata,"cellpose_labels",ax=axes3[int(row),3],labels=True,alpha=0.6,save_show_or_return="return",)
            st.pl.imshow(seg_adata,'X', ax=axes3[int(row),4],vmax = 10,vmin =2,save_show_or_return="return")
            st.pl.imshow(seg_adata,"cell_labels",ax=axes3[int(row),4],labels=True,alpha=0.6,save_show_or_return="return",)
            
            ax = axes3[int(row),0]
            ax.set_title(project)
            
            #####################################################
            ###################### cellbin ######################
            #####################################################
            cell_layer = "cell_labels"
            # cell labels
            cells = pd.DataFrame(seg_adata.layers[cell_layer])
            cells["x"] = cells.index
            cells = pd.melt(cells, id_vars=["x"])
            cells = cells[cells["value"] != 0]
            cells.columns = ["x", "y", cell_layer]
            cells.sort_values(by=[cell_layer, "x", "y"], inplace=True)
            cells.to_csv(seg_gem_save_02,sep="\t",index=False,)

            # cell bin
            adata = st.io.read_bgi(
                path=new_bin1_save,
                segmentation_adata=seg_adata,
                labels_layer=cell_layer,
                seg_binsize=1,
            )
            adata.obs["slices"] = str(project)
            adata.uns['raw_min'] = {project+'_x':a,project+'_y':b}
            adata.write(seg_data_save_02, compression="gzip")
            
            # qc_cellbin
            qc_cellbin(
                project = project,
                adata = adata,
                save = seg_out_qc_02,
            )

            row = row+1
        
        pdf.savefig(fig1)
        plt.close()
        pdf.savefig(fig2)
        plt.close()
        pdf.savefig(fig3)
        plt.close()
        print('===== all slices segmentation finished ! =====')


#### start ####

min_txt_save = os.path.join(aout,'min_x_y.txt')
k = {}
for root,dirs,files in os.walk(lasso_path):
    for file in files:
        k.update({file[:-7]:os.path.join(root,file)})

tif_list = os.listdir(tifdir_path)
tif_list.sort()

num = len(tif_list)

#01
tif_list1 = tif_list
comparison_save_1 = os.path.join(aout,'comparison_of_all_01.pdf')

seg(k = k,
    tifdir_path = tifdir_path,
    aout = aout,
    min_txt_save = min_txt_save,
    tiflist=tif_list1,
    pdf_save = comparison_save_1)

