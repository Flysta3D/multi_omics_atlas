#before_PS
import os,re,sys,getopt

#input
def getargv(argv):
    try:
        opts, args = getopt.getopt(argv,"-h-l:-o:",["help","lasso_from_web_folder=","save_folder="])
    except getopt.GetoptError:
        print("01_drosophila_pipe_v1_before_PS.py -l <lasso_from_web_folder>  -o <save_folder>")
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print("01_drosophila_pipe_v1_before_PS.py -l <lasso_from_web_folder>  -o <save_folder>")
            sys.exit()
        elif opt in ("-l", "--lasso_from_web_folder"):
            lasso_path = arg
        elif opt in ("-o", "--save_folder"):
            save_folder = arg 
    return lasso_path,save_folder

indir,outdir = getargv(sys.argv[1:])


if not os.path.exists(outdir):
    os.makedirs(outdir)

import spateo as st
import matplotlib.pyplot as plt
from PIL import Image

k = os.listdir(indir)
k.sort()

for i in k:
    project = i[:-7]
    raw_bin1_path = indir+'/'+i
    print('======',project,'start ... ======')
    
    #output
    bin1_img_save = outdir+'/'+project+'_bin1_beforePS.tiff'
    adata = st.io.read_bgi_agg(raw_bin1_path)

    img = adata.X.todense()
    img[img>=10]=10
    img = img/10*255
    img = Image.fromarray(img)
    img = img.convert("L")

    img.save(bin1_img_save)
