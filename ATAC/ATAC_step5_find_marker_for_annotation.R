library(Seurat)
library(Signac)
library(ggplot2)

##########################Create gene activity matrix###########################
data_dir <- "/workdir"
df <- readRDS(paste0(data_dir,"/",file))
gene.activities <- GeneActivity(df)
df[['RNA']] <- CreateAssayObject(counts = gene.activities)
df <- NormalizeData(df)
saveRDS(paste0(data_dir,"/Gene_activity_whole_embeding.Rds"))

######################find marker for tissue annotation#########################
res <- 0.3
df <- FindNeighbors(df, reduction = "pca.l2", nn.eps=0.5, dims=2:50)
df <- FindClusters(df,n.start=20, resolution=res)
DE.genes <- FindAllMarkers(df, assay = 'RNA', min.pct = 0.1,
                           logfc.threshold = 0.1)
write.table(DE.gene,paste0(data_dir,"/Marker_genes_for_tissue_anno.tsv"),
            quote = F,row.names = T,col.names = T,sep = "\t")

##########################add tissue annotation#################################
tissue.anno <- read.delim(paste0(data_dir,"/tissue_anno.tsv"))
df$tissue <- as.numeric(as.character(df$seurat_clusters))
for (i in 0:unique(df$tissue)) {
    df$tissue[which(df$tissue == i)] <- tissue.anno$tissue[which(tissue.anno$cluster == i)]
}

#########################find marker for sub-clusters###########################
res <- 0.1
#res <- 0.9
dir.create(paste0(data_dir,"/total/sub_clusters/res",res))
clusters <-  unique(df$tissue)
color <- c("#ea5545", "#27aeef", "#ef9b20", "#b33dc6", "#edbf33", "#87bc45", 
           "#ede15b", "#1a53ff", "#bdcf32", "#f46a9b",  "#7c1158", "#b30000",
           "#4421af", "#0d88e6", "#00b7c7", "#5ad45a", "#8be04e", "#ebdc78", 
           "#F08080FF","#ef92b5" ,"#6495EDFF", "#2196F3FF" ,"#29B6F6FF", 
           "#87CEEBFF" ,"#8CBC68FF" ,"#A6BE54FF","#BEBC48FF" ,"#D1B541FF" ,
           "#DDAA3CFF" ,"#E49C39FF" ,"#E78C35FF" , "#E4632DFF" ,"#DF4828FF" ,
           "#DA2222FF","#ecd452","#ffee6f"  , "#FEE08B" ,  "#E6F598","#b32142" ,
           "#d12820","#ea5514","#f88922","#DDA0DDFF", "#FF69B4FF" ,"#5F9EA0FF", 
           "#FFDAB9FF", "#FFA07AFF","#721E17FF" ,"#521A13FF","#DDD8EFFF", 
           "#D1C1E1FF" ,"#C3A8D1FF","#A778B4FF","#8C4E99FF" ,"#9C27B0FF",
           "#BA55D3FF" , "#6F4C9BFF" ,"#6059A9FF" ,"#5568B8FF" ,"#4E96BCFF",
           "#59A5A9FF","#69B190FF" ,"#77B77DFF", "#d3a237","#B8221EFF",
           "#95211BFF" , "#9E0142"  , "#ABDDA4" ,  "#66C2A5" ,"#006064FF",  
           "#4CAF50FF" , "#FFEB3BFF", "#FF9800FF","#795548FF", "#9E9E9EFF", 
           "#607D8BFF","#32CD32FF"  ,"#4682B4FF" ,"#9ACD32FF" ,"#40E0D0FF",
           "#F0E68CFF" ,"#D2B48CFF","#c0d695" )
lapply(clusters,function(clusters){
  df.sub <- subset(df, tissue == clusters)
  df.sub <- RunUMAP(df.sub, reduction = 'pca.l2', dims = 2:50)
  df.sub <- FindNeighbors(df.sub, reduction = "pca.l2", nn.eps=0.5, dims=2:50)
  df.sub <- FindClusters(df.sub,n.start=20, resolution=res)
  umap = df.sub@reductions$umap@cell.embeddings %>% 
    as.data.frame()
  ext_len1 <- (max(umap$UMAP_1)-min(umap$UMAP_1))/10
  ext_len2 <- (max(umap$UMAP_2)-min(umap$UMAP_2))/10
  p <- DimPlot(df.sub, reduction = "umap", group.by = "seurat_clusters", 
               cols = color, raster=FALSE)+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.border = element_blank(), 
          axis.title = element_blank(),  
          axis.text = element_blank(), 
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          plot.title= element_blank(),
          panel.background = element_rect(fill = 'white'), 
          plot.background=element_rect(fill="white"))
  ggsave(p, filename = paste0(clusters,"_sub_clusters_res",res,".pdf"), 
         path = paste0(data_dir,"/total/sub_clusters/res",res),width = 10, 
         height = 7, units = "in",dpi = 300)
  DE.genes <- FindAllMarkers(df.sub,assay = 'RNA', min.pct = 0.1, 
                             logfc.threshold = 0.1)
  write.table(DE.genes,paste0(data_dir,"/total/sub_clusters/res",res,"/",
                              clusters,"_res",res,"_DE_genes.tsv"), sep = "\t")
})

