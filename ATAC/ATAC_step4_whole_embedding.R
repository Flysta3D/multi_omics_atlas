library(harmony)
library(Seurat)
library(Matrix)
library(Signac)
library(GenomeInfoDb)
library(annotables)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(gridExtra)
library(reticulate) 
library(parallel)

############################necessary information###############################
data_dir <- '/workdir'
bdgp_df_subset = bdgp6[bdgp6$chr %in% c('2R','2L','3L','3R','4'),]
bdgp_sub <- bdgp_df_subset %>% mutate(strand=case_when(strand==1~"+",strand==-1~"-"))
colnames(bdgp_sub)[colnames(bdgp_sub) %in% c("ensgene","symbol","biotype")]=c("gene_id","gene_name","gene_biotype")
annotations <- makeGRangesFromDataFrame(bdgp_sub,strand.field = "strand",keep.extra.columns = T) 
seqlevelsStyle(annotations) <- 'NCBI'
genome(annotations) <- "dm6"
timewindow_dirs <-  list.dirs(path = data_dir, recursive = F)
timewindow_dirs <- timewindow_dirs[grep("E",timewindow_dirs)]
timewindows <- str_remove(timewindow_dirs,pattern = paste0(data_dir,"/"))


#######################first round merge of chromatin assays####################  
count <- list()
for(i in 1:length(timewindows)){
    timewindow_dir <- timewindow_dirs[i]
    timewindow <- timewindows[i]
    lanes <- list.dirs(path = timewindow_dir, recursive = F)
    lane.list <- str_remove(lanes,pattern = paste0(timewindow_dir,"/"))
    lane.list <- sapply(strsplit(lane.list,"_"),"[[",1)
    exclude <- c("analysis","merge-peaks-idr")
    excludedir <- sapply(exclude, function(exclude){paste0(timewindow_dir,"/",exclude)})
    lane.list <- lane.list[which(!(lane.list %in% exclude))]
    lanes <- lanes[which(!(lanes %in% excludedir))]
    clnum<-20
    cl <- makeCluster(getOption("cl.cores", clnum),type="FORK")
    clusterExport(cl, varlist=c("timewindow_dir","annotations","data_dir"),envir=environment())
    largelist <- parLapply(cl,lanes,function(lanes){
        lane.id <- strsplit(sub(paste0(timewindow_dir,"/"),"" ,lanes),"_")[[1]][1]
        barcode.path <- paste0(timewindow_dir,"/analysis_new/03.peak_by_clusters_and_reclustering/",lane.id,"_kept_cells.tsv")
        peak.path <- paste0(data_dir,"/total.new.merged.bed")
        peak.file <- read.delim(peak.path, sep = "\t" ,header = F, stringsAsFactors = F)
        colnames(peak.file) <- c("seqnames","start","end")
        peak_ranges <- makeGRangesFromDataFrame(peak.file)
        barcode.names <- read.delim(barcode.path,header = F, stringsAsFactors = F)
        fragment.name <- list.files(path = lanes, pattern = "fragments.tsv.gz$")
        fragment.path <- paste(lanes,fragment.name, sep = "/")
        kept_cell <- barcode.names$V1
        kept_cell <- sub(paste0(lane.id,"_"),"",kept_cell)
        fragment.object <- CreateFragmentObject(fragment.path, cells = kept_cell,validate.fragments = T)
        peak_matrix <- FeatureMatrix(
            fragments = fragment.object,
            features = peak_ranges
        )
        count.small <- CreateChromatinAssay(
            counts = peak_matrix,
            fragments = fragment.path
        )
        metadata <- read.delim(paste0(timewindow_dir,"/analysis_new/03.peak_by_clusters_and_reclustering/",lane.id,"_metadata.tsv"), header = TRUE, row.names = 1, stringsAsFactors = F)
        metadata <- as.data.frame(metadata)
        count.small <- CreateSeuratObject(
            counts = count.small,
            assay = "peak",
            project = lane.id
        )
        Annotation(count.small) <- annotations
        count.small$blacklist_region_fragments <- metadata[paste0(lane.id,"_",colnames(count.small)),"blacklist_region_fragments"]
        count.small$blacklist_ratio <- metadata[paste0(lane.id,"_",colnames(count.small)),"blacklist_ratio"]
        count.small$high.tss  <- metadata[paste0(lane.id,"_",colnames(count.small)),"high.tss"]
        count.small$fragments <- metadata[paste0(lane.id,"_",colnames(count.small)),"fragments"]
        count.small$nucleosome_signal <- metadata[paste0(lane.id,"_",colnames(count.small)),"nucleosome_signal"]
        count.small$TSS.enrichment <- metadata[paste0(lane.id,"_",colnames(count.small)),"TSS.enrichment"]
        count.small$doublet_score <- metadata[paste0(lane.id,"_",colnames(count.small)),"doublet_score"]
        count.small$predicted_doublet <- metadata[paste0(lane.id,"_",colnames(count.small)),"predicted_doublet"]
        count.small$time <- timewindow
        return(count.small)
    })
    stopCluster(cl)
    count[[i]]<- merge(largelist[[1]],
              y = largelist[2:length(largelist)],
              add.cell.ids = lane.list,
              project = timewindow)
    rm(largelist)
}
total <- merge(count[[1]], y = count[2:length(count)],project = timewindows)
rm(count)
total <- FRiP(object = total,
              assay = 'peak',
              total.fragments = 'fragments'
)
total$pct_reads_in_peaks <- total$FRiP*100
total$peak_region_fragments <- total$nCount_peak
total$reads_per_peak <- total$nCount_peak/total$nFeature_peak
saveRDS(total,file=paste0(data_dir,"/total/scATAC_all_time_window_raw_df_with_qc.Rds"))
print("total cell numbers:")
print(length(colnames(total)))

########first round dimensional reduction,clustering and remove clusters########
total <- ScaleData(NormalizeData(total))
total <- RunTFIDF(total, method = 3)
total <- FindTopFeatures(total, min.cutoff = 'q0')
total <- RunSVD(total)
total <- RunPCA(total,npcs = 50, min.dist=0.3)
total <- L2Dim(total, reduction = "pca")
total <- RunUMAP(total, reduction= "pca.l2",dims=2:50)
total <- FindNeighbors(total, reduction = "pca.l2", nn.eps=0.5, dims=2:50)
total <- FindClusters(total,n.start=20, resolution=0.3)
p1 <- DimPlot(total, reduction = "umap", group.by = 'orig.ident',raster=FALSE)
p11 <- p1+theme(legend.position = 'none')
p2 <- DimPlot(total, reduction = "umap",raster=FALSE)
p21 <- p2+theme(legend.position = 'none')
p3 <- DimPlot(total, reduction = "umap",group.by = 'time', raster=FALSE)
p31 <-p3+theme(legend.position = 'none')
clusters <- unique(total$seurat_clusters)
discard_cluster <- ''
total_cell_num <- length(colnames(total))
for (i in 1:length(clusters)) {
  total2 <- subset(total, seurat_clusters == clusters[i])
  x <- total2$orig.ident
  x <- as.data.frame(x)
  a <- summary(factor(x$x))
  b <- sum(a)
  lane_faction <- a/b
  if(max(lane_faction) >= 0.8 | length(colnames(total2)) <100){
    discard_cluster = c(discard_cluster, as.character(clusters[i]))
  }
}
kept_clusters <-as.character(clusters[which(!(clusters %in% discard_cluster))])
print("1st round kept clusters:")
print(sort(kept_clusters))
total2 <- subset(total,seurat_clusters %in% kept_clusters)
print("1st round kept cell numbers:")
print(length(colnames(total2)))
total$kept_clusters <- total$seurat_clusters
total$kept_clusters[total$kept_clusters %in% discard_cluster] <- NA
p4 <- DimPlot(total2, reduction = "umap",raster=FALSE)
p41 <- p4+theme(legend.position = 'none')
p5 <- DimPlot(total, reduction = "umap",group.by = "kept_clusters",raster=FALSE)
p51 <- p5+theme(legend.position = 'none')

pdf(paste0(data_dir,"/total/whole_embeding_umap_1st.pdf"))
print(p11)
print(p21)
print(p31)
print(p41)
print(p51)
dev.off()
umapname <- "whole_embeding_umap_1st_round.png"
ggsave(grid.arrange(p1,p2,p3,p4,p5, ncol=2,nrow=3), filename = umapname, path = paste0(data_dir,'/total'),width = 22,height = 36,units = "in",dpi = 300)

total2@assays$peak@scale.data <- as.matrix(0)
saveRDS(total2, file=paste0(data_dir, '/total/whole_embeding_1st_round.Rds'))

######################call peak by clusters and re-merging######################
peaks <- CallPeaks(
  object = total2,
  group.by = "seurat_clusters",
  macs2.path = "/pathtomacs2",
  broad = FALSE,
  extsize = 150,
  shift = 75,
  effective.genome.size = 1.2e8,
  outdir = paste0(data_dir,"/total"),
  fragment.tempdir = paste0(data_dir,"/total")
)
peaks.sub <- subset(peaks,seqnames %in% c('2R','2L','3L','3R','4'))
count <- list()
for(i in 1:length(timewindows)){
    timewindow_dir <- timewindow_dirs[i]
    timewindow <- timewindows[i]
    lanes <- list.dirs(path = timewindow_dir, recursive = F)
    lane.list <- str_remove(lanes,pattern = paste0(timewindow_dir,"/"))
    lane.list <- sapply(strsplit(lane.list,"_"),"[[",1)
    exclude <- c("analysis","merge-peaks-idr")
    excludedir <- sapply(exclude, function(exclude){paste0(timewindow_dir,"/",exclude)})
    lane.list <- lane.list[which(!(lane.list %in% exclude))]
    lanes <- lanes[which(!(lanes %in% excludedir))]
    clnum<-20
    cl <- makeCluster(getOption("cl.cores", clnum),type="FORK")
    clusterExport(cl, varlist=c("timewindow_dir","annotations","data_dir","peaks.sub"),envir=environment())
    largelist <- parLapply(cl,lanes,function(lanes){
        lane.id <- strsplit(sub(paste0(timewindow_dir,"/"),"" ,lanes),"_")[[1]][1]
        barcode.path <- paste0(timewindow_dir,"/analysis_new/03.peak_by_clusters_and_reclustering/",lane.id,"_kept_cells.tsv")
        barcode.names <- read.delim(barcode.path,header = F, stringsAsFactors = F)
        fragment.name <- list.files(path = lanes, pattern = "fragments.tsv.gz$")
        fragment.path <- paste(lanes,fragment.name, sep = "/")
        kept_cell <- barcode.names$V1
        kept_cell <- sub(paste0(lane.id,"_"),"",kept_cell)
        fragment.object <- CreateFragmentObject(fragment.path, cells = kept_cell,validate.fragments = T)
        peak_matrix <- FeatureMatrix(
            fragments = fragment.object,
            features = peaks.sub
        )
        count.small <- CreateChromatinAssay(
            counts = peak_matrix,
            fragments = fragment.path
        )
        metadata <- read.delim(paste0(timewindow_dir,"/analysis_new/03.peak_by_clusters_and_reclustering/",lane.id,"_metadata.tsv"), header = TRUE, row.names = 1, stringsAsFactors = F)
        metadata <- as.data.frame(metadata)
        count.small <- CreateSeuratObject(
            counts = count.small,
            assay = "peak",
            project = lane.id
        )
        Annotation(count.small) <- annotations
        count.small$blacklist_region_fragments <- metadata[paste0(lane.id,"_",colnames(count.small)),"blacklist_region_fragments"]
        count.small$blacklist_ratio <- metadata[paste0(lane.id,"_",colnames(count.small)),"blacklist_ratio"]
        count.small$high.tss  <- metadata[paste0(lane.id,"_",colnames(count.small)),"high.tss"]
        count.small$fragments <- metadata[paste0(lane.id,"_",colnames(count.small)),"fragments"]
        count.small$nucleosome_signal <- metadata[paste0(lane.id,"_",colnames(count.small)),"nucleosome_signal"]
        count.small$TSS.enrichment <- metadata[paste0(lane.id,"_",colnames(count.small)),"TSS.enrichment"]
        count.small$doublet_score <- metadata[paste0(lane.id,"_",colnames(count.small)),"doublet_score"]
        count.small$predicted_doublet <- metadata[paste0(lane.id,"_",colnames(count.small)),"predicted_doublet"]
        count.small$time <- timewindow
        return(count.small)
    })
    stopCluster(cl)
    count[[i]]<- merge(largelist[[1]],
              y = largelist[2:length(largelist)],
              add.cell.ids = lane.list,
              project = timewindow)
    rm(largelist)
}
total <- merge(count[[1]], y = count[2:length(count)],project = timewindows)
rm(count)
total <- FRiP(object = total,
              assay = 'peak',
              total.fragments = 'fragments'
)
total$pct_reads_in_peaks <- total$FRiP*100
total$peak_region_fragments <- total$nCount_peak
total$reads_per_peak <- total$nCount_peak/total$nFeature_peak
print("2nd total cell numbers:")
print(length(colnames(total)))

#######second round dimensional reduction,clustering and remove clusters########
total <- ScaleData(NormalizeData(total))
total <- RunTFIDF(total, method = 3)
total <- RunSVD(total)
total <- RunPCA(total,npcs = 50, min.dist=0.3)
total <- L2Dim(total, reduction = "pca")
total <- RunUMAP(total, reduction= "pca.l2",dims=2:50)
total <- FindNeighbors(total, reduction = "pca.l2", nn.eps=0.5, dims=2:50)
total <- FindClusters(total,n.start=20, resolution=0.3)
p1 <- DimPlot(total, reduction = "umap", group.by = 'orig.ident',raster=FALSE)
p11 <- p1+theme(legend.position = 'none')
p2 <- DimPlot(total, reduction = "umap",raster=FALSE)
p21 <- p2+theme(legend.position = 'none')
p3 <- DimPlot(total, reduction = "umap",group.by = 'time', raster=FALSE)
p31 <-p3+theme(legend.position = 'none')
clusters <- unique(total$seurat_clusters)
discard_cluster <- ''
total_cell_num <- length(colnames(total))
for (i in 1:length(clusters)) {
  total2 <- subset(total, seurat_clusters == clusters[i])
  x <- total2$orig.ident
  x <- as.data.frame(x)
  a <- summary(factor(x$x))
  b <- sum(a)
  lane_faction <- a/b
  if(max(lane_faction) >= 0.8 | length(colnames(total2)) <100){
    discard_cluster = c(discard_cluster, as.character(clusters[i]))
  }
}
kept_clusters <-as.character(clusters[which(!(clusters %in% discard_cluster))])
print("2nd round kept clusters:")
print(sort(kept_clusters))
total2 <- subset(total,seurat_clusters %in% kept_clusters)
print("2nd round kept cell numbers:")
print(length(colnames(total2)))
total$kept_clusters <- total$seurat_clusters
total$kept_clusters[total$kept_clusters %in% discard_cluster] <- NA
p4 <- DimPlot(total2, reduction = "umap",raster=FALSE)
p41 <- p4+theme(legend.position = 'none')
p5 <- DimPlot(total, reduction = "umap",group.by = "kept_clusters",raster=FALSE)
p51 <- p5+theme(legend.position = 'none')

pdf(paste0(data_dir,"/total/whole_embeding_umap_final.pdf"))
print(p11)
print(p21)
print(p31)
print(p41)
print(p51)
dev.off()
umapname <- "whole_embeding_umap_final.png"
ggsave(grid.arrange(p1,p2,p3,p4,p5, ncol=2,nrow=3), filename = umapname, path = paste0(data_dir,'/total'),width = 22,height = 36,units = "in",dpi = 300)
total2@assays$peak@scale.data <- as.matrix(0)
saveRDS(total2, file=paste0(data_dir, '/total/whole_embeding_final.Rds'))
