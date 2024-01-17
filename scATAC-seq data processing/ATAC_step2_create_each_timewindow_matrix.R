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
source_python("/workdir/script/scrublet_helper.py")
library(methods)
library(cowplot)
library(stringr)
library(ggrastr)
library(data.table)

##########################necessary information#################################
info <- commandArgs(trailingOnly = TRUE)
timewindow <- info[1]###time slot to analyze
cutoff <- info[2]###reads per peak cutoff
cutoff <- as.numeric(cutoff)
data_dir <- info[3]###data dir 
timewindow_dir <- paste(data_dir,timewindow, sep = "/")
lanes <- list.dirs(path = timewindow_dir, recursive = F)
lane.list <- str_remove(lanes,pattern = paste0(timewindow_dir,"/"))
lane.list <- sapply(strsplit(lane.list,"_"),"[[",1)
exclude <- c("merge-peaks-idr","analysis")
excludedir <- sapply(exclude, function(exclude){paste0(timewindow_dir,"/",exclude)})
lane.list <- lane.list[which(!(lane.list %in% exclude))]
lanes <- lanes[which(!(lanes %in% excludedir))]
bdgp_df_subset = bdgp6[bdgp6$chr %in% c('2R','2L','3L','3R','4'),]
bdgp_sub <-bdgp_df_subset %>% mutate(strand=case_when(strand==1~"+",strand==-1~"-"))
colnames(bdgp_sub)[colnames(bdgp_sub) %in% c("ensgene","symbol","biotype")]=c("gene_id","gene_name","gene_biotype")
annotations <- makeGRangesFromDataFrame(bdgp_sub,strand.field = "strand",keep.extra.columns = T) 
seqlevelsStyle(annotations) <- 'NCBI'
genome(annotations) <- "dm6"
dir.create(paste0(timewindow_dir,"/analysis"))
dir.create(paste0(timewindow_dir,"/analysis/01.qc_filtering"))
dir.create(paste0(timewindow_dir,"/analysis/02.clustering"))
dir.create(paste0(timewindow_dir,"/analysis/03.peak_by_clusters_and_reclustering"))
dir.create(paste0(timewindow_dir,"/analysis/03.peak_by_clusters_and_reclustering/peaks"))

####################first round generation of chromatin assay###################  
clnum<-10
cl <- makeCluster(getOption("cl.cores", clnum),type="FORK")
clusterExport(cl, varlist=c("timewindow_dir","annotations"),envir=environment())
largelist <- parLapply(cl,lanes,function(lanes){
  barcode.path <-paste0(lanes,"/barcodes.tsv")
  peak.path <- paste0(timewindow_dir,"/merge-peaks-idr/total.merged.narrowPeak")
  peak.file <- read.delim(peak.path, sep = "\t" ,header = F, stringsAsFactors = F)
  peak.file <- peak.file[peak.file$V1 %in% c('2R','2L','3L','3R','4'),]
  colnames(peak.file) <- c("seqnames","start","end")
  peak_ranges <- makeGRangesFromDataFrame(peak.file)
  barcode.names <- read.delim(barcode.path,header = F, stringsAsFactors = F)
  fragment.name <- list.files(path = lanes, pattern = "fragments.tsv.gz$") 
  fragment.path <- paste(lanes,fragment.name, sep = "/")
  fragment.object <- CreateFragmentObject(fragment.path, cells = barcode.names$V1,validate.fragments = T)
  peak_matrix <- FeatureMatrix(
    fragments = fragment.object,
    features = peak_ranges
  )
  mincell=dim(peak_matrix)[2]%/%100
  count.small <- CreateChromatinAssay(
    counts = peak_matrix, 
    min.cells = mincell, 
    fragments = fragment.path
  )
  lane.id <- strsplit(sub(paste0(timewindow_dir,"/"),"" ,lanes),"_")[[1]][1]
  count.small <- CreateSeuratObject(counts = count.small,
                                    assay = "peak",
                                    project = lane.id
  )
  Annotation(count.small) <- annotations
  count.small <- TSSEnrichment(object = count.small, fast = T)
  count.small$high.tss <- ifelse(count.small$TSS.enrichment > quantile(count.small$TSS.enrichment)[2], 'High', 'Low')
  fragments_info<- CountFragments(fragments = fragment.path)
  rownames(fragments_info)<-fragments_info[,1]
  count.small$fragments <- fragments_info[colnames(count.small), "reads_count"]
  scrublet_res <- do.call(scrublet_py,
                          list(data=as(Matrix::t(count.small@assays[['peak']]@counts), "TsparseMatrix"),
                               expected_doublet_rate=0.06, min_counts=2, min_cells=3,
                               min_gene_variability_pctl=85, n_prin_comps=50, sim_doublet_ratio=2,
                               n_neighbors=round(0.5*sqrt(ncol(count.small)))))
  count.small$doublet_score <- scrublet_res$scores
  doublet_cutoff = quantile(scrublet_res$scores,0.95)
  scrublet_res$is_doublets <- scrublet_res$scores>doublet_cutoff
  count.small$predicted_doublet <- ifelse(scrublet_res$is_doublets, 'Doublet', 'Singlet')
  new_thresh <- 0.95
  
  a_df <- data.frame(score=scrublet_res$scores, class='Observed') %>%
    rbind(data.frame(score=scrublet_res$sim_scores, class='Simulated'))
  p <- a_df %>%
    ggplot(aes(x=score)) + geom_histogram() +
    facet_wrap(class~., ncol=2, scales='free') +
    ggtitle(paste0('Scrublet threshold: ', round(new_thresh, 2))) +
    geom_vline(xintercept=new_thresh, linetype='dashed') +
    xlab("Scrublet-based doublet score") +
    ylab("Count") +
    scale_y_log10()
  doublet_his <- paste0(timewindow,"_",lane.id,"_doublet_histogram_idr.png")
  ggsave(p, filename = doublet_his, path = paste0(timewindow_dir,"/analysis/01.qc_filtering"), height=4, width=8)
  return(count.small)
})
stopCluster(cl)
count<- merge(largelist[[1]], 
              y = largelist[2:length(largelist)], 
              add.cell.ids = lane.list, 
              project = timewindow)
rm(largelist)

############################first round qc checking#############################
count <- NucleosomeSignal(object = count)
blacklist <- blacklist_dm6
seqlevelsStyle(blacklist) <- 'NCBI'
count$blacklist_region_fragments <-CountsInRegion(object = count,
                                                  assay = 'peak',
                                                  regions = blacklist
)
count$blacklist_ratio <- FractionCountsInRegion(object = count,
                                                assay = 'peak',
                                                regions = blacklist
)
count <- FRiP(object = count,
              assay = 'peak',
              total.fragments = 'fragments'
)
count$pct_reads_in_peaks <- count$FRiP*100
count$peak_region_fragments <- count$nCount_peak
count$reads_per_peak <- count$nCount_peak/count$nFeature_peak
saveRDS(count,paste0(timewindow_dir, '/analysis/01.qc_filtering/raw_df_with_qc_idr_', timewindow, '.Rds'))
p <- VlnPlot(count,features =  c('pct_reads_in_peaks', 'peak_region_fragments',
                                 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal','reads_per_peak'),
             pt.size = 0.1,
             ncol = 6
)
vlnname <- paste0(timewindow,"_violin_plot.png")
ggsave(p, filename = vlnname,path = paste0(timewindow_dir,'/analysis/01.qc_filtering'),width = length(lane.list)*10,height = 7,units = "in",dpi = 300, limitsize =F)

###########################first round qc filtering#############################
count.sub <- subset(
  x = count,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 40000 &
    pct_reads_in_peaks > 20 &
    blacklist_ratio < 0.05 &
    TSS.enrichment > quantile(count$TSS.enrichment)[2] &
    nucleosome_signal < 4 &
    reads_per_peak >cutoff
)
print("1st round subset cell numbers:")
print(length(colnames(count.sub)))

######first round dimensional reduction, remove batch effect and clustering#####
count.sub <- ScaleData(NormalizeData(count.sub))
count.sub <- RunTFIDF(count.sub, method = 3)
count.sub <- RunSVD(count.sub)
count.sub <- RunPCA(count.sub,npcs = 50, min.dist=0.3)
count.sub <- L2Dim(count.sub, reduction = "pca")
count.sub2 <- RunUMAP(count.sub, reduction= "pca.l2",dims=2:50)
p1 <- DimPlot(count.sub2, reduction = "umap", group.by = 'predicted_doublet',raster=FALSE)
count.sub <- subset(  x = count.sub,
                      subset =  predicted_doublet == 'Singlet')
saveRDS(count.sub, file=paste0(timewindow_dir, '/analysis/02.clustering/dedoublet_', timewindow, '.Rds'))

count.sub2 <- RunUMAP(count.sub, reduction= "pca.l2",dims=2:50)
p2 <- DimPlot(count.sub2, reduction = "umap",raster=FALSE)
count.sub <- RunHarmony(count.sub, reduction="pca.l2", group.by.vars="orig.ident")
count.sub <- RunUMAP(count.sub, reduction= "harmony",dims=2:50)
p3 <- DimPlot(count.sub, reduction = "umap",raster=FALSE)
count.sub <- FindNeighbors(count.sub, reduction = "harmony", nn.eps=0.5, dims=2:50)
count.sub <- FindClusters(count.sub,n.start=20, resolution=0.3)
p4 <- DimPlot(count.sub, reduction = "umap",raster=FALSE)

##########################first round filter clusters###########################
clusters <- unique(count.sub$seurat_clusters)
discard_cluster <- ''
total_cell_num <- length(colnames(count.sub))
for (i in 1:length(clusters)) {
  count.sub2 <- subset(count.sub, seurat_clusters == clusters[i])
  x <- count.sub2$orig.ident
  x <- as.data.frame(x)
  a <- summary(factor(x$x))
  b <- sum(a)
  lane_faction <- a/b
  if(max(lane_faction) >= 0.8 | length(colnames(count.sub2)) < total_cell_num*0.005 | length(colnames(count.sub2)) <100){
    discard_cluster = c(discard_cluster, as.character(clusters[i]))
  }
}
kept_clusters <-as.character(clusters[which(!(clusters %in% discard_cluster))])
print("1st round kept clusters:")
print(sort(kept_clusters))
count.sub2 <- subset(count.sub,seurat_clusters %in% kept_clusters)
print("kept cell numbers:")
print(length(colnames(count.sub2)))
count.sub$kept_clusters <- count.sub$seurat_clusters
count.sub$kept_clusters[count.sub$kept_clusters %in% discard_cluster] <- NA
p5 <- DimPlot(count.sub2, reduction = "umap",raster=FALSE)
p6 <- DimPlot(count.sub, reduction = "umap",group.by = "kept_clusters",raster=FALSE)
umapname <- paste0(timewindow, "_dedoublet_batch_harmony_cluster_filtered_umap.png")
ggsave(grid.arrange(p1,p2,p3,p4,p5,p6, ncol=2,nrow=3), filename = umapname, path = paste0(timewindow_dir,'/analysis/02.clustering/'),width = 22,height = 24,units = "in",dpi = 300)
saveRDS(count.sub2, file=paste0(timewindow_dir, '/analysis/02.clustering/clustered_de_batch_effect_filtered_df_', timewindow, '.Rds'))

#####################call peak by clusters and re-clustering####################
count <- count.sub2
peaks <- CallPeaks(
  object = count,
  group.by = "seurat_clusters",
  macs2.path = "/pathtomacs2",
  broad = FALSE,
  extsize = 150,
  shift = 75,
  effective.genome.size = 1.2e8,
  outdir = paste0(timewindow_dir,"/analysis/03.peak_by_clusters_and_reclustering/peaks"),
  fragment.tempdir = paste0(timewindow_dir,"/analysis/03.peak_by_clusters_and_reclustering/peaks")
)
peaks.sub <- subset(peaks,seqnames %in% c('2R','2L','3L','3R','4'))
clnum<-10
cl <- makeCluster(getOption("cl.cores", clnum),type="FORK")
clusterExport(cl, varlist=c("timewindow_dir","annotations","peaks.sub"),
              envir=environment())
largelist <- parLapply(cl,lanes, function(lanes){
  barcode.path <-paste0(lanes,"/barcodes.tsv")
  barcode.names <- read.delim(barcode.path,header = F, stringsAsFactors = F)
  fragment.name <- list.files(path = lanes, pattern = "fragments.tsv.gz$") 
  fragment.path <- paste(lanes,fragment.name, sep = "/")
  fragment.object <- CreateFragmentObject(fragment.path, cells = barcode.names$V1,validate.fragments = T)
  peak_matrix <- FeatureMatrix(
    fragments = fragment.object,
    features = peaks.sub
  )
  mincell=dim(peak_matrix)[2]%/%100
  count.small <- CreateChromatinAssay(
    counts = peak_matrix, 
    min.cells = mincell, 
    fragments = fragment.object
  )
  lane.id <- strsplit(sub(paste0(timewindow_dir,"/"),"" ,lanes),"_")[[1]][1]
  count.small <- CreateSeuratObject(counts = count.small,
                                    assay = "peak",
                                    project = lane.id
  )
  fragments_info<- CountFragments(fragments = fragment.path)
  rownames(fragments_info)<-fragments_info[,1]
  count.small$fragments <- fragments_info[colnames(count.small), "reads_count"]
  Annotation(count.small) <- annotations
  count.small <- TSSEnrichment(object = count.small, fast = T)
  count.small$high.tss <- ifelse(count.small$TSS.enrichment > quantile(count.small$TSS.enrichment)[2], 'High', 'Low') 
  scrublet_res <- do.call(scrublet_py,
                          list(data=as(Matrix::t(count.small@assays[['peak']]@counts), "TsparseMatrix"),
                               expected_doublet_rate=0.06, min_counts=2, min_cells=3,
                               min_gene_variability_pctl=85, n_prin_comps=50, sim_doublet_ratio=2,
                               n_neighbors=round(0.5*sqrt(ncol(count.small)))))
  count.small$doublet_score <- scrublet_res$scores
  doublet_cutoff = quantile(scrublet_res$scores,0.95)
  scrublet_res$is_doublets <- scrublet_res$scores>doublet_cutoff
  count.small$predicted_doublet <- ifelse(scrublet_res$is_doublets, 'Doublet', 'Singlet')
  return(count.small)
})
stopCluster(cl)
count<- merge(largelist[[1]],
              y = largelist[2:length(largelist)],
              add.cell.ids = lane.list,
              project = timewindow)
rm(largelist)

##########################second round qc checking##############################
count <- NucleosomeSignal(object = count)
blacklist <- blacklist_dm6
seqlevelsStyle(blacklist) <- 'NCBI'
count$blacklist_region_fragments <-CountsInRegion(object = count,
                                                  assay = 'peak',
                                                  regions = blacklist
)
count$blacklist_ratio <- FractionCountsInRegion(object = count,
                                                assay = 'peak',
                                                regions = blacklist
)
count <- FRiP(object = count,
              assay = 'peak',
              total.fragments = 'fragments'
)
count$pct_reads_in_peaks <- count$FRiP*100
count$peak_region_fragments <- count$nCount_peak
count$reads_per_peak <- count$nCount_peak/count$nFeature_peak
p <- VlnPlot(
  object = count,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal','reads_per_peak'),
  pt.size = 0.1,
  ncol = 6
)
vlnname <- paste0(timewindow,"_recall_peak_violin_plot_filtered.png")
ggsave(p, filename = vlnname,path = paste0(timewindow_dir,'/analysis/03.peak_by_clusters_and_reclustering'),width = length(lane.list)*10,height = 7,units = "in",dpi = 300, limitsize =F)
saveRDS(count, file=paste0(timewindow_dir, '/analysis/03.peak_by_clusters_and_reclustering/recall_peak_raw_df_with_qc_filtered_', timewindow, '.Rds'))

#######################second round qc filtering################################
count.sub <- subset(
  x = count,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 40000 &
    pct_reads_in_peaks > 20 &
    blacklist_ratio < 0.05 &
    TSS.enrichment > quantile(count$TSS.enrichment)[2] &
    nucleosome_signal < 4&
    reads_per_peak > cutoff
)
print("2nd round subset cell numbers:")
print(length(colnames(count.sub)))

#####second round dimensional reduction, remove batch effect and clustering#####
count.sub <- ScaleData(NormalizeData(count.sub))
count.sub <- RunTFIDF(count.sub, method = 3)
count.sub <- RunSVD(count.sub)
count.sub <- RunPCA(count.sub,npcs = 50, min.dist=0.3)
count.sub <- L2Dim(count.sub, reduction = "pca")
count.sub2 <- RunUMAP(count.sub, reduction= "pca.l2",dims=2:50)
p1 <- DimPlot(count.sub2, reduction = "umap", group.by = 'predicted_doublet',raster=FALSE)
count.sub <- subset(  x = count.sub,
                      subset =  predicted_doublet == 'Singlet')
count.sub2 <- RunUMAP(count.sub, reduction= "pca.l2",dims=2:50)
p2 <- DimPlot(count.sub2, reduction = "umap",raster=FALSE)
count.sub <- RunHarmony(count.sub, reduction="pca.l2", group.by.vars="orig.ident")
count.sub <- RunUMAP(count.sub, reduction= "harmony",dims=2:50)
p3 <- DimPlot(count.sub, reduction = "umap",raster=FALSE)
count.sub <- FindNeighbors(count.sub, reduction = "harmony", nn.eps=0.5, dims=2:50)
count.sub <- FindClusters(count.sub,n.start=20, resolution=0.3)
p4 <- DimPlot(count.sub, reduction = "umap",raster=FALSE)

#########################second round filter clusters###########################
clusters <- unique(count.sub$seurat_clusters)
discard_cluster <- ''
total_cell_num <- length(colnames(count.sub))
for (i in 1:length(clusters)) {
  count.sub2 <- subset(count.sub, seurat_clusters == clusters[i])
  x <- count.sub2$orig.ident
  x <- as.data.frame(x)
  a <- summary(factor(x$x))
  b <- sum(a)
  lane_faction <- a/b
  if(max(lane_faction) >= 0.8 | length(colnames(count.sub2)) < total_cell_num*0.005 | length(colnames(count.sub2)) <100){
    discard_cluster = c(discard_cluster, as.character(clusters[i]))
  }
}
kept_clusters <-as.character(clusters[which(!(clusters %in% discard_cluster))])
print("2nd round kept clusters:")
print(sort(kept_clusters))
count.sub2 <- subset(count.sub,seurat_clusters %in% kept_clusters)
print("kept cell numbers:")
print(length(colnames(count.sub2)))
count.sub$kept_clusters <- count.sub$seurat_clusters
count.sub$kept_clusters[count.sub$kept_clusters %in% discard_cluster] <- NA
p5 <- DimPlot(count.sub2, reduction = "umap",raster=FALSE)
p6 <- DimPlot(count.sub, reduction = "umap",group.by = "kept_clusters",raster=FALSE)
umapname <- paste0(timewindow, "_recall_peak_batch_harmony_cluster_filtered_umap.png")
ggsave(grid.arrange(p1,p2,p3,p4,p5,p6, ncol=2,nrow=3), filename = umapname, path = paste0(timewindow_dir,'/analysis/03.peak_by_clusters_and_reclustering'),width = 22,height = 24,units = "in",dpi = 300)
saveRDS(count.sub2, file=paste0(timewindow_dir, '/analysis/03.peak_by_clusters_and_reclustering/Reclustered_de_batch_effect_filtered_df_', timewindow, '.Rds'))

#################save necessary information for whole embedding#################
peak <- rownames(count.sub2)
peak <- str_split(peak,"-",simplify = T)
write.table(peaks,file = paste0(timewindow_dir, '/analysis/03.peak_by_clusters_and_reclustering/kept_peaks.bed'),quote = F,sep = "\t",col.names = F,row.names = F)
lanes <- unique(mat$orig.ident)
lapply(lanes,function(lanes){
  sub <- subset(count.sub2, orig.ident == lanes)
  barcode <- colnames(sub)
  write.table(barcode,file = paste0(timewindow_dir, '/analysis/03.peak_by_clusters_and_reclustering/',lanes,'_kept_cells.tsv'),quote = F,sep = "\t",col.names = F,row.names = F)
  meta <- sub@meta.data[,1:15]
  write.table(meta,file = paste0(timewindow_dir, '/analysis/03.peak_by_clusters_and_reclustering/',lanes,'_metadata.tsv'),quote = F,sep = "\t",col.names = T,row.names = T)
})