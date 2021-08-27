---
title: "Modified ScRNA script for Wang et al 2020 Nat Commun"
author: "Foong Lian Chee"
date: "8/25/2021"
output: html_document
---

I reproduced the ScRNA analysis for a Nat Comm paper: https://www.nature.com/articles/s41467-020-14296-y#Sec2


Wang, X., Xu, H., Cheng, C., Ji, Z., Zhao, H., Sheng, Y., ... & Zhu, H. H. (2020). Identification of a Zeb1 expressing basal stem cell subpopulation in the prostate. Nature communications, 11(1), 1-16.


Here's my modified and updated script to reproduce the works. 


```{css, echo=FALSE}
.scroll-100 {
  max-height: 100px;
  overflow-y: auto;
  background-color: inherit;
}
```
```{r, class.output="scroll-100", message=FALSE, fig.width= 10, fig.height=6 }
#setwd("R_HelenHeZhu_StemCell")

###############################################################################
### Step01 removal of poor quality cells and contaminated non-epithelial cells 
###############################################################################

library(Seurat)
library(dplyr)
library(Matrix)
library(ggpubr)
#

data <- Read10X(data.dir = "GSE111429_RAW")
PG_all    <- CreateSeuratObject(count = data, min.cells = 0, min.genes = 0, project = "PG")
PG_all <- PercentageFeatureSet(PG_all, pattern = "^mt-", col.name = "percent.mt")
before <- VlnPlot(PG_all, features =c('nFeature_RNA','nCount_RNA','percent.mt'))

# Poor-quality cells with less than 1000 genes detected, less than 5000 UMIs or more than 5% UMI mapped to mitochondria genes were removed. 
PG_all <- subset(PG_all, subset = nFeature_RNA > 1400 & nCount_RNA > 5000 & percent.mt < 5)
after <- VlnPlot(PG_all, features =c('nFeature_RNA','nCount_RNA','percent.mt'))
ggarrange(before,after, ncol = 2, nrow = 1)
```
```{r, class.output="scroll-100", message=FALSE, out.width="50%"}
#
PG_all <- NormalizeData(PG_all, normalization.method='LogNormalize', scale.factor=1e4)
PG_all <- FindVariableFeatures(PG_all,
                               nfeatures = 1606,
                               mean.cutoff = c(0.0125,3),
                               dispersion.cutoff = c(0.5, Inf))
length(PG_all@assays$RNA@var.features)
summary(PG_all@assays$RNA@meta.features)
#
PG_all <- ScaleData(PG_all, vars.to.regress=c('nCount_RNA','percent.mt'))
PG_all <- RunPCA(PG_all, features = PG_all@assays$RNA@var.features, verbose = TRUE, 
                 ndims.print=1:5, nfeatures.print = 5)
PCAPlot(PG_all, dims = c(1, 2))
# 
DimHeatmap(PG_all, dims = 1:12, cells = 500, balanced = TRUE)
#
PG_all <- FindNeighbors(PG_all, dims = 1:12)
PG_all <- FindClusters(PG_all, resolution=c(0.2,0.3,0.4,0.5))
PG_all <- RunTSNE(PG_all, dims = 1:12, do.fast=TRUE)
Idents(object = PG_all) <- PG_all$RNA_snn_res.0.2
DimPlot(PG_all, reduction = "tsne", label = TRUE, pt.size = 1, label.size = 6)
PCAPlot(PG_all)
####
```
```{r, class.output="scroll-100", message=FALSE, out.width="50%"}
####
#### retrieve non-epithelial cells
FeaturePlot(PG_all, features =c('Cd74','Cd72','Cd54') , cols =c('yellow','red'), 
            reduction ='tsne', pt.size = 0.7,ncol = 3, label = T)
FeaturePlot(PG_all, features = c('Eng','S1pr1','Emcn'), cols =c('yellow','red'), 
            reduction = 'tsne',pt.size = 0.7,ncol=3, label = T)
#
#Immu.Cluster <- FindMarkers(PG_all, ident.1='9', only.pos=TRUE, min.pct=0.25)
#Endo.Cluster <- FindMarkers(PG_all, ident.1='4', only.pos=TRUE, min.pct=0.25)
#Immu.Cluster$Symbol <- rownames(Immu.Cluster)
#Endo.Cluster$Symbol <- rownames(Endo.Cluster)
```
```{r, fig.height=14, class.output="scroll-100", message=FALSE}
# 
v1 <- VlnPlot(PG_all, features = c("Epcam","Cdh1","Krt5", "Krt14"), ncol = 4) #Epithelial markers
v2 <- VlnPlot(PG_all, features = c("Igfbp4","Fn1","Gng11"), ncol = 3) # stromal markers
v3 <- VlnPlot(PG_all, features = c("Eng","S1pr1","Emcn"), ncol = 3) # endothelial markers
v4 <- VlnPlot(PG_all, features = c("Cd74","Cd72"), ncol = 2) #immune markers
ggarrange(v1, v2, v3, v4, ncol = 1)
#
```
```{r, out.width="50%", class.output="scroll-100", message=FALSE}
###
### Draw heatmap
### 
library(ComplexHeatmap)
library(circlize)
#
geneSets                 <- c('Cd74','Cd72','Cd54','Eng','S1pr1','Emcn')
cellRanks                <- PG_all@meta.data[order(PG_all@meta.data$RNA_snn_res.0.2),]
PartialMatrix            <- PG_all@assays$RNA@scale.data[which(rownames(PG_all@assays$RNA@scale.data) %in% geneSets), rownames(cellRanks)]
cellRanks$col            <- rainbow(max(as.numeric(cellRanks$RNA_snn_res.0.2))+1)[as.numeric(cellRanks$RNA_snn_res.0.5)+1]
#
ha_column <- HeatmapAnnotation(
  df  = data.frame(
    ClusterID = as.numeric(cellRanks$RNA_snn_res.0.2)
  ),
  col = list(
    ClusterID = colorRamp2(unique(as.numeric(cellRanks$RNA_snn_res.0.2)), 
                           rainbow(max(as.numeric(cellRanks$RNA_snn_res.0.2))))
  )
)
ht1 = Heatmap(PartialMatrix, name = "Scaled\nUMI", column_title = "non-epithelial cells markers", 
              top_annotation = ha_column, show_column_names=FALSE, cluster_columns = FALSE,
              row_names_gp = gpar(fontsize = 10))
draw(ht1)
```
```{r, out.width="50%", class.output="scroll-100", message=FALSE}
# seurat function
### heatmap for non-epithelial cells markers
geneSets                 <- c('Cd74','Cd72','Cd54','Eng','S1pr1','Emcn')
DoHeatmap(PG_all, features = geneSets)

```
```{r, out.width="50%", class.output="scroll-100", message=FALSE}
################################################################################
### Step02 clustering using Seurat #############################################
################################################################################
#
PG_all <- subset(PG_all, idents = c("4","8"), invert = TRUE)
DefaultAssay(PG_all) <- "RNA"
PG_all <- NormalizeData(PG_all, normalization.method='LogNormalize', scale.factor=1e4)
PG_all <- FindVariableFeatures(PG_all, mean.function = ExpMean, 
                               dispersion.function= LogVMR,
                               mean.cutoff = c(0.0125, 3),
                               dispersion.cutoff = c(0.5, Inf))
length(PG_all@assays$RNA@var.features)
#
PG_all <- ScaleData(PG_all, vars.to.regress=c('nCount_RNA','percent.mt'))
PG_all <- RunPCA(PG_all, features = PG_all@assays$RNA@var.features, verbose = TRUE, 
                 ndims.print=1:5, nfeatures.print = 5)
PG_all <- FindNeighbors(PG_all, dims = 1:12)
PG_all <- FindClusters(PG_all, resolution=c(0.2,0.3,0.4,0.5))
PG_all <- RunTSNE(PG_all, dims = 1:12, do.fast=TRUE)
#saveRDS(PG_all, file = "PG_all.RDS")
```
```{r, out.width="50%", class.output="scroll-100", message=FALSE}
################################################
PG_all <- readRDS("PG_all.RDS")
Idents(object = PG_all) <- PG_all$RNA_snn_res.0.5
DimPlot(PG_all, label = T, pt.size = 1, label.size = 6)

plot <- VlnPlot(PG_all, features = c("Krt5","Krt14"), combine = FALSE, 
                fill.by = c("feature","ident")) 
plot <-  lapply(
  X = plot,
  FUN = function(p) p + aes(color= PG_all$RNA_snn_res.0.5)
)
CombinePlots(plots = plot, legend = 'none')
```
```{r, out.width="50%", class.output="scroll-100"}
################################################
# Doheatmap
listsA <- c('Cdh1','Epcam','Cldn1','Ocln','Vim','Fn1','S100a4','Zeb1','Zeb2',
            'Twist1','Twist2','Snai1','Snai2','Prrx1','Prrx2','Foxc2','Cdh2','Acta2')
# Option 1
DoHeatmap(PG_all, features = listsA)

# Option 2: customise
DoHeatmap(PG_all, features = listsA, disp.min = -1,
          slot = 'scale.data', 
          size = 3,
          group.colors = rainbow(9)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",midpoint = 1)
#cluster 7 showed both epithelial and stromal markers

```
```{r, out.width="50%", class.output="scroll-100"}
################################################################################
# Step 03 DEGs and GSEA analysis
################################################################################
dev.off()
PG.markers.7 <- FindMarkers(PG_all, ident.1 = "7", min.pct = 0.20)
PG.markers.7$gene <- rownames(PG.markers.7)
library(biomaRt)
gene_id <- getBM(attributes = c("ensembl_gene_id","external_gene_name","entrezgene_id"),
                 mart = useDataset("mmusculus_gene_ensembl", useMart("ensembl")))
PG.markers.7 <- merge(PG.markers.7, gene_id[,c(2,3)], by.x = "gene", by.y = "external_gene_name")
PG.markers.7.filt <- PG.markers.7
PG.markers.7.filt <- PG.markers.7.filt[!duplicated(PG.markers.7.filt$entrezgene_id),]
#PG.markers.7.filt <- PG.markers.7[which(PG.markers.7$p_val_adj < 0.05),]
genelist <- PG.markers.7$avg_log2FC
names(genelist) <-  PG.markers.7$entrezgene_id
genelist <- sort(genelist, decreasing = TRUE)

library(fgsea)
library(msigdbr)
library(data.table)
library(ggplot2)
m_df <- msigdbr(species = "Mus musculus", category = "C2") %>% dplyr::select(gs_name, entrez_gene)
m_df_fgsea <- split(x = m_df$entrez_gene, f = m_df$gs_name)
Zhang_BasHi <- read.csv("Zhang_BasHi.csv", header = FALSE, fileEncoding="UTF-8-BOM") # add other cell markers from the publication
#Zhang_LumHi <- read.csv("Zhang_LumHi.csv") # add other cell markers from the publication
Zhang_BasHi <- data.frame(Zhang_BasHi)
colnames(Zhang_BasHi)[1] <- "gene"
library(tidyverse)
Zhang_BasHi$gene <- str_to_sentence(Zhang_BasHi$gene)
colnames(Zhang_BasHi)[1] <- "gene"
Zhang_BasHi <- merge(Zhang_BasHi, gene_id[,c(2,3)],by.x = "gene", by.y = "external_gene_name")
Zhang_BasHi <- Zhang_BasHi[,2]
m_df_fgsea[["Zhang_BasHi"]] <- Zhang_BasHi
  
fgseaRes <- fgsea(pathways = m_df_fgsea, 
                  stats    = genelist,
                  eps      = 0.05,
                  minSize  = 15,
                  maxSize  = 800)
selected <- fgseaRes[c(979, 329, 443, 445, 907),]
selected <- selected$pathway
plotGseaTable(m_df_fgsea[selected], genelist, fgseaRes, 
              gseaParam=0.5)
```
```{r, out.width="50%", class.output="scroll-100"}
#### Pathview ####
dev.off()
library(pathview)
dme <- pathview(gene.data= genelist, pathway.id="mmu04310", species = "mmu")
knitr::include_graphics("mmu04310.pathview.png")
```
```{r, out.width="50%", class.output="scroll-100"}
################################################################################
#### Step 04 Trajectory analysis ########
################################################################################

library(monocle3)
library(dplyr)
library(SeuratWrappers)
library(ggplot2)
data <- GetAssayData(PG_all)[VariableFeatures(PG_all),]
pd <- data.frame(PG_all@meta.data)
fData <- data.frame(gene_short_name = rownames(data), row.names = row.names(data))
cds <- new_cell_data_set(expression_data = data,
                         cell_metadata = pd,
                         gene_metadata = fData)
cds <- preprocess_cds(cds, num_dim = 12)
#plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds,
                        preprocess_method = "PCA",
                        reduction_method= "tSNE",
                        verbose = T)
cds <- cluster_cells(cds = cds , reduction_method = "tSNE", verbose = TRUE,
                     resolution = NULL)
p1 <- plot_cells(cds, reduction_method = "tSNE")
#compare with seurat data, res = 0.5
cds@clusters$tSNE$clusters <- PG_all$RNA_snn_res.0.5
p2 <- plot_cells(cds, reduction_method = "tSNE")
p1 + p2
```
```{r}
# learn graph using UMAP reduction method
cds <- reduce_dimension(cds,
                        reduction_method= "UMAP",
                        verbose = T)
cds <- cluster_cells(cds = cds , reduction_method = "UMAP", verbose = TRUE)
plot_cells(cds, reduction_method = "UMAP", show_trajectory_graph = FALSE, group_label_size = 5)
cds <- learn_graph(cds, use_partition = TRUE)
# save pre-selected cluster 7 cells as the root cells
pre_selected_cells <- colnames(subset(PG_all, idents = "7"))
cds <- order_cells(cds, reduction_method = "UMAP", root_cells = pre_selected_cells)
```
```{r}
# plot trajectories colored by pseudotime
plot_cells(
  cds = cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE,
  group_label_size = 5,
  graph_label_size = 3
)

cds@clusters$UMAP$clusters <- PG_all$RNA_snn_res.0.5
plot_cells(
  cds = cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = F,
  group_label_size = 5,
  graph_label_size = 3
)
plot_cells(
  cds = cds,
  color_cells_by = "cluster",
  show_trajectory_graph = FALSE,
  group_label_size = 5,
  graph_label_size = 3
)

plot_cells(
  cds = cds,
  color_cells_by = "cluster",
  group_label_size = 5,
  graph_label_size = 4,
  label_leaves = F,
  label_branch_points = F,
  show_trajectory_graph = TRUE,
  scale_to_range = T
)
```
```{r, out.width="50%", class.output="scroll-100"}
# Add pseudotime value to the seurat object
PG_all <- AddMetaData(object = PG_all, 
                      metadata = cds@principal_graph_aux@listData$UMAP$pseudotime,
                      col.name =  "pseudotime")
FeaturePlot(PG_all, "pseudotime", reduction = "tsne", label = T, label.size = 8) + scale_color_viridis_c()
F1 <- FeaturePlot(PG_all, "pseudotime", reduction = "pca", label = T, label.size = 6, repel = T) + scale_color_viridis_c()
F2 <- DimPlot(PG_all, reduction = "pca", group.by = "seurat_clusters", label = T, repel = T)
ggarrange(F1, F2, ncol = 2)
```
```{r}
################################################################################
#### Diffusion Map ####
################################################################################
# package 'destiny' is not available for Bioconductor version '3.13'
```
```{r, class.output="scroll-100"}
################################################################################
#### Slingshot ####
################################################################################
library(slingshot)

pd <- data.frame(PG_all@meta.data)
fData <- data.frame(gene_short_name = rownames(data), row.names = row.names(data))
counts <- GetAssayData(PG_all)[VariableFeatures(PG_all),]
sce <- SingleCellExperiment(list(counts = counts),
                            metadata = pd)


# use PCA
reducedDims(sce) <- list(pca = PG_all@reductions$pca@cell.embeddings)
reducedDims(sce)$pca <- reducedDims(sce)$pca[, 1:12]
sce
sce <- slingshot(sce,
                 reducedDim = "pca",
                 clusterLabels = sce@metadata$RNA_snn_res.0.5,
                 start.clus = 7)
sds <- SlingshotDataSet(sce)
sds  # has 5 or 6 lineages

```
```{r, class.output="scroll-100"}
##########################################################################
############## visualization for pca results ##############################
#######################################################################
# Plot clusters with identified lineages overlayed
rd <- data.frame(reducedDim(sce))
cl <- unlist(sce@metadata$RNA_snn_res.0.5)
names(cl) <- rownames(PG_all@meta.data)

par(xpd=TRUE)
par(mar=c(4.5,5.5,2,7))
plot(rd[,1:2], col = rainbow(10)[cl], asp = 0.5, pch = 20)
lines(SlingshotDataSet(sce),  lwd = 2, col = 'black')
legend(x = 30, y = 50, legend= order(unique(cl)), 
       fill=rainbow(10)[order(unique(cl))])
```
```{r, class.output="scroll-100"}
# Add into Seurat object
PG_all$slingPseudotime_1 <- sce$slingPseudotime_1
PG_all$slingPseudotime_2 <- sce$slingPseudotime_2
PG_all$slingPseudotime_3 <- sce$slingPseudotime_3
PG_all$slingPseudotime_4 <- sce$slingPseudotime_4
PG_all$slingPseudotime_5 <- sce$slingPseudotime_5
PG_all$slingPseudotime_6 <- sce$slingPseudotime_6
 
FeaturePlot(PG_all, c("slingPseudotime_1","slingPseudotime_2","slingPseudotime_3",
                      "slingPseudotime_4","slingPseudotime_5","slingPseudotime_6"), 
            label = T, label.size = 5, cols = rainbow(12))

```
```{r, out.width="50%", class.output="scroll-100"}
##########################################################################
### Visualize pseudotime using pca data
##########################################################################
## scatter3d plot, plotly 
# refer to https://github.com/Dragonmasterx87/Interactive-3D-Plotting-in-Seurat-3.0.0/blob/master/3D%20UMAP%20Plotting%20v1.3.R
library(rgl)
library(plotly)


#Embeddings(object = PG_all, reduction = "pca")
#summary(Embeddings(object = PG_all, reduction = "pca"))
plot.data <- FetchData(object = PG_all, vars = c("PC_1", "PC_2", "PC_3", "seurat_clusters", "slingPseudotime_1","slingPseudotime_2","slingPseudotime_3","slingPseudotime_4","slingPseudotime_5","slingPseudotime_6"))
plot.data$label <- paste(rownames(plot.data))
fig <- plot_ly(data = plot.data, 
               x = ~PC_1, y = ~PC_2, z = ~PC_3, 
               color = ~seurat_clusters, 
               colors = rainbow(10),
               type = "scatter3d", 
               mode = "markers", 
               marker = list(size = 2, width=2), # controls size of points
               text=~seurat_clusters, #This is that extra column we made earlier for which we will use for cell ID
               hoverinfo="text") #When you visualize your plotly object, hovering your mouse pointer over a point shows cell names
               
fig <- fig %>% layout(scene = list(xaxis= list(title = "PC_1"),
                                   yaxis= list(title = "PC_2"),
                                   zaxis= list(title = "PC_3")))
fig 
```
```{r}
## Failed to add lineage
# <- fig %>% add_trace(type = 'scatter3d', data = plot.data,
               #          x = ~slingPseudotime_1, y = ~slingPseudotime_2, y = ~slingPseudotime_3,
                #         line = list(color = 'black', width = 1))

#fig
```
```{r, out.width="50%", class.output="scroll-100"}
##########################################################################
### Helen's script for the scatter3d plot
###########################################################################
Dims.pca          <- sds@reducedDim[,1:3]         # pca dimension
clusterLabels.pca <- sds@clusterLabels        # cluster labels
connectivity.pca  <- sds@adjacency         # 
clusters.pca      <- rownames(connectivity.pca)       #
nclus.pca         <- nrow(connectivity.pca)          #
centers.pca       <- t(sapply(clusters.pca,function(x){
  x.sub <- Dims.pca[rownames(clusterLabels.pca[which(clusterLabels.pca[, x] == "1"),]),]
  return(colMeans(x.sub))
}))    
rownames(centers.pca) <- clusters.pca                                              # add row names
#Dims.pca              <- Dims.pca[ clusterLabels.pca %in% clusters.pca, ]          # 
#clusterLabels.pca     <- clusterLabels.pca[clusterLabels.pca %in% clusters.pca]
#
sds@lineages
xs <- c(NULL, centers.pca[, 1 ])
ys <- c(NULL, centers.pca[, 2 ])
zs <- c(NULL, centers.pca[, 3 ])
xs <- c( xs, as.numeric(sapply( sds@curves, function(c){     c$s[, 1 ] }) ))
ys <- c( ys, as.numeric(sapply( sds@curves, function(c){     c$s[, 2 ] }) ))
zs <- c( zs, as.numeric(sapply( sds@curves, function(c){     c$s[, 3 ] }) ))

# straight line 3d plot
rgl::plot3d(x = NULL, y = NULL, z = NULL, aspect = 'iso', 
            xlim = range(xs), ylim = range(ys), zlim = range(zs), 
            box=FALSE, axes=FALSE, xlab = '', ylab = '', zlab = '' )
colpal = rainbow(11)
rgl::plot3d(Dims.pca, col= rainbow(10)[cl], 
             add=TRUE, type='p', size=4, pch=20,alpha=I(1/8), 
            box=FALSE, axes=FALSE)
rgl::abclines3d(max(Dims.pca[,1]),max(Dims.pca[,2]),max(Dims.pca[,3]), 
                a = diag(3), col = "black", lwd=2)
rgl::plot3d(centers.pca, size = 10, add = TRUE, pch = 17,
            col = colpal[as.numeric(rownames(centers.pca))+1], alpha=1)
for (i in 1:(nclus.pca-1)){
  for (j in (i+1):nclus.pca){
    if (connectivity.pca[i,j]==1){
      rgl::lines3d(x=centers.pca[c(i,j),1], y=centers.pca[c(i,j),2], 
                   z=centers.pca[c(i,j),3], 
                   col=colpal[as.numeric(rownames(connectivity.pca))+1], 
                   lwd=2)
    }
  }
}
rgl::rgl.snapshot('Slingshot.Cluster9.branchs.straightLine2.png', fmt='png',top=TRUE)


```
```{r, out.width="50%", class.output="scroll-100"}
##############################################################################
#### draw multiple lineage inference ###
#############################################################################
## order the cells on each branch
library(dplyr)
for (i in 1:ncol(slingPseudotime(sds))){
  linedf <- data.frame(pseudotime=slingPseudotime(sds)[,i], 
                       clus.labels = cl, 
                       samples=rownames(slingPseudotime(sds)))
  linedf <- linedf[order(linedf$pseudotime),]
  medoids <- sapply(levels(linedf$clus.labels), function(clID){
    x.sub <- linedf %>% filter(linedf$clus.labels == clID) %>% dplyr::select(pseudotime)
    col <- colpal[as.numeric(as.character(linedf$clus.labels))+1][which.max(linedf$clus.labels == clID)]
    return( list(means=mean(x.sub$pseudotime, na.rm=TRUE), sdev= sd(x.sub$pseudotime, na.rm=TRUE),col=col) )
  })

  means = unlist(medoids$means)
  sdev  = unlist(medoids$sdev)
  col   = unlist(medoids$col)  
  #pdf(paste('2021.08.25.Slingshot.Cluster09.Pseudo.branch_',i,'.pdf',sep=''), 
   #   width=6,height=3)
  plot(linedf$pseudotime,rep(0, length(linedf$pseudotime)),cex=3,axes=F, 
       pch=16, xlab='', ylab='', 
       col=colpal[as.numeric(as.character(linedf$clus.labels))+1], 
       ylim=c(-0.1, 0.1), xlim = range(linedf$pseudotime, na.rm=TRUE)); 
  abline(h=0, col="black")
  points(x=means,y=rep(0.07, length(means)), col=col, pch=19)
  arrows(means-sdev, rep(0.07, length(means)), means+sdev, rep(0.07, length(means)), 
         length=0.05, angle=90, code=3, col=col)
  #dev.off()
}



```
