# Load packages

library(SeuratObject)
library(tibble)
library(metap)
library(Seurat)
library(hdf5r)
library(stringr)
library(dynwrap)
library(tidyverse)
library(plyr)
library(dplyr)
library(sctransform)
library(readr)
library(readxl)
library(DoubletFinder)
library(patchwork)
library(ggplot2)
library(dyno)
library(pheatmap)
library(cowplot)
library(ggsci)
library(paletteer)
library(RColorBrewer)
library(copykat)
library(heatmap3)
library(EnhancedVolcano)
library(org.Hs.eg.db)
library(scales) 
library(survival)
library(survminer)
library(SummarizedExperiment)
library(venn)
library(clusterProfiler)
library(msigdbr)
library(GSVA)
library(fgsea)
library(presto)
library(enrichplot)
library(circlize)
library(ComplexHeatmap)
library(ggstatsplot)
library(ggpubr)
library(ggthemes)
library(ggheatmap)
library(UCell)
library(irGSEA)
library(NMF)
library(Biobase)
library(harmony)
library(tidydr)
library(monocle)
library(SCENIC)
library(dorothea)
library(SCEVAN)
library(infercnv)
library(msigdf)
library(limma)
library(GseaVis)
library(viridis)
library(ggnewscale)
library(clustree)
library(ggalluvial)
library(scRepertoire)
library(grafify)
library(ggprism)
library(rstatix)
library(ggsignif)
library(reshape2)
library(CellChat)



## Initialize the Seurat object with the raw (non-normalized data)

MM_list <- list()

samples = basename(list.dirs("data", recursive = F))

for (sample in samples) {
  
  file =str_c("data/", sample)
  scrna_data <- Read10X(
    data.dir = file)
  
  seob <- CreateSeuratObject(
    counts = scrna_data,
    project = sample,
    min.cells = 3,
    min.features = 200)
  
  seob[['sample']] <- sample
  seob[['patient']] <- unlist(strsplit(sample, "_"))[1]
  seob[['period']] <- unlist(strsplit(sample, "_"))[2]

  MM_list[[sample]] = seob
}



## Merge Seurat Object

MM <- merge(x = MM_list[[1]],
            y = MM_list[-1],
            add.cell.ids = names(MM_list))



## QC and selecting cells for further analysis

MM[["percent.mt"]] = PercentageFeatureSet(object = MM, pattern = "^MT-")
MM[["percent.rb"]]<-PercentageFeatureSet(MM,pattern="^RP[SL]")

MM <- subset(MM,
             subset = nFeature_RNA > 200 &
               nFeature_RNA < 5000 &
               percent.mt < 10)

VlnPlot(MM, features = c("nCount_RNA"), ncol = 1,y.max = 25000, pt.size = 0,cols = pal)+  
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        #aspect.ratio = 1,
        strip.background = element_rect(colour = NA,fill = NA),  
        axis.text.x=element_text(angle=45,hjust = 1,vjust=1),
        plot.title = element_text(hjust = 0.5))

VlnPlot(MM, features = c("nFeature_RNA"), ncol = 1,y.max = 5000, pt.size = 0,cols = pal)+  
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        #aspect.ratio = 1,
        strip.background = element_rect(colour = NA,fill = NA),  
        axis.text.x=element_text(angle=45,hjust = 1,vjust=1),
        plot.title = element_text(hjust = 0.5))

VlnPlot(MM, features = c("percent.mt"), ncol = 1,y.max = 10, pt.size = 0,cols = pal)+  
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        #aspect.ratio = 1,
        strip.background = element_rect(colour = NA,fill = NA),  
        axis.text.x=element_text(angle=45,hjust = 1,vjust=1),
        plot.title = element_text(hjust = 0.5))

VlnPlot(MM, features = c("percent.rb"), ncol = 1,y.max = 40, pt.size = 0,cols = pal)+  
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        #aspect.ratio = 1,
        strip.background = element_rect(colour = NA,fill = NA),  
        axis.text.x=element_text(angle=45,hjust = 1,vjust=1),
        plot.title = element_text(hjust = 0.5))



## Normalizing the data

# split the dataset into a list of seurat objects

MM.list <- SplitObject(MM, split.by = "sample")

# Perform integration with SCTransform-normalized datasets

MM.list <- lapply(X = MM.list, FUN = SCTransform, method = "glmGamPoi")

features <- SelectIntegrationFeatures(object.list = MM.list, nfeatures = 3000)

MM.list <- PrepSCTIntegration(object.list = MM.list, anchor.features = features)

MM.list <- lapply(X = MM.list, FUN = RunPCA, features = features)

anchors <- FindIntegrationAnchors(object.list = MM.list, normalization.method = "SCT", 
                                  anchor.features = features, dims = 1:30, reduction = "rpca")

MM <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:30)

DefaultAssay(MM) <- "integrated"



## Perform PCA

MM <- RunPCA(MM)

ElbowPlot(MM, ndims = 50)



## UMAP visualization

MM <- RunUMAP(MM, dims = 1:30)

pal <- paletteer_d("ggsci::category20_d3")[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)]

DimPlot(MM,reduction = 'umap',group.by = "sample",label = F,pt.size = 0.1,cols = pal)+ 
  labs(title = 'Samples',size= 60)+  
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        aspect.ratio = 1,
        strip.background = element_rect(colour = NA,fill = NA),
        axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5))

DimPlot(MM,reduction = 'umap',group.by = "sample",split.by = "sample",label = F,pt.size = 0.1,cols = pal,ncol = 5)& 
  labs(title = 'Samples',size= 60)&  
  theme_bw(base_size = 14) &
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        aspect.ratio = 1,
        strip.background = element_rect(colour = NA,fill = NA),
        axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5))


pal <- paletteer_d("ggsci::nrc_npg")[c(1,2,3,4,5,6,7,8,9,10)]

DimPlot(MM,reduction = 'umap',group.by = "patient",label = F,pt.size = 0.1,cols = pal)+ 
  labs(title = 'Patients',size= 60)+  
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        aspect.ratio = 1,
        strip.background = element_rect(colour = NA,fill = NA),
        axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5))

DimPlot(MM,reduction = 'umap',group.by = "patient",split.by = "patient",label = F,pt.size = 0.1,cols = pal,ncol = 4)& 
  labs(title = 'Patients',size= 60)&  
  theme_bw(base_size = 14) &
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        aspect.ratio = 1,
        strip.background = element_rect(colour = NA,fill = NA),
        axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5))


pal <- paletteer_d("ggsci::category20c_d3")[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)]

DimPlot(MM,reduction = 'umap',group.by = "period",label = F,pt.size = 0.1,cols = pal)+ 
  labs(title = 'Periods',size= 60)+  
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        aspect.ratio = 1,
        strip.background = element_rect(colour = NA,fill = NA),
        axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5))

DimPlot(MM,reduction = 'umap',group.by = "period",split.by = "period",label = F,pt.size = 0.1,cols = pal,ncol = 2)& 
  labs(title = 'Periods',size= 60)&  
  theme_bw(base_size = 14) &
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        aspect.ratio = 1,
        strip.background = element_rect(colour = NA,fill = NA),
        axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5))



## Identify Doublets

sweep.res.list <- paramSweep_v3(MM, PCs = 1:50, sct = T)

sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE) 

bcmvn <- find.pK(sweep.stats)

pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()

DoubletRate = 0.1

homotypic.prop <- modelHomotypic(MM@meta.data$seurat_clusters)

nExp_poi <- round(DoubletRate*ncol(MM))

nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

MM <- doubletFinder_v3(MM, PCs = 1:50, pN = 0.25, pK = pK_bcmvn, 
                       nExp = nExp_poi.adj, reuse.pANN = F, sct = T)

head(MM@meta.data)

# Visualization of inferred Singlet and Doublet

MM@meta.data$DF.classifications_0.25_0.27_8892<-factor(MM@meta.data$DF.classifications_0.25_0.27_8892,levels = c("Singlet","Doublet"))

table(MM$DF.classifications_0.25_0.27_8892)

pal <- paletteer_d("ggsci::category20c_d3")[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)]

DimPlot(MM, reduction = "umap", group.by = "DF.classifications_0.25_0.27_8892",label = F,pt.size = 0.1,cols = pal)+ 
  labs(title = 'DoubletFinder',size= 60)+  
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        aspect.ratio = 1,
        strip.background = element_rect(colour = NA,fill = NA),
        axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5))

# Subset inferred Singlet

MM <- subset(MM, cells= rownames(MM@meta.data[MM@meta.data$DF.classifications_0.25_0.27_8892=="Singlet",]))



## Find cell clusters

MM <- FindNeighbors(MM,
                    dims = 1:50)

MM <- FindClusters(MM,
                   resolution =0.4,
                   random.seed = 1)

pal <- paletteer_d("ggsci::default_igv")[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,
                                           31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51)]

DimPlot(MM,reduction = 'umap',group.by = "seurat_clusters",label = T,pt.size = 0.1,cols = pal)+ 
  labs(title = 'Cell Clusters',size= 60)+  
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        aspect.ratio = 1,
        strip.background = element_rect(colour = NA,fill = NA),
        axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5))

DimPlot(MM,reduction = 'umap',group.by = "seurat_clusters",label = F,pt.size = 0.1,cols = pal)+ 
  labs(title = 'Cell Clusters',size= 60)+  
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        aspect.ratio = 1,
        strip.background = element_rect(colour = NA,fill = NA),
        axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5))



## Cell type annotation 

DefaultAssay(MM) <- "SCT"

Idents(MM) <- MM@meta.data$seurat_clusters

DotPlot(MM, features = c('CD3D','CD3E','CD3G','GNLY','NKG7',
                         'CD79A','CD79B','MZB1',"TNFRSF17",
                         'CD14','CD163','CD68','FCGR3A',"MPO","ELANE",
                         'CLEC10A','CD1C',"LILRA4","CLEC4C",
                         'PPBP','PF4'))+
  coord_flip()+
  theme_bw(base_size = 14)+  
  theme(panel.grid = element_blank())+ 
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))+  
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3)) 

cluster2type <- c(
  "0"="MNCs",
  "1"="MNCs",
  "2"="T cells",
  "3"="T cells",
  "4"="T cells",
  "5"="NKs",
  "6"="T cells",
  "7"="MNCs", 
  "8"="PCs/MMs", 
  "9"="GMPs", 
  "10"="T cells", 
  "11"="T cells", 
  "12"="T cells", 
  "13"="T cells", 
  "14"="GMPs", 
  "15"="B cells", 
  "16"="T cells", 
  "17"="MKs", 
  "18"="DCs", 
  "19"="pDCs")

MM[['cell_type']] = unname(cluster2type[MM@meta.data$seurat_clusters])

MM@meta.data$cell_type<-factor(MM@meta.data$cell_type,levels = c("T cells","NKs","B cells","PCs/MMs","MNCs","GMPs","DCs","pDCs","MKs"))

pal <- paletteer_d("ggsci::category20c_d3")[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)]

DimPlot(MM,reduction = 'umap',group.by = "cell_type",label = F,pt.size = 0.1,cols = pal)+ 
  labs(title = 'Cell Types',size= 60)+  
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        aspect.ratio = 1,
        strip.background = element_rect(colour = NA,fill = NA),
        axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5))

DimPlot(MM,reduction = 'umap',group.by = "cell_type",split.by = "period",label = F,pt.size = 0.1,cols = pal)+ 
  labs(title = 'Cell Types',size= 60)+  
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        aspect.ratio = 1,
        strip.background = element_rect(colour = NA,fill = NA),
        axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5))



## Add major cell type information

cluster2type <- c(
  "T cells"="Lymphocytes",
  "NKs"="Lymphocytes",
  "B cells"="Lymphocytes",
  "PCs/MMs"="Lymphocytes",
  "MNCs"="Myeloid cells",
  "GMPs"="Myeloid cells",
  "DC"="Myeloid cells", 
  "pDC"="Myeloid cells",
  "MKs"="Myeloid cells")

MM[['cell_type1']] = unname(cluster2type[MM@meta.data$cell_type])

DimPlot(MM,reduction = 'umap',group.by = "cell_type1",label = F,pt.size = 0.1,cols = pal)+ 
  labs(title = 'Cell Types',size= 60)+  
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        aspect.ratio = 1,
        strip.background = element_rect(colour = NA,fill = NA),
        axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5))



## Add minor cell type information

cluster2type <- c(
  "T cells"="T/NK",
  "NKs"="T/NK",
  "B cells"="B",
  "PCs/MMs"="B",
  "MNCs"="MNCs",
  "GMPs"="MNCs",
  "DC"="DCs", 
  "pDC"="DCs",
  "MKs"="MKs")

MM[['cell_type2']] = unname(cluster2type[MM@meta.data$cell_type])

MM@meta.data$cell_type2<-factor(MM@meta.data$cell_type2,levels = c("T/NK","B","MNCs","DCs","MKs"))

DimPlot(MM,reduction = 'umap',group.by = "cell_type2",label = F,pt.size = 0.1,cols = pal)+ 
  labs(title = 'Cell Types',size= 60)+  
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        aspect.ratio = 1,
        strip.background = element_rect(colour = NA,fill = NA),
        axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5))



## Save data

save(MM, file = "data/MM.rdata")