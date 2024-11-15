# Aim: Functions useful for scRNA-seq analysis
# ref: https://gitee.com/dawnEve/neu_colon/blob/master/base/functions.R

if(0){
	cat("v0.2 支持 html 二级标题，写法见上文")
	cat("v0.5 支持 html 一级标题，写法见上文")
}


#{**Colors**}#


##{**colorset**}##

# gene cluster 3
colors.cluster = c('#E64B35','#4DBBD5','#3C5488')


# 4 sample: 2 pairs
# colorset: AML025
colorset.AML = c(AML027="#F8766D", AML035="#FF0000", ctrl1="#A4DDDF", ctrl2="#00BFC4") #light before, dark after
barplot( rep(1,4), col = colorset.AML )


# 3 time points, 8 samples
colorset.sample=c( "orange", "orange4",
                   # "hotpink","deeppink", "red",  
                   "#ff9999", "deeppink", "#9F0429",
                   "cadetblue3","royalblue1", "navy")

# 3 time points;
colorset.time=c("0h"="#FFA500", "18h"="#FF1493", "48h"="#4876FF");
# scale_color_manual(values=c("#FFA500", "#FF1493", "#4876FF") ) #for 0h, 18h ,48h
colour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",  
         "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
         "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
         "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
barplot( rep(1, length(colour)), col = colour, border = NA, space = 0)


# 9 clusters
library(scales)
color.set1=scales::hue_pal()(9); color.set1
# "#F8766D" "#D39200" "#93AA00" "#00BA38" "#00C19F" "#00B9E3" "#619CFF" "#DB72FB" "#FF61C3"
barplot( rep(1, length(color.set1)), col = color.set1, border = NA, space = 0, yaxt="n")


# G1, S, G2M
colorset.cycle=c(G1="#3AC96D", S="#FDE725", G2M='deeppink')
barplot( rep(1,3), col = colorset.cycle )
# scale_color_manual(breaks = c("G1", "S", "G2M"), values=colorset.cycle)

# 3 color dot plot 
barplot( rep(1,3), col = c("#FDE725", "#3AC96D", "#440154") )


# 14 colors
colors <-c("#FED439FF","#709AE1FF",
           "#D5E4A2FF","#197EC0FF","#F05C3BFF","#46732EFF",
           "#71D0F5FF","#370335FF","#075149FF","#C80813FF","#91331FFF",
           "#1A9993FF","#FD8CC1FF", "deeppink")
barplot(rep(1,length(colors)), col=colors)





==> 获取 Seurat 中默认颜色配方/ 配色提取

1. DimPlot中默认的配色方案
library(scales)
show_col(hue_pal()(16))


2. 提取DimPlot中画聚类时用到的颜色
library(scales)
p1 <- DimPlot(pbmc.nmf, group.by = "celltype_assign")
x<-ggplot_build(p1)
info = data.frame(colour = x$data[[1]]$colour, group = x$data[[1]]$group)
info <- unique((arrange(info, group)))
cols <- as.character(info$colour)

> cols
 [1] "#F8766D" "#EC823C" "#DD8D00" "#CA9700" "#B3A000" "#97A900" "#71B000" "#2FB600" "#00BB4B" "#00BF76" "#00C098" "#00C0B7"
[13] "#00BDD1" "#00B7E8" "#00AEFA" "#3DA1FF" "#8F91FF" "#BE80FF" "#DE71F9" "#F265E7" "#FE61CF" "#FF64B3" "#FF6C92"






==> 查看颜色
scales::show_col( paletteer_d("ggsci::nrc_npg")[1:8] )















##{**show_colorset**}##

#' 一套颜色在单细胞点图的模拟效果
#' 
#' version:0.2 https://blog.csdn.net/wangjunliang/article/details/143275486
#' Seurat colorset: https://zhuanlan.zhihu.com/p/541666692 #2(9)
#'
#' @param colors2 颜色向量
#' @param pt.size  点的大小
#' @param dot.per.cluster 每个类产生颜色数量
#' @param cluster_number 主类大小，不设置，则表示和颜色总数一致
#' @param radius 主图半径，默认即可
#' @param zoom.factor 针对主图放大倍数
#' @param scale.factor 随机点沿着核心点间距的缩放倍数
#' @param shuffle 是否对颜色随机打乱，默认不随机
#'
#' @return 无返回值，就一个绘图效果
#' @export
#'
#' @examples
show_colorset=function(
  colors2,
  pt.size=1,
  dot.per.cluster=100,
  cluster_number=0,
  radius=10,
  zoom.factor=2, #绘制核心点时，整体放大倍数，方便个后续随机点留下空间
  scale.factor=6,
  shuffle=F,
  main=""
    ){
  
  message(length(colors2))
  if(cluster_number<=0){
    cluster_number = length(colors2)
  }
  
  if(shuffle){
    colors2=sample(colors2)
  }
  
  #1.确定几个核心点
  arr_x= radius * cos(2*pi / cluster_number* (1:cluster_number))
  arr_y= radius * sin(2*pi / cluster_number* (1:cluster_number))
  
  #2.计算两点的距离
  dot_dist = sqrt( (arr_x[1]-arr_x[2])**2 +  (arr_y[1]-arr_y[2])**2); dot_dist
  
  
  #3.噪音点，随机分布在核心点周围，距离大概是：核心点距离/scale.factor
  noiseX = dot_dist*rnorm(n=dot.per.cluster)/scale.factor
  noiseY = dot_dist*rnorm(n=dot.per.cluster)/scale.factor
  
  #3. 绘制空坐标轴
  main=ifelse(main=="", "Color test", main)
  plot(arr_x*zoom.factor, arr_y*zoom.factor, col="white", xlab="UMAP_1", ylab="UMAP_2", main=main, mgp=c(2,1,0))
  #4. 绘制噪音点
  for(i in 1:cluster_number){
    points(arr_x[i] + sample(noiseX), arr_y[i]+sample(noiseY), col=colors2[i], pch=19, cex=pt.size)
  }
}
if(0){
  show_colorset( DiscretePalette(26, palette='alphabet')[1:10], dot.per.cluster=500, zoom.factor = 1.2, pt.size = 2, main="alphabet")
  show_colorset( c("red", "orange", "blue", "navy", "cyan", "grey"), dot.per.cluster=2000, zoom.factor = 2)
  show_colorset( c("red", "orange", "blue", "navy", "cyan", "grey"), dot.per.cluster=2000, zoom.factor = 2, shuffle = T )
}











#{**1. Loading Data**}#

##{**load 10x with prefix**}##

# Part1: 非标准格式，cell barcode 和 gene 都有header
dat.data=Matrix::readMM(file = "/data/wangjl/scPolyA-seq2/rawData/NG2022/matrix.mtx.gz")
cell.barcodes <- as.data.frame(data.table::fread("/data/wangjl/scPolyA-seq2/rawData/NG2022/barcodes.tsv.gz", header = TRUE))
feature.names <- as.data.frame(data.table::fread("/data/wangjl/scPolyA-seq2/rawData/NG2022/features.tsv.gz", header = TRUE))

dim(dat.data) #[1] 713403  23537
dim(cell.barcodes) #713403     31
dim(feature.names) #23537    15

rownames(cell.barcodes)=cell.barcodes$index

dat=Matrix::t(dat.data)
dim(dat) #23537 713403
rownames(dat)=feature.names[,1] #使用第一列 feature name
colnames(dat)=cell.barcodes[,1] #使用第一列作为 cell id
dat[1:2, 1:3]

scCD4T <- CreateSeuratObject(counts = dat, project = "CD4T", meta.data = cell.barcodes, min.cells = 3, min.features = 200)



# Part2: 标准格式，仅是加了前缀 prefix，cell bacode 和 gene 都没有header:
Read10X_2=function (data.dir, gene.column = 2, cell.column = 1, unique.features = TRUE, 
          strip.suffix = FALSE, prefix="") 
{
  full.data <- list()
  has_dt <- requireNamespace("data.table", quietly = TRUE) && 
    requireNamespace("R.utils", quietly = TRUE)
  for (i in seq_along(along.with = data.dir)) {
    run <- data.dir[i]
    if (!dir.exists(paths = run)) {
      stop("Directory provided does not exist")
    }
    barcode.loc <- file.path(run, paste0(prefix, "barcodes.tsv") )
    gene.loc <- file.path(run, paste0(prefix, "genes.tsv"))
    features.loc <- file.path(run, paste0(prefix, "features.tsv.gz"))
    matrix.loc <- file.path(run, paste0(prefix, "matrix.mtx"))
    pre_ver_3 <- file.exists(gene.loc)
    if (!pre_ver_3) {
      addgz <- function(s) {
        return(paste0(s, ".gz"))
      }
      barcode.loc <- addgz(s = barcode.loc)
      matrix.loc <- addgz(s = matrix.loc)
    }
    
    if (!file.exists(barcode.loc)) {
      stop("Barcode file missing. Expecting ", basename(path = barcode.loc))
    }
    if (!pre_ver_3 && !file.exists(features.loc)) {
      stop("Gene name or features file missing. Expecting ", 
           basename(path = features.loc))
    }
    if (!file.exists(matrix.loc)) {
      stop("Expression matrix file missing. Expecting ", 
           basename(path = matrix.loc))
    }
    data <- readMM(file = matrix.loc)
    if (has_dt) {
      cell.barcodes <- as.data.frame(data.table::fread(barcode.loc, 
                                                       header = FALSE))
    }
    else {
      cell.barcodes <- read.table(file = barcode.loc, header = FALSE, 
                                  sep = "\t", row.names = NULL)
    }
    if (ncol(x = cell.barcodes) > 1) {
      cell.names <- cell.barcodes[, cell.column]
    }
    else {
      cell.names <- readLines(con = barcode.loc)
    }
    if (all(grepl(pattern = "\\-1$", x = cell.names)) & strip.suffix) {
      cell.names <- as.vector(x = as.character(x = sapply(X = cell.names, 
                                                          FUN = ExtractField, field = 1, delim = "-")))
    }
    if (is.null(x = names(x = data.dir))) {
      if (length(x = data.dir) < 2) {
        colnames(x = data) <- cell.names
      }
      else {
        colnames(x = data) <- paste0(i, "_", cell.names)
      }
    }
    else {
      colnames(x = data) <- paste0(names(x = data.dir)[i], 
                                   "_", cell.names)
    }
    if (has_dt) {
      feature.names <- as.data.frame(data.table::fread(ifelse(test = pre_ver_3, 
                                                              yes = gene.loc, no = features.loc), header = FALSE))
    }
    else {
      feature.names <- read.delim(file = ifelse(test = pre_ver_3, 
                                                yes = gene.loc, no = features.loc), header = FALSE, 
                                  stringsAsFactors = FALSE)
    }
    if (any(is.na(x = feature.names[, gene.column]))) {
      warning("Some features names are NA. Replacing NA names with ID from the opposite column requested", 
              call. = FALSE, immediate. = TRUE)
      na.features <- which(x = is.na(x = feature.names[, 
                                                       gene.column]))
      replacement.column <- ifelse(test = gene.column == 
                                     2, yes = 1, no = 2)
      feature.names[na.features, gene.column] <- feature.names[na.features, 
                                                               replacement.column]
    }
    if (unique.features) {
      fcols = ncol(x = feature.names)
      if (fcols < gene.column) {
        stop(paste0("gene.column was set to ", gene.column, 
                    " but feature.tsv.gz (or genes.tsv) only has ", 
                    fcols, " columns.", " Try setting the gene.column argument to a value <= to ", 
                    fcols, "."))
      }
      rownames(x = data) <- make.unique(names = feature.names[, 
                                                              gene.column])
    }
    if (ncol(x = feature.names) > 2) {
      data_types <- factor(x = feature.names$V3)
      lvls <- levels(x = data_types)
      if (length(x = lvls) > 1 && length(x = full.data) == 
          0) {
        message("10X data contains more than one type and is being returned as a list containing matrices of each type.")
      }
      expr_name <- "Gene Expression"
      if (expr_name %in% lvls) {
        lvls <- c(expr_name, lvls[-which(x = lvls == 
                                           expr_name)])
      }
      data <- lapply(X = lvls, FUN = function(l) {
        return(data[data_types == l, , drop = FALSE])
      })
      names(x = data) <- lvls
    }
    else {
      data <- list(data)
    }
    full.data[[length(x = full.data) + 1]] <- data
  }
  list_of_data <- list()
  for (j in 1:length(x = full.data[[1]])) {
    list_of_data[[j]] <- do.call(cbind, lapply(X = full.data, 
                                               FUN = `[[`, j))
    list_of_data[[j]] <- as.sparse(x = list_of_data[[j]])
  }
  names(x = list_of_data) <- names(x = full.data[[1]])
  if (length(x = list_of_data) == 1) {
    return(list_of_data[[1]])
  }
  else {
    return(list_of_data)
  }
}

if(0){
	txt_filenames=dir("/home/wangjl/data4/others/hanlu/raw/GSE274187/", pattern = "*.gz"); txt_filenames
	getGSM=function(txt_filenames){
	  gsm=c()
	  for(i in 1:length(txt_filenames)){
		# ignore TCR
		if( 0 != length( grep("contig", txt_filenames[i], value=T) ) ){
		  next;
		}
		t1=strsplit(txt_filenames[i], "_")[[1]]
		gsm=c(gsm, t1[1])
	  }
	  gsm
	}
	gsm.arr=getGSM(txt_filenames) |> unique()
	gsm.arr #19 total; 
	gsm_index=1; 
	gsm = gsm.arr[gsm_index];
	title=paste0(strsplit(grep(gsm.arr[gsm_index], txt_filenames, value=T)[1], "_")[[1]][1:2], collapse = "_")
	
	dat=Read10X_2("/home/wangjl/data4/others/hanlu/raw/GSE274187/", prefix = paste0(title, "_"))
	names(dat) #"Gene Expression"  "Antibody Capture" #两个矩阵：RNA和 22个蛋白
	
	
	scRNA=CreateSeuratObject(counts = dat$`Gene Expression`, project = title)
}

# https://gitee.com/dawnEve/others/raw/master/hanlu/script/a02_BCMA-2B_load.Rmd






##{**load h5 format**}##
h5_filename="/datapool/wangjl/others/hanlu/raw/GSE210079/GSM6459763_32-3mo_raw_feature_bc_matrix.h5"
dat=Read10X_h5(h5_filename)
str(dat)
names(dat) #"Gene Expression"  "Antibody Capture" #两个矩阵：RNA和 55个蛋白
str(dat$`Gene Expression`)

# make sure cell id are the same
all.equal(colnames(dat[["Gene Expression"]]), colnames(dat[["Antibody Capture"]])) #T

# (2). use RNA data to create Obj
scRNA=CreateSeuratObject(counts = dat$`Gene Expression`, project = "A1")

# (3). add protein mat
# https://zhuanlan.zhihu.com/p/567253121
adt_assay <- CreateAssayObject(counts = dat$`Antibody Capture`)
scRNA[["ADT"]] <- adt_assay

# https://gitee.com/dawnEve/others/blob/master/hanlu/script/a01_BCMA-1_load.h5.R

















#{**2. Main steps**}#



##{**Feature_scale_PCA**}##
# 1. HVG, Scale(all), PCA
Feature_scale_PCA=function(object, nfeatures, scale.all=F){
  object <- FindVariableFeatures(object = object, nfeatures=nfeatures)
  if(scale.all){
	genelist=rownames(object)
  }else{
	genelist=VariableFeatures(object)
  }
  object <- ScaleData(object = object, features = genelist)
  object <- RunPCA(object=object, features=VariableFeatures(object))
  return(object)
}


##{**Neighbor_cluster_umap**}##
# 2. cellCluster,UMAP
Neighbor_cluster_umap = function(object, dims, resolution){
  object <- FindNeighbors(object = object, dims=dims)
  object <- FindClusters(object = object, resolution = resolution)
  object <- RunUMAP(object = object, dims=dims)
  return(object)
}






##{**addQC**}##
AddQC=function(obj){
	obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
	obj[["percent.rp"]] <- PercentageFeatureSet(obj, pattern = "^RP[SL]")
	obj[["percent.hb"]] <- PercentageFeatureSet(obj, pattern = "^HB[^(P)]")
	obj
}

VlnPlot(obj, 
        features = c("nFeature_RNA",
                     "nCount_RNA", 
                     "percent.mt",
                     "percent.rp",
                     "percent.hb"),
        ncol = 3,pt.size = 0.1, group.by = "orig.ident")
#




##{**FeaturePlot2**}##

library(ggplot2)

FeaturePlot2=function(obj, features, cols = c("lightgrey", 'red'), ncol=1, ...){
  FeaturePlot(obj, features =features, 
              cols=cols,
              ncol = ncol, ... ) & NoLegend() & NoAxes() & theme(
                panel.border = element_rect(color = "black", size = 1)
              )
}
# 效果 https://blog.csdn.net/wangjunliang/article/details/134716690
if(0){
	FeaturePlot2(scObj, features =c(
	  "CD3D", "CD4", "CD8A", "CD79A", "NKG7", "GZMB", "CD14", "CD68", "CCR7", 
	  "KRT8", "EPCAM", "COL1A2", "COL3A1", "VWF", "MKI67", "CCNB1", "CD34", "SOCS2"
	), ncol=5)
}









##{**DimPlot2**}##

#' DimPlot with 缩小的坐标轴
#'
#' @param scObject 
#' @param reduction 
#' @param group.by 
#' @param label 
#' @param raster 
#' @param legend.position 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
DimPlot2=function(scObject, reduction="umap", group.by = "seurat_clusters", 
                  label=F, raster=F, legend.position="right", #"none"
                  ...){
  # https://stackoverflow.com/questions/78667978/plotting-only-half-length-axis-lines
  #1) init plot
  p1 <- DimPlot(scObject, 
                reduction = reduction, 
                group.by = group.by, 
                label = label, 
                raster = raster)#+ 
  #scale_colour_manual(values = group2.cols) #+ 
  #labs(title = "10x RNA", x = "UMAP_1", y = "UMAP_2")
  #2) get range
  getRange=function(x){
    min(x) + 0.25 * diff(x)
  }
  #3) set range
  p1 + 
    scale_x_continuous(breaks = getRange(p1$data[,1]), guide = guide_axis(cap = 'upper')) +
    #scale_y_continuous(breaks = quantile(p1$data[,2], prob = 0.20), guide = guide_axis(cap = 'upper')) +
    scale_y_continuous(breaks = getRange(p1$data[,2]), guide = guide_axis(cap = 'upper')) +
    
    theme(aspect.ratio = 1,
      panel.border = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(arrow = arrow(type = "closed", length = unit(0.2, "cm"))),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_text(hjust = 0.025),
      legend.position = legend.position
    )
}

if(0){
	DimPlot2(scObj, reduction = "umap", label=T)
	DimPlot2(scObj, reduction = "umap", label=T, legend.position = "none")
}

















##{**ggplot theme**}##

# large dot in legend
LargeLegend = function(size=5){
  guides(colour = guide_legend(override.aes = list(alpha = 1,size=size)))
}

# add panel box for ggplot
AddBox = function(...){
  theme(
    panel.border = element_rect(color = "black", size = 1),
    validate = TRUE, 
    ...
  )
}










##{**StackVlnPlot**}##

#' Stack VlnPlot
#'
#' @param obj Seurat object
#' @param features gene list
#' @param flip default: cluster on x, gene symbol on y and right side
#' @param ... 
#'
#' @return ggplot2 obj
#' @export
#'
#' @examples
StackVlnPlot=function(obj, features, flip=T, ...){
  VlnPlot(obj, features = features, pt.size = 0,
          stack =T, flip = flip, ...) &
    theme(
      legend.position = "none",
      axis.ticks.y = element_blank(),
      axis.text.y = element_text(size=10, vjust = 1, hjust = 1),
      strip.text.y = element_text(size=12, hjust = 0),
    )
}








##{**ShowCluster**}##

#' Show the position of a cluster on umap
#'
#' @param scObj 
#' @param clusterId 
#' @param slot default seurat_clusters
#' @param color highlight color, default darkred
#'
#' @return
#' @export
#'
#' @examples
#' showCluster(scObj, 17)
#' showCluster(scObj, "BL1Y", slot="sample")
#' showCluster(scObj, "BL1Y", slot="sample", color="purple")
ShowCluster = function(scObj, clusterId, slot="seurat_clusters", color="darkred"){
  df2 = as.data.frame(scObj@meta.data);
  DimPlot(scObj, label = F, #group.by = "sample", 
          #cells.highlight = df2[which(df2[, slot] == clusterId), "cell"],
          cells.highlight = rownames( df2[which(df2[, slot] == clusterId), ]),
          cols.highlight = color, cols = "grey")+ #c("darkred", "darkblue")
    labs(title=paste0(slot, ": ", clusterId))+
    theme(
      legend.position = "none"
    )
}














#{**3. Adv funcs **}#



##{**table2barplot**}##

#' draw barplot, from a data.frame generated by table(para1, para2), colored by 1st parameter of table()
#' v2.1
#' v2.2 第一参数的空值跳过
#' v2.3 legendTitle
#'
#' @param tbl1 data.frame, draw by each column
#' @param colors colors of each row(default NULL, auto-color)
#' @param scale whether scale to 1 or not(default T)
#' @param title the main title of the figure,(is main in the function)
#' @param legendY the position y of legend, adjust as needed, default -0.25
#' @param omit whether to omit some columns
#' @param ... 
#'
#' @return NULL
#' @export
#'
#' @examples
#' 
table2barplot=function(tbl1, colors=NULL,levels=NULL, scale=T, title="", 
                       omit=NULL, xlab="", ylab="", legendTitle=NULL, width=0.9){
  tbl1= tbl1[, which(colSums(tbl1)>0)] #remove all 0 columns
  tbl1= tbl1[which(rowSums(tbl1)>0), ] #remove all 0 rows
  # remove some columns by column names
  if(!is.null(omit)){
    tbl1=tbl1[, setdiff( colnames(tbl1), as.character( omit ) ) ]
  }
  # table to data.frame(wide to long)
  df2=as.data.frame(tbl1)
  if(!is.null(levels)){
    df2$Var1=factor(df2$Var1, levels =levels )#change order
  }
  if(""==ylab){ ylab=ifelse(scale, "Freq", "Count") }
  if(""==xlab){ xlab="Index" }
  
  # draw
  g1=ggplot(df2, aes(x=Var2, y=Freq, fill=Var1))+
    geom_bar(stat="identity", position=ifelse(scale, "fill","stack"), width=width )+
    labs(x=xlab, y=ylab, title=title)+
    theme_classic(base_size = 14)+
    theme(axis.text.x=element_text(angle=60, hjust=1,size=rel(1.2)) )
  
  legendTitle=ifelse(is.null(legendTitle), "", legendTitle)
  if(is.null(colors)){
    return (g1 + scale_fill_discrete( legendTitle ))
  }else{
    return( g1+scale_fill_manual(legendTitle, values = colors) )
  }
}
if(0){
  table2barplot(
    as.table(t(
      apply(
        table(scObj$time, scObj$seurat_clusters), 
        1, 
        function(x){ x/sum(x) * 1e4})
    )),
    colors = c("0h"="#FFA500", "18h"="#FF1493", "48h"="#4876FF"),
    title="time"
  )
}









##{**VolcanoPlot**}##

#' Classical VolcanoPlot
#' @version v0.2 add more para to set: fig title, x axis label
#' @version v0.3 change parameter xlab to xlab2, to avoid ggplot2::xlab() | ggplot2 3.4.1
#'
#' @param dif a data frame with 3 column: symbol, log2FoldChange, padj
#' @param log2FC threshold, default log2(1.5)
#' @param padj threshold, default 0.05
#' @param label.symbols label dots on the plot, default NA
#' @param label.max label how many dots, default max=30
#' @param cols colors for dots: sig-down(blue), sig-up(red)
#' @param title title of the fig, default ""
#' @param title.base first part of title, default "DEG", then :
#' @param xlab2 x axis label, default NA
#'
#' @return ggplot obj
#' @export
#'
#' @examples
VolcanoPlot=function(dif, log2FC=log2(1.5), padj=0.05, 
                     label.symbols=NULL, label.max=30,
                     cols=c("#497aa2", "#ae3137"), title="",title.base="DEG", xlab2=NA){
  if( all( !c("log2FoldChange", "padj", "symbol") %in% colnames(dif) )){
    stop("Colnames must include: log2FoldChange, padj, symbol")
  }
  rownames(dif)=dif$symbol
  
  # (1) define up and down
  dif$threshold="ns";
  dif[which(dif$log2FoldChange > log2FC & dif$padj < padj),]$threshold="up";
  dif[which(dif$log2FoldChange < (-log2FC) & dif$padj < padj),]$threshold="down";
  dif$threshold=factor(dif$threshold, levels=c('down','ns','up'))
  #head(dif)
  #
  tb2=table(dif$threshold); print(tb2)
  library(ggplot2)
  # (2) plot
  g1 = ggplot(data=dif, aes(x=log2FoldChange, y=-log10(padj), color=threshold)) +
    geom_point(alpha=0.8, size=0.8) +
    geom_vline(xintercept = c(-log2FC, log2FC), linetype=2, color="grey")+
    geom_hline(yintercept = -log10(padj), linetype=2, color="grey")+
    labs(title= ifelse(""==title, "", paste(title.base,":", title)))+
    xlab(if(is.na(xlab2)) bquote(Log[2]*FoldChange) else xlab2 )+
    ylab(bquote(-Log[10]*italic(P.adj)) )+
    theme_classic(base_size = 14) +
    theme(legend.box = "horizontal",
          legend.position="top",
          legend.spacing.x = unit(0, 'pt'),
          legend.text = element_text( margin = margin(r = 20) ),
          legend.margin=margin(b= -10, unit = "pt"),
          plot.title = element_text(hjust = 0.5, size=10)
    ) +
    scale_color_manual('',labels=c(paste0("down(",tb2[[1]],')'),'ns',
                                   paste0("up(",tb2[[3]],')' )),
                       values=c(cols[1], "grey", cols[2]) )+
    guides(color=guide_legend(override.aes = list(size=3, alpha=1))); g1;
  # (3)label genes
  if(is.null(label.symbols)){
    dif.sig=dif[which(dif$threshold != "ns" ), ]
    len=nrow(dif.sig)
    if(len<label.max){
      label.symbols=rownames(dif.sig)
    }else{
      dif.sig=dif.sig[order(dif.sig$log2FoldChange), ]
      dif.sig= rbind(dif.sig[1:(label.max/2),], dif.sig[(len-label.max/2):len,])
      label.symbols=rownames(dif.sig)
    }
  }
  dd_text = dif[label.symbols, ]
  print(dd_text)
  # add text
  library(ggrepel)
  g2=g1 + geom_text_repel(data=dd_text,
                       aes(x=log2FoldChange, y=-log10(padj), label=row.names(dd_text)),
                       #size=2.5, 
                       colour="black",alpha=1)
  g2
}


if(0){
  # 0.
  deg_all=FindMarkers(scObj, ident.1 = "DSS", ident.2 = "WT", group.by="origin", min.pct = 0.001)
  dim(deg_all) #2137    5
  head(deg_all)
  
  #1. data frame with 3 columns: symbol, log2FoldChange, padj
  dif=data.frame(
    symbol=rownames(deg_all),
    log2FoldChange=deg_all$avg_log2FC,
    padj=deg_all$p_val_adj
  )
  
  #2. plot: give DEG max number, select genes order by FC max and min
  VolcanoPlot(dif, padj=0.05, title="DSS vs WT", label.max = 50)
  # dot color for down and up genes
  VolcanoPlot(dif, padj=0.05, title="DSS vs WT", label.max = 50, cols=c("blue", "red"))
  
  
  # Or spicify gene names
  VolcanoPlot(dif, padj=1e-10, title="DSS vs WT -2", 
              label.symbols=dif[ ((abs(dif$log2FoldChange) > 2) & (dif$padj < 1e-50) ) | 
                                   abs(dif$log2FoldChange) > 4,]$symbol )
}












##{**multiVolcanoPlot**}##

# > head(DEG)
#              p_val avg_log2FC pct.1 pct.2     p_val_adj     cluster  gene label
#RPS12 1.273332e-143  0.7298951 1.000 0.991 1.746248e-139 Naive CD4 T RPS12    Up
#RPS6  6.817653e-143  0.6870694 1.000 0.995 9.349729e-139 Naive CD4 T  RPS6    Up
#


color.pals = c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
               "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
               "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
               "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

#
#' multi volcano plot for scRNA-seq
#' @version 0.2 change legend order
#' @version 0.3 add max_overlaps for annotation
#' @version 0.4 label.size = NA,      # No border
#'
#' @param dat Seurat FindAllMarkers returns, must set only.pos = F;
#' @param color.arr color list, default same as Seurat
#' @param onlyAnnotateUp only annote gene symbols for up genes
#' @param log2Foldchang threshold for annotation
#' @param adjp  threshold for annotation
#' @param top_marker gene number for annotation
#' @param max_overlaps annotation label overlapping
#'
#' @return ggplot2 obj
#' @export
#'
#' @examples
multiVolcanoPlot = function(dat, color.arr=NULL, onlyAnnotateUp=T,
                            log2Foldchang=0.58, adjp=0.05, top_marker=5, 
                            max_overlaps=10, width=0.9){
  library(dplyr)
  library(ggrepel)
  # set default color list
  if(is.null(color.arr)){
    len = length(unique(dat$cluster))
    color.arr=scales::hue_pal()(len)
  }
  
  dat.plot = dat %>% mutate(
    "significance"=case_when(p_val_adj < adjp & avg_log2FC >= log2Foldchang  ~ 'Up',
                             p_val_adj < adjp & avg_log2FC <= -log2Foldchang  ~ 'Down',
                             TRUE ~ 'None'))
  tbl = table(dat.plot$significance)
  print( tbl )
  background.dat = data.frame(
    dat.plot %>% group_by(cluster) %>% filter(avg_log2FC>0) %>%
      summarise("y.localup"=max(avg_log2FC)),
    dat.plot %>% group_by(cluster) %>% filter(avg_log2FC<=0) %>%
      summarise("y.localdown"=min(avg_log2FC)),
    x.local=seq(1:length(unique(dat.plot$cluster)))
  ) %>% select(-cluster.1)
  #names(background.dat)
  #head(background.dat)
  #dim(background.dat)
  
  #
  x.number = background.dat %>% select(cluster, x.local)
  dat.plot = dat.plot%>% left_join(x.number,by = "cluster")
  #names(dat.plot)
  #head(dat.plot)
  
  #selecting top-up and top-down proteins
  dat.marked.up = dat.plot %>% filter(significance=="Up") %>%
    group_by(cluster) %>% arrange(-avg_log2FC) %>%
    top_n(top_marker,abs(avg_log2FC))
  dat.marked.down = dat.plot %>% filter(significance=="Down") %>%
    group_by(cluster) %>% arrange(avg_log2FC) %>%
    top_n(top_marker,abs(avg_log2FC))
  dat.marked = dat.marked.up %>% bind_rows(dat.marked.down)
  #referring group information data
  dat.infor = background.dat %>%
    mutate("y.infor"=rep(0,length(cluster)))
  #names(dat.infor)
  #dim(dat.infor)
  #head(dat.infor)
  
  ##plotting:
  #setting color by loading local color schemes
  vol.plot = ggplot()+
    # background
    geom_col(background.dat,mapping=aes(x.local, y.localup),
             fill="grey80", alpha=0.2, width=0.9, just = 0.5)+
    geom_col(background.dat,mapping=aes(x.local,y.localdown),
             fill="grey80", alpha=0.2, width=0.9, just = 0.5)+
    # point plot
    geom_jitter(dat.plot, mapping=aes(x.local, avg_log2FC, #x= should be number, Not string or factor
                                      color=significance),
                size=0.8, width = 0.4, alpha= 1)+
    scale_color_manual(name="significance", 
                       breaks = c('Up', 'None', 'Down'),
                       values = c("#d56e5e","#cccccc", "#5390b5")) + #set color for: Down None   Up
    geom_tile(dat.infor, mapping=aes(x.local, y.infor), #x axis color box
              height = log2Foldchang*1.3,
              fill = color.arr[1:length(unique(dat.plot$cluster))],
              alpha = 0.5,
              width=width) +
    labs(x=NULL,y="log2 Fold change")+
    geom_text(dat.infor, mapping=aes(x.local,y.infor,label=cluster))+
    # Down is not recommend, not meaningful, hard to explain; so prefer dat.marked.up to dat.marked
    ggrepel::geom_label_repel(data=if(onlyAnnotateUp) dat.marked.up else dat.marked, #gene symbol, of up group default
                              mapping=aes(x=x.local, y=avg_log2FC, label=gene),
                              force = 2, #size=2,
                              max.overlaps = max_overlaps,
                              label.size = NA, #no border
                              fill="#00000000", #box fill color
                              seed = 233,
                              min.segment.length = 0,
                              force_pull = 2,
                              box.padding = 0.1,
                              segment.linetype = 3,
                              #segment.color = 'black',
                              #segment.alpha = 0.5,
                              #direction = "x", #line direction
                              hjust = 0.5)+
    annotate("text", x=1.5, y=max(background.dat$y.localup)+1,
             label=paste0("|log2FC|>=", log2Foldchang, " & FDR<", adjp))+
    theme_classic(base_size = 12)+
    
    theme(
      axis.title = element_text(size = 13, color = "black"),
      axis.text = element_text(size = 15, color = "black"),
      axis.line.y = element_line(color = "black", size = 0.8),
      #
      axis.line.x = element_blank(), #no x axis line
      axis.ticks.x = element_blank(), #no x axis ticks
      axis.title.x = element_blank(), #
      axis.text.x = element_blank(),
      #
      legend.spacing.x = unit(0.1,'cm'),
      legend.key.width = unit(0.5,'cm'),
      legend.key.height = unit(0.5,'cm'),
      legend.background = element_blank(),
      legend.box = "horizontal",
      legend.position = c(0.13, 0.77),legend.justification = c(1,0)
    )+
    guides( #color = guide_legend( override.aes = list(size=5) ), #legend circle size
      color=guide_legend( override.aes = list(size=5), title="Change")
    )
  #guides(fill=guide_legend(title="Change"))+ #change legend title
  vol.plot
}
if(0){
	#multiVolcanoPlot(DEG, color.pals)
	multiVolcanoPlot(scObj.markers.time)
	multiVolcanoPlot(scObj.markers.time, onlyAnnotateUp = F)
}
