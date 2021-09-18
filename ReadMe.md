# Seurat_demo 
目的：Seurat 示例网站。
My another 单细胞代码库: https://github.com/DawnEve/scRNAseqCode

因为Seurat官方网站最近常常打不开，这里编写几个示例网页，方便代码查询。不定期更新。


## 更新方式
- 从具体repo下载Rmd和原始示例数据，在服务器上用Rnotebook Knit成html格式，这种格式把图片转为 base64 内置到html中。
- Chrome打开html后，查看页面显示是否符合预期。合理命名后收录到该 Repo 中。
- 更新版本信息 changelog.txt



## 现有页面


### 1. Seurat 基本示例

[seurat_3k_demo.html](./seurat_3k_demo.html)

来源：https://github.com/satijalab/seurat/tree/master/vignettes



### 2. sctransform 

In this vignette, we demonstrate how using sctransform based normalization enables recovering sharper biological distinction compared to log-normalization.

[sctransform.html](./sctransform.html) 

来源: https://raw.githubusercontent.com/satijalab/seurat/master/vignettes/sctransform_vignette.Rmd


### 3. cell cycle

dataset: https://github.com/DawnEve/Seurat_demo/issues/2

来源: https://github.com/satijalab/seurat/blob/master/vignettes/cell_cycle_vignette.Rmd

[cell_cycle_vignette](./cell_cycle_vignette.html)






