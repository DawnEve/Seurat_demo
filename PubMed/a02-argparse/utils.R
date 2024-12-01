library(RColorBrewer)
sample_cols = unique(c(brewer.pal(8,"Dark2"),brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),brewer.pal(12,"Set3")))
defined_cols = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', sample_cols)

#qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#qual_col_pals = qual_col_pals[c(3,1,2,4,5,6,7,8),]
#qual_cols = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

Read_SC <- function(input.file, 
	input.format="h5", 
	assay="RNA", 
	GENENAME=NULL, 
	image.file=NULL, 
	image.name="image",
	use.names=FALSE, 
	filter.count=FALSE, 
	project="SeuratProject", ...
){
	if(!is.null(GENENAME)){
		use.names=FALSE
	}
	## input.format: h5, 10x, csvd, csv
	if(input.format == "h5"){
		reader_h5 = Read10X_h5(input.file, use.names = use.names)
		count_data = as(reader_h5, "matrix")
	}else if(input.format == "10x" | input.format == "10X"){
		reader_h5 = Read10X_h5(file.path(input.file,'filtered_feature_bc_matrix.h5'), use.names = use.names)
		count_data = as(reader_h5, "matrix")
		if(is.null(image.file)){
			image.file = file.path(input.file,"spatial")
		}
	}else if(input.format == "csvd"){
		gene.col = 1
		if(use.names) {gene.col = 2}
		reader_10x = Read10X(data.dir=input.file, gene.column = gene.col)
		count_data = as(reader_10x, "matrix")
	}else if(input.format == "csv"){
		count_data = read.file(input.file, row.names=1, check.names = FALSE, sep=",", head=TRUE)
	}else if(input.format == "tsv"){
		count_data = read.file(input.file, row.names=1, check.names = FALSE, sep="\t", head=TRUE)
	}
	if(!is.null(GENENAME)){
		rownames(count_data) = map.genename(rownames(count_data), GENENAME)
	}
	if(filter.count){
		indices = which(rowSums(count_data)>0)
		count_data = count_data[indices,]
	}
	if(assay == "RNA"){
		seurat.obj = CreateSeuratObject(counts = count_data, assay = assay, project = project, ...)
	}else{
		if(is.null(image.file)){
			stop("image file need for Spatial transcriptome!!!")
		}
		seurat.obj = CreateSeuratObject(counts = count_data, assay = assay, project = project, ...)
		seurat.obj$slice = '1'
		seurat.obj$region = "region"
		image.obj = Read10X_Image(image.file, filter.matrix = TRUE)
		DefaultAssay(object = image.obj) = 'Spatial'
		image.obj = image.obj[colnames(seurat.obj)]
		seurat.obj[[image.name]] = image.obj
	}
	return(seurat.obj)
}

read.file <- function(filename, ...){
	if(grepl(".gz$",filename)){
		return(read.table(gzfile(filename), ...))
	}else{
		return(read.table(filename, ...))
	}
}

Read_genename <- function(genename.file, geneid.col=1, genename.col=2){
	dat = read.file(genename.file, sep="\t", colClasses="character")
	indices = order(dat[,geneid.col])
	GENENAME = data.frame(gid = dat[indices, geneid.col], gname = dat[indices, genename.col], stringsAsFactors=FALSE)
	GENENAME$uname = make.unique(GENENAME$gname)
	return(GENENAME)
}

Read_GeneName <- function(input.file, input.format, genename.file=NULL){
	if(!is.null(genename.file)){
		GENENAME = Read_genename(genename.file)
		return (GENENAME)
	}
	if(input.format == "h5"){
		h5.infile <- hdf5r::H5File$new(filename = input.file, mode = "r")
		GENENAME = data.frame(gid = h5.infile[["matrix/features/id"]][], 
			gname = h5.infile[["matrix/features/name"]][], stringsAsFactors=FALSE)
		GENENAME = GENENAME[order(GENENAME$gid),]
		GENENAME$uname = make.unique(GENENAME$gname)
		return (GENENAME)
	}else if(input.format == "10x" | input.format == "10X"){
		h5.infile <- hdf5r::H5File$new(filename = file.path(input.file,'filtered_feature_bc_matrix.h5'), mode = "r")
		GENENAME = data.frame(gid = h5.infile[["matrix/features/id"]][], 
			gname = h5.infile[["matrix/features/name"]][], stringsAsFactors=FALSE)
		GENENAME = GENENAME[order(GENENAME$gid),]
		GENENAME$uname = make.unique(GENENAME$gname)
		return (GENENAME)
	}else if(input.format == "csvdir"){
		feature.gzfile = file.path(input.file, "features.tsv.gz")
		GENENAME = Read_genename(feature.gzfile)
		return (GENENAME)
	}else{
		return (NULL)
	}
}



PercentageTopN <- function(seurat_obj, n = 20){
	#count_data = as(seurat_obj@assays$RNA@counts,"matrix")
	count_data = GetAssayData(seurat_obj)
	perN = foreach(i=1:ncol(count_data), .combine=c) %do% {
		sort_column = sort(count_data[,i],decreasing=T)
		if(n<=0){
			n = 20
		}
		sum(sort_column[1:n])/sum(sort_column)*100
	}
	perN.df = data.frame(topN = perN, row.names = colnames(count_data))
	return(perN.df)
}

determine_threshold <- function(x, prob, min.x=NULL, max.x=NULL, sep=1){
	if(is.null(prob)){
		return (NULL)
	}
	q.x = round(quantile(x, prob))
	if(is.null(min.x)) { min.x = min(x) }
	if(is.null(max.x)) { max.x = max(x) }
	if(q.x > max.x){
		print("percentage > max threshold\n")
	}
	for(i in seq(min.x, max.x, sep)){
		if(q.x < i+sep){
			q.x = i+sep
			break
		}
	}
	return (q.x)
}

determine_PCnum <- function(pca_fitted, pc.min = 10, pc.max = 20){
	pc.max = min(pc.max, nrow(pca_fitted))
	pc = pc.max
	for (i in 2:pc.max){
		pc.diff = pca_fitted[i-1] - pca_fitted[i]
		if(pc.diff < 0.05){
			pc = i
			break
		}
	}
	return(max(pc.min, pc))
}

dual.plot <- function(fig, file.prefix, w=7, h=7, res=75){
	pdf(paste(file.prefix,".pdf",sep=""), width = w, height = h)
	print(fig)
	dev.off()
	png(paste(file.prefix,".png",sep=""), width = w*res, height = h*res, res = res, type="cairo-png")
	print(fig)
	dev.off()
}

dual.plot_list <- function(plot.list, file.prefix, title="", scale=0.9, ncol=NULL, w.each=7, h.each=7, res=75, ...){
	n = length(plot.list)
	w.num = 3
	h.num = ceiling(n/3)
	if(n == 1){
		w.num = 1
		h.num = 1
	}else if(n == 2){
		w.num = 2
		h.num = 1
	}else if(n == 4){
		w.num = 2
		h.num = 2
	}
	if(!is.null(ncol)){
		w.num = ncol
		h.num = floor(n/ncol)
	}
	fig = plot_grid(plotlist = plot.list, ncol = w.num, scale=scale, labels=title, ...)
	w = w.each * w.num
	h = h.each * h.num
	pdf(paste(file.prefix,".pdf",sep=""), width = w, height = h)
	print(fig)
	dev.off()
	png(paste(file.prefix,".png",sep=""), width = w*res, height = h*res, res = res, type="cairo-png")
	print(fig)
	dev.off()
}

write.xls <- function(df, xls.file, first.colname="Gene", sep="\t", quote=FALSE, ...) {
	dat = data.frame(rownames(df), df, check.names=F)
	colnames(dat)[1] = first.colname
	write.table(dat, file = xls.file, sep=sep, quote=quote, row.names=F, ...)
}

write.xls.genename <- function(df, xls.file, GENENAME=NULL, sep="\t", quote=FALSE, ...) {
	if(is.null(GENENAME)){
		write.xls(df, xls.file, first.colname="Gene", sep="\t", quote=FALSE, ...)
	}else{
		dat = data.frame(Gene = map.genename(rownames(df), GENENAME, from=3, to=1),
			GeneName = map.genename(rownames(df), GENENAME, from=3, to=2), df, check.names=F)
		write.table(dat, file = xls.file, sep=sep, quote=quote, row.names=F, ...)
	}
}

## genename colomn names: gene.id, gene.name, unique.name
map.genename <- function(gid, GENENAME=NULL, from=1, to=3){
	if(is.null(GENENAME)){
		return (gid)
	}
	map.vector = GENENAME[,to]
	names(map.vector) = GENENAME[,from]
	x = map.vector[gid]
	indices.na = which(is.na(x))
	x[indices.na] = gid[indices.na]
	names(x) = NULL
	return (x)
}


get_mtgenes <- function(genenames, value=T){
	mt.genes = grep("^MT-", genenames, ignore.case = TRUE, value=value)
	return (mt.genes)
}

find_target_genes <- function(SCrna, genelist, GENENAME=NULL, use.name=FALSE, mt=FALSE){
	all.genes = rownames(SCrna)
	genes_in_scrna = genelist[toupper(genelist) %in% toupper(all.genes)]
	if(length(genes_in_scrna) == 0){
		genes_in_scrna = NULL
	}
	if(use.name){
		if(!is.null(GENENAME)){
			indices.id = which(toupper(GENENAME[,1]) %in% toupper(genelist))
			indices.name = which(toupper(GENENAME[,2]) %in% toupper(genelist))
			indices = indices.name %||% indices.id
			if(mt){
				indices = unique(c(indices, get_mtgenes(GENENAME[,3],value=F)))
			}
			target.genes = GENENAME[indices,3]
			target.genes = unique(c(target.genes, genes_in_scrna))
			target.genes = all.genes[toupper(all.genes) %in% toupper(target.genes)]
			return (target.genes)
		}else{
			target.genes = genes_in_scrna
			if(mt){
				mt_in_scrna = get_mtgenes(all.genes)
				target.genes = unique(c(target.genes, mt_in_scrna))
			}
			target.genes = all.genes[toupper(all.genes) %in% toupper(target.genes)]
			return (target.genes)
		}
	}else{
		if(!is.null(GENENAME)){
			indices.id = which(toupper(GENENAME[,1]) %in% toupper(genelist))
			indices.name = which(toupper(GENENAME[,2]) %in% toupper(genelist))
			indices = indices.id %||% indices.name
			if(mt){
				indices = unique(c(indices, get_mtgenes(GENENAME[,3],value=F)))
			}
			target.genes = GENENAME[indices,1]
			target.genes = unique(c(target.genes, genes_in_scrna))
			target.genes = all.genes[toupper(all.genes) %in% toupper(target.genes)]
			return (target.genes)
		}else{
			target.genes = genes_in_scrna
			if(mt){
				mt_in_scrna = get_mtgenes(all.genes)
				target.genes = unique(c(target.genes, mt_in_scrna))
			}
			target.genes = all.genes[toupper(all.genes) %in% toupper(target.genes)]
			return (target.genes)
		}
	}
	return (NULL)
}

filterPlot <- function(SCrna, feature, low.threshold = 0, high.threshold = Inf, label=NULL, log_y = FALSE, coord_flip = TRUE, assay="Spatial", ...){
	value = SCrna@meta.data[,feature]
    value = sort(value)
    total_num = length(value)
    filter_num = sum(value > low.threshold & value < high.threshold)
    filter_percent = round(filter_num/total_num*100,2)
	if(assay == "Spatial"){
	    filter_text = paste("Total spots: ", total_num, "\nFiltered spots: ", filter_num,"\nPercent: ", filter_percent," %",sep="")
	}else{
	    filter_text = paste("Total cells: ", total_num, "\nFiltered cells: ", filter_num,"\nPercent: ", filter_percent," %",sep="")
	}
    x_text = length(value)*0.3
    y_text = max(value) 

    y.label = feature
    if(!is.null(label)){
        y.label = label
    }
    d = data.frame(x = 1:length(value), y = value, filter=factor(as.numeric(value > low.threshold & value < high.threshold),levels=c(0,1)))
    p = ggplot(data = d, aes(x = x, y =y, colour=filter)) +
        geom_point(show.legend=FALSE) +
        scale_colour_manual(values=c("0"="#B6B6B6","1"="#142346"))+
        ylab(y.label) +
        theme(plot.title = element_text(hjust = 0.5)) +
        annotate("text", label = filter_text, x = x_text, y = y_text, adj = 1.1, size = 3.5)
	if(assay == "Spatial"){
		p = p + xlab("Spot number") + ggtitle(paste("Filtering spots with ", y.label, sep=""))
	}else{
		p = p + xlab("Cell number") + ggtitle(paste("Filtering cells with ", y.label, sep=""))
	}
    if(!is.null(low.threshold) & !is.infinite(low.threshold) & low.threshold>0){
        y.low = low.threshold + max(value) * 0.05
        p = p + geom_hline(yintercept = low.threshold, col="#509628") +
            annotate("text", label = paste("low = ", low.threshold, sep=""), x = max(d$x), y = low.threshold, adj = -0.1, size=4)
    }
    if(!is.null(high.threshold) & !is.infinite(high.threshold) & high.threshold>0){
        y.high = high.threshold - max(value) * 0.05
        p = p + geom_hline(yintercept = high.threshold, col="#509628") +
            annotate("text", label = paste("high = ", high.threshold, sep=""), x = 0, y = high.threshold, adj = 1.1, size=4)
    }
    if(log_y){
        p = p + scale_y_log10()
    }
    if(coord_flip){
        p = p + coord_flip()
    }
    return (p)
}



`%||%` <- function(lhs, rhs) {
  if (!is.null(x = lhs) & length(lhs) > 0) {
    return(lhs)
  } else {
    return(rhs)
  }
}

`%iff%` <- function(lhs, rhs) {
  if (!is.null(x = lhs) & length(lhs) > 0) {
    return(rhs)
  } else {
    return(lhs)
  }
}

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
geom_spatial <-  function(mapping = NULL,
                         data = NULL,
                         stat = "identity",
                         position = "identity",
                         na.rm = FALSE,
                         show.legend = NA,
                         inherit.aes = FALSE,
                         ...) {
  
  GeomCustom <- ggproto(
    "GeomCustom",
    Geom,
    setup_data = function(self, data, params) {
      data <- ggproto_parent(Geom, self)$setup_data(data, params)
      data
    },
    
    draw_group = function(data, panel_scales, coord) {
      vp <- grid::viewport(x=data$x, y=data$y)
      g <- grid::editGrob(data$grob[[1]], vp=vp)
      ggplot2:::ggname("geom_spatial", g)
    },
    
    required_aes = c("grob","x","y")
    
  )
  
  layer(
    geom = GeomCustom,
    mapping = mapping,
    data = data,
    stat = stat,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}


## Plotting

plot_spatial_metric <- function(image.basic, 
                            stat.vector, 
                            image.tibble, 
                            sample.name, 
                            legend.text,
                            point.size=1.75,
                            color.palette=myPalette)
{
   stat.matrix = data.frame(barcode = names(stat.vector), value = stat.vector)
   data.plot = merge(image.basic, stat.matrix)
   p = data.plot %>% 
		ggplot(aes(x=imagecol,y=imagerow,fill=value))
   if(!is.null(image.tibble)){
     p = p + geom_spatial(data=image.tibble, aes(grob=grob), x=0.5, y=0.5, crop=TRUE)
   }
   p = p + geom_point(shape = 21, colour = "black", size = point.size, stroke = 0)+
		coord_cartesian(expand=FALSE)+ 
		scale_fill_gradientn(colours = color.palette(100))+
		xlim(0, max(data.plot %>% select(width)))+
		ylim(max(data.plot %>% select(height)),0)+
		xlab("") + ylab("") + labs(fill = legend.text) + ggtitle(sample.name) + 
		theme_set(theme_bw(base_size = 10))+
		theme(plot.title = element_text(hjust = 0.5),
		   panel.grid.major = element_blank(),
		   panel.grid.minor = element_blank(),
		   panel.background = element_blank(), 
		   axis.line = element_line(colour = "black"),
		   axis.text = element_blank(),
		   axis.ticks = element_blank())
  return(p)
}		


plot_mosaic_metric <- function(image.basic, 
                            stat.vector, 
                            sample.name, 
                            legend.text, 
                            plot.func="tile", 
                            lwd.tile = 0.5, 
                            size.point = 1.75, 
                            color.palette=myPalette)
{
   stat.matrix = data.frame(barcode = names(stat.vector), value = stat.vector)
   data.plot = merge(image.basic, stat.matrix)
   data.plot$mosaicrow = data.plot$row
   if(plot.func=="tile"){
	   data.plot$mosaiccol = data.plot$col - data.plot$row %% 2
   }else{
	   data.plot$mosaiccol = data.plot$col
   }
   p = data.plot %>% 
    ggplot(aes(x=mosaiccol,y=mosaicrow,fill=value)) +
    geom_tile(aes(fill = value),colour = "white",lwd = lwd.tile) +
    scale_fill_gradientn(colours = color.palette(100))+
	xlim(0, max(data.plot %>% select(mosaiccol)))+
	ylim(max(data.plot %>% select(mosaicrow)),0)+
    xlab("") + ylab("") + labs(fill = legend.text) + ggtitle(sample.name) + 
    theme_set(theme_bw(base_size = 10))+
    theme(plot.title = element_text(hjust = 0.5),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       panel.background = element_blank(), 
       axis.line = element_line(colour = "black"),
       axis.text = element_blank(),
       axis.ticks = element_blank())
   if(plot.func=="tile"){
      p = p + geom_tile(aes(fill = value),colour = "white",lwd = lwd.tile)
   }else{
      p = p + geom_point(shape = 21, colour="white", size = size.point, stroke = 0.1)
   }
  return(p)
}   



plot_spatial_cluster <- function(image.basic, 
                            cluster.vector, 
                            image.tibble, 
                            sample.name, 
                            legend.text, 
                            point.size=1.75,
                            colours = NULL)
{
   cluster.matrix = data.frame(barcode = names(cluster.vector), cluster = cluster.vector)
   data.plot = merge(image.basic, cluster.matrix)
   clusters = sort(unique(cluster.vector))
   p = data.plot %>% 
		ggplot(aes(x=imagecol, y=imagerow, fill=factor(cluster)))
   if(!is.null(image.tibble)){
     p = p + geom_spatial(data=image.tibble, aes(grob=grob), x=0.5, y=0.5)
   }  
   p = p + geom_point(shape = 21, colour = "black", size = point.size, stroke = 0.2)+
		coord_cartesian(expand=FALSE)+ 
		xlim(0, max(data.plot %>% select(width)))+
		ylim(max(data.plot %>% select(height)),0)+
		xlab("") + ylab("") + labs(fill = legend.text) + ggtitle(sample.name) + 
		guides(fill = guide_legend(override.aes = list(size=3)))+
		theme_set(theme_bw(base_size = 10))+
		theme(plot.title = element_text(hjust = 0.5),
		   panel.grid.major = element_blank(),
		   panel.grid.minor = element_blank(),
		   panel.background = element_blank(), 
		   axis.line = element_line(colour = "black"),
		   axis.text = element_blank(),
		   axis.ticks = element_blank())
   if(!is.null(colours)){
     p = p + scale_fill_manual(values = colours)
   }
   return(p)
}		

plot_mosaic_cluster <- function(image.basic, 
                            cluster.vector, 
                            sample.name, 
                            legend.text, 
                            plot.func="tile",
                            lwd.tile=0.5,
                            size.point=1.75,
                            colours = NULL)
{
   cluster.matrix = data.frame(barcode = names(cluster.vector), cluster = cluster.vector)
   data.plot = merge(bc, cluster.matrix)
   data.plot$mosaicrow = data.plot$row
   if(plot.func=="tile"){
	   data.plot$mosaiccol = data.plot$col - data.plot$row %% 2
   }else{
	   data.plot$mosaiccol = data.plot$col
   }
   clusters = sort(unique(cluster.vector))
   p = data.plot %>% 
    ggplot(aes(x=mosaiccol, y=mosaicrow, fill=factor(cluster))) +
	xlim(0, max(data.plot %>% select(mosaiccol)))+
	ylim(max(data.plot %>% select(mosaicrow)),0)+
    xlab("") + ylab("") + labs(fill = legend.text) + ggtitle(sample.name) + 
    guides(fill = guide_legend(override.aes = list(size=3)))+
    theme_set(theme_bw(base_size = 10))+
    theme(plot.title = element_text(hjust = 0.5),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       panel.background = element_blank(), 
       axis.line = element_line(colour = "black"),
       axis.text = element_blank(),
       axis.ticks = element_blank())
   if(plot.func=="tile"){
      p = p + geom_tile(aes(fill = cluster),colour = "white",lwd = lwd.tile)
   }else{
      p = p + geom_point(shape = 21, colour="white", size = size.point, stroke = 0.1)
   }
  if(!is.null(colours)){
    p = p + scale_fill_manual(values = colours)
  }
  return(p)
}   

jamPDF = function(in.files,
                  out.file = 'jamPDF.pdf',
                  layout = '3x3',
                  delete.original = TRUE,
                  crop.white = TRUE,
                  page = 'none',
                  hide.output = TRUE,
                  toJam = FALSE,
                  ignore.stderr = TRUE)
{
  # Global option to control for this untill I remove it.
  if(!toJam)
  {
    message("jamPDF disabled")
    return(invisible(1))
  }
  ######################################################
  
  in.files = in.files[sapply(in.files, file.exists)]

  if(length(in.files) == 0) stop("All the input files that you asked to jam are missing -- check input!")

  pio::pioHdr(paste("REVOLVER jamPDF to", out.file),
              toPrint = c(
                `Input files` = paste(in.files, collapse = ', '),
                `PDF layout` = layout,
                `PDF page type` = page,
                `PDF crop white margins` = crop.white,
                `PDF delete input files` = delete.original),
              prefix = '\t -'
              )

  stopifnot(!is.na(out.file))

  cmd = paste(
    'pdfjam ',
    paste(in.files, collapse = ' '),
    ' --nup ',
    layout,
    ' --landscape --outfile ',
    out.file,
    sep = ' '
  )


  # pio::pioTit("Assembling PDFs")
  aa = system(cmd, intern = hide.output, ignore.stderr = ignore.stderr)

  if (crop.white)
  {
    # pio::pioTit("Cropping white margins (3 3 3 3)")

    aa = system(
      paste('pdfcrop --margins "3 3 3 3"', out.file, out.file),
      intern = hide.output,
      ignore.stderr = ignore.stderr
    )
  }


  if (page != 'none')
  {
    page = R.utils::capitalize(page)

    # pio::pioTit(paste("Resizing pages to", page))


    installation = find.package('revolver')
    cmd = paste(paste0(installation, '/bin/pdfScale.sh'),
                '-v -r',
                page,
                out.file)


    aa = system(cmd, intern = hide.output, ignore.stderr = ignore.stderr)


    f = gsub(out.file, pattern = '.pdf', replacement = '')
    file.rename(paste(f, page, 'pdf', sep = '.'), paste(f, 'pdf', sep =
                                                          '.'))
  }

  if (delete.original)
    file.remove(setdiff(in.files, out.file))

  invisible(NULL)
}

FilterSeuratObject = function(seurat.obj, 
	assay = "RNA",
	low.feature = NULL, 
	high.feature = NULL, 
	mt.genes = NULL, 
	high.mt = NULL, 
	hb.genes = NULL,
	high.hb = NULL,
	plot.file = NULL
){
	variable.name = paste("nFeature_",assay,sep="")
	if(is.null(low.feature)){
		low.feature = 0
	}
	if(low.feature>0 & low.feature<1){
		low.feature = determine_threshold(seurat.obj@meta.data[,variable.name], low.feature, 500, 25000, 500)
	}
	if(is.null(high.feature)){
		high.feature = Inf
	}
	if(high.feature>0 & high.feature<1){
		high.feature = determine_threshold(seurat.obj@meta.data[,variable.name], high.feature, 500, 25000, 500)
	}
	cell.list = list()
	plot.list = list()
	freq = 1
	## feature/gene number filtration
	cell.list[[1]] = rownames(seurat.obj@meta.data)[which(seurat.obj@meta.data[,variable.name]>low.feature & 
		seurat.obj@meta.data[,variable.name]<high.feature)]
	plot.list[[freq]] = filterPlot(seurat.obj, variable.name, low.threshold = low.feature, 
		high.threshold = high.feature, label = "nFeature", coord_flip=TRUE, assay=assay)
#	if(assay == "RNA"){
#		cell.list[[1]] = WhichCells(seurat.obj,expression = nFeature_RNA>low.feature & nFeature_RNA<high.feature)
#	}else{
#		cell.list[[1]] = WhichCells(seurat.obj,expression = nFeature_Spatial>low.feature & nFeature_Spatial<high.feature)
#	}
	## MT genes filtration
	if(length(mt.genes)>0){
		#seurat.obj[['percent_MT']] = PercentageFeatureSet(seurat.obj, features = mt.genes)
		if(!is.null(high.mt)){
			if(high.mt>0 & high.mt<1){
				high.mt = determine_threshold(seurat.obj@meta.data$percent_MT, high.mt, 0, 30, 1)
			}
#			cell.list[[2]] = WhichCells(seurat.obj, expression = percent_MT < high.mt)
			cell.list[[2]] = rownames(seurat.obj@meta.data)[which(seurat.obj@meta.data[,"percent_MT"] < high.mt)]
			freq = freq + 1
			plot.list[[freq]] = filterPlot(seurat.obj, "percent_MT", low.threshold = 0, 
				high.threshold = high.mt, label = "percent_MT (%)", coord_flip=TRUE, assay=assay)
		}
	}
	## HB genes filtration
	if(length(hb.genes)>0){
		#seurat.obj[['percent_HB']] = PercentageFeatureSet(seurat.obj, features = hb.genes)
		if(!is.null(high.hb)){
			if(high.hb>0 & high.hb<1){
				high.hb = determine_threshold(seurat.obj@meta.data$percent_HB, high.hb, 0, 30, 1)
			}
			#cell.list[[2]] = WhichCells(seurat.obj, expression = percent_HB < high.hb)
			cell.list[[3]] = rownames(seurat.obj@meta.data)[which(seurat.obj@meta.data[,"percent_HB"] < high.hb)]
			freq = freq + 1
			plot.list[[freq]] = filterPlot(seurat.obj, "percent_HB", low.threshold = 0, 
				high.threshold = high.hb, label = "percent_HB (%)", coord_flip=TRUE, assay=assay)
		}
	}
	select.cells = names(which(table(unlist(cell.list))==freq))
	filter.genes = c(mt.genes, hb.genes)
	cat(filter.genes)
	select.genes = rownames(seurat.obj)[!rownames(seurat.obj) %in% filter.genes]
	print(length(select.cells))
	print(length(select.genes))
	cat(length(plot.list),"\n")
	if(!is.null(plot.file)){
		p = CombinePlots(plot.list, ncol=length(plot.list), scale=0.9)
		dual.plot(p, plot.file, w=length(plot.list)*6, h=6, res=200)
	}
	#return (list(cells = select.cells, features = select.genes))
	seurat.obj = seurat.obj[select.genes, select.cells]  ## SubsetData(seurat.obj, cells = select.cells)
	#seurat.obj = subset(seurat.obj, cells = select.cells, features = select.genes)
	return(seurat.obj)
}

plot_SpatialCluster <- function(SCrna,
                            sample.name, 
                            legend.text, 
                            image.name="image", 
                            cluster.label="seurat_clusters", 
                            point.size=1.75,
                            stroke = 0.2,
                            alpha = 1,
                            legend = TRUE,
							width = NULL,
                            plot.image = TRUE,
                            colours = NULL)
{
   active.image = SCrna[[image.name]]
   height = nrow(active.image@image)
   width = ncol(active.image@image)
   image.grob = rasterGrob(active.image@image, width=unit(1,"npc"), height=unit(1,"npc"))
   image.tibble= tibble(sample=factor(sample.name), grob=list(image.grob))
   image.tibble$height = height
   image.tibble$width = width

   image.basic = active.image@coordinates
   image.basic$barcode = rownames(image.basic)
   image.basic$sample = sample.name
   image.basic$imagerow = image.basic$imagerow * active.image@scale.factors$lowres
   image.basic$imagecol = image.basic$imagecol * active.image@scale.factors$lowres
   image.basic$height = height
   image.basic$width = width

   cluster.matrix = data.frame(barcode = rownames(SCrna@meta.data), cluster = SCrna@meta.data[,cluster.label])
   data.plot = merge(image.basic, cluster.matrix)
   p = data.plot %>% 
		ggplot(aes(x=imagecol, y=imagerow, color=factor(cluster)))
   if(!is.null(image.tibble) & plot.image){
     p = p + geom_spatial(data=image.tibble, aes(grob=grob), x=0.5, y=0.5)
   }  
   #p = p + geom_point(shape = 21, colour = "#FFFFFF", alpha=alpha, size = point.size, stroke = stroke)+
   p = p + geom_point(shape = 19, alpha=alpha, size = point.size, stroke = stroke)+
		coord_cartesian(expand=FALSE)+ 
		ylim(max(data.plot %>% select(height)),0)+
		xlab("") + ylab("") + labs(color = legend.text) + ggtitle(sample.name) + 
		guides(color = guide_legend(override.aes = list(size=3)))+
		theme_set(theme_bw(base_size = 10))+
		theme(plot.title = element_text(hjust = 0.5,size=20),
		   panel.border = element_blank(),
		   panel.grid.major = element_blank(),
		   panel.grid.minor = element_blank(),
		   panel.background = element_blank(), 
		   axis.line = element_blank(),
		   axis.text = element_blank(),
		   axis.ticks = element_blank())
   if(!is.null(width)){
     p = p + xlim(0, width)
   }else{
     p = p + xlim(0, max(data.plot %>% select(width)))
   }
   if(!is.null(colours)){
     p = p + scale_color_manual(values = colours)
   }
   if(!legend){
     p = p + guides(color=FALSE)
   }
   return(p)
}		


clustering_barplot = function(seurat.obj, cluster.cols=NULL, group.by="orig.ident", title="Clustering", flip=FALSE, item=NULL){
    if(is.null(item)){
        cluster.matrix = data.frame(cluster=Idents(seurat.obj), sample=seurat.obj[[group.by]][,1])
    }else{
        cluster.matrix = data.frame(cluster=seurat.obj[[item]][,1], sample=seurat.obj[[group.by]][,1])
    }
    clust_stat = data.frame(table(cluster.matrix))
    cluster.bar = ggplot(data=clust_stat,aes(sample,weight=Freq,fill=cluster))+
        geom_bar(position="fill", width=0.7) +
        ggtitle(title) +
        ylab("Relative proportion")
        theme(plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank())
    if(!is.null(cluster.cols)){
        cluster.bar = cluster.bar + scale_fill_manual(values=cluster.cols)
    }
    if(flip){
        cluster.bar = cluster.bar + coord_flip()
    }
    return(cluster.bar)
}

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
plot_SpatialFeature <- function(SCrna,
                            sample.name, 
                            feature, 
                            image.name="image", 
                            point.size=1.75,
                            stroke = 0.2,
                            alpha = 1,
                            legend = TRUE,
							width = NULL,
                            plot.image = TRUE,
                            color.palette=myPalette)
{
   active.image = SCrna[[image.name]]
   height = nrow(active.image@image)
   width = ncol(active.image@image)
   image.grob = rasterGrob(active.image@image, width=unit(1,"npc"), height=unit(1,"npc"))
   image.tibble= tibble(sample=factor(sample.name), grob=list(image.grob))
   image.tibble$height = height
   image.tibble$width = width

   image.basic = active.image@coordinates
   image.basic$barcode = rownames(image.basic)
   image.basic$sample = sample.name
   image.basic$imagerow = image.basic$imagerow * active.image@scale.factors$lowres
   image.basic$imagecol = image.basic$imagecol * active.image@scale.factors$lowres
   image.basic$height = height
   image.basic$width = width

   expr = data.frame(value=GetAssayData(SCrna,slot="data")[feature,])
   expr$barcode = rownames(expr)
   data.plot = merge(image.basic, expr)
   p = data.plot %>% 
		ggplot(aes(x=imagecol, y=imagerow, color=value))
   if(!is.null(image.tibble) & plot.image){
     p = p + geom_spatial(data=image.tibble, aes(grob=grob), x=0.5, y=0.5)
   }  
   #p = p + geom_point(shape = 21, colour = "#FFFFFF", alpha=alpha, size = point.size, stroke = stroke)+
   p = p + geom_point(shape = 19, alpha=alpha, size = point.size, stroke = stroke)+
		coord_cartesian(expand=FALSE)+ 
		scale_color_gradientn(colours = color.palette(100))+
		ylim(max(data.plot %>% select(height)),0)+
		xlim(0, max(data.plot %>% select(width))) + 
		xlab("") + ylab("") + ggtitle(sample.name) + 
		theme_set(theme_bw(base_size = 10))+
		theme(plot.title = element_text(hjust = 0.5,size=20),
		   panel.border = element_blank(),
		   panel.grid.major = element_blank(),
		   panel.grid.minor = element_blank(),
		   panel.background = element_blank(),
		   legend.position = "top", 
		   axis.line = element_blank(),
		   axis.text = element_blank(),
		   axis.ticks = element_blank())
   if(!legend){
     p = p + guides(color=FALSE)
   }
   return(p)
}		


