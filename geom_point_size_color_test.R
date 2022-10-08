library(dplyr)
library(ggplot2)
library(cowplot)
library(grid)
library(readbitmap)
library(Seurat)
library(dplyr)
library(data.table)
expr=readRDS('combin.data.RDS')
expr[['proportation']]=runif(dim(expr)[2])  ####测试数据，可以是其他的信息
outs=list()
for(slice in names(expr@images)){
  print(slice)
  expr_sub=subset(expr,orig.ident==gsub('\\.','-',slice))
  grobs=list()
  grobs[[1]] <- rasterGrob(expr_sub@images[[slice]]@image, width=unit(1,"npc"), height=unit(1,"npc"))
  images_tibble <- tibble(sample=factor(slice), grob=grobs)
  images_tibble$height <- nrow(expr_sub@images[[slice]]@image)
  images_tibble$width <- ncol(expr_sub@images[[slice]]@image)
  spatial_coord <- data.frame(expr_sub@images[[slice]]@coordinates) %>%
    tibble::rownames_to_column("barcodeID") %>%
    dplyr::mutate(imagerow_scaled = imagerow * expr_sub@images[[slice]]@scale.factors$lowres,
                  imagecol_scaled = imagecol * expr_sub@images[[slice]]@scale.factors$lowres)
  spatial_coord$proportation=expr_sub$proportation  
  print(head(spatial_coord))
  p1=ggplot(data=spatial_coord,mapping=aes(x=imagecol_scaled,y=imagerow_scaled)) +
    geom_spatial(data=images_tibble[1,], aes(grob=grob),image.alpha = 0.3, x=0.5, y=0.5)+
    geom_point(aes(size=proportation,color=proportation), stroke = 0.5)+
    coord_cartesian(expand=FALSE)+
    xlim(0,images_tibble$height)+ylim(images_tibble$width,0)+
    xlab("") +
    ylab("") +
    ggtitle(slice)+
    theme_set(theme_bw(base_size = 10))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text = element_blank(),
          axis.ticks = element_blank())+
    scale_colour_continuous(type='gradient')+
    scale_size_continuous(range=c(0,2))+
    guides(
      color= guide_legend()
      # size=guide_legend(override.aes = list(size = c(0,0.25,0.5,0.75,1)))
      # scale_alpha_continuous
    )+scale_alpha(range=c(0,1))+scale_colour_gradient(low='yellow',high='red')
  outs[[slice]]=p1
}
p2=plot_grid(plotlist = outs,ncol=4)
ggsave('test1.png',p2,width=25,height=6)
ggsave('test1.pdf',p2,width=25,height=6)

geom_spatial <-  function(mapping = NULL,
                          data = NULL,
                          stat = "identity",
                          position = "identity",
                          na.rm = FALSE,
                          show.legend = NA,
                          inherit.aes = FALSE,
                          image.alpha = 0.2,
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
      # img <- editGrob(grob = img.grob, vp = vp)
      g$raster = as.raster(
        x = matrix(
          data = alpha(colour = g$raster, alpha = image.alpha),
          nrow = nrow(x = g$raster),
          ncol = ncol(x = g$raster),
          byrow = TRUE)
      )
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