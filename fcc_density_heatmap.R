# Rscript fcc_density_heatmap.R

# each column => dataset; color-coded by density

options(scipen = 100)

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")

script_name <- "fcc_density_heatmap.R"
cat("> START ", script_name, "\n")
startTime <- Sys.time()

script0_name <- "0_prepGeneData"

require(foreach)
require(doMC)
require(ggpubr)
require(ggplot2)
registerDoMC(40)

plotType <- "png"
myHeight <- 400
myWidth <- 400
myHeightGG <- 7
myWidthGG <- 12

plotCex <- 1.4

source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("../FIGURES_V2_YUANLONG/settings.R")

pipFolder<- file.path("../v2_Yuanlong_Cancer_HiC_data_TAD_DA/")
stopifnot(dir.exists(pipFolder))

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

pipOutFolder <- file.path(pipFolder, "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipOutFolder))

outFolder <- "FCC_DENSITY_HEATMAP" 
dir.create(outFolder, recursive = TRUE)

all_hicds <- list.files(pipOutFolder)

all_hicds <- all_hicds[!(grepl("RANDOM", all_hicds) | grepl("PERMUT", all_hicds))]
# all_hicds <- all_hicds[! (grepl("RANDOM", all_hicds) | grepl("PERMUT", all_hicds)) ]
# all_hicds <- all_hicds[grepl("ENCSR489OCU_NCI-H460_40kb", all_hicds)]
# all_hicds = "ENCSR489OCU_NCI-H460_40kb"
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipOutFolder, x)))


rd_patterns <- c("RANDOMMIDPOS", "RANDOMNBRGENES", "RANDOMMIDPOSDISC" , "RANDOMMIDPOSSTRICT", "RANDOMSHIFT", "PERMUTG2T")


# retrieve dataset order
tmp <- get(load("../FIGURES_V3_YUANLONG/BARPLOT_FCC_AUC_RATIO/all_dt.Rdata"))
tmp <- tmp[order(tmp$fcc_auc, decreasing = TRUE),]
ds_levels <- file.path(tmp$hicds, tmp$exprds)

buildData <- TRUE

if(buildData){
  hicds = all_hicds[1]
  hicds = all_hicds[2]
  all_result_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
    exprds = all_exprds[[paste0(hicds)]][1]
    hicds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      
      cat("... start ", hicds, " - ", exprds, " \n")
      
      fcc_file <- file.path(pipOutFolder, hicds, exprds, "8cOnlyFCC_runAllDown", "all_obs_prodSignedRatio.Rdata")
      
      if(!file.exists(fcc_file)) return(NULL)
      stopifnot(file.exists(fcc_file))
      all_fcc <- as.numeric(get(load(fcc_file)))
      # if(!file.exists(fcc_file)) {
      #   fcc_file <- file.path(pipOutFolder, hicds, exprds, "8cOnlyFCConlyObs_runAllDown", "all_obs_prodSignedRatio.Rdata")  
      # }
      
      data.frame(
        hicds = hicds,
        exprds=exprds,
        FCC = all_fcc,
        stringsAsFactors = FALSE
      )
    }
    hicds_dt
  }
  outFile <- file.path(outFolder, "all_result_dt.Rdata")
  save(all_result_dt, file = outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  
} else {
  outFile <- file.path(outFolder, "all_result_dt.Rdata")
  all_result_dt <- get(load(outFile))
}  

all_result_dt$dataset <- file.path(all_result_dt$hicds, all_result_dt$exprds)
all_result_dt$dataset <- factor(all_result_dt$dataset, levels=ds_levels)
stopifnot(!is.na(all_result_dt$dataset))
all_result_dt$heatmap_x <- as.numeric(all_result_dt$dataset) 

all_result_dt$ds_col <- all_cols[all_cmps[basename(all_result_dt$exprds)]]

all_result_dt$hicds_lab <- gsub(".+RANDOMSHIFT_40kb", "RANDOMSHIFT", 
                                gsub(".+RANDOMMIDPOSDISC_40kb", "RANDOMMIDPOSDISC",
                                     gsub(".+RANDOMMIDPOSSTRICT_40kb", "RANDOMMIDPOSSTRICT",
                                          gsub(".+RANDOMNBRGENES_40kb", "RANDOMNBRGENES",
                                               gsub(".+PERMUTG2T_40kb", "PERMUTG2T", 
                                                    gsub(".+RANDOMMIDPOS_40kb", "RANDOMMIDPOS", all_result_dt$hicds))))))
all_result_dt$hicds_lab[! all_result_dt$hicds_lab %in% rd_patterns] <- "OBSERVED"

all_types <- unique(all_result_dt$hicds_lab)

for(a_t in all_types) {

  plot_dt <- all_result_dt[all_result_dt$hicds_lab == a_t,]
  
  nDS <- length(unique(as.character(plot_dt$dataset)))
  
  density_plot <- ggplot(plot_dt, aes(x = heatmap_x, y = FCC))+ 
    stat_density2d(geom="tile",   aes(fill = ..density..), contour = FALSE) +
    theme(
      plot.background = element_blank()
    ) + 
    ggtitle(paste0("FCC score distribution -", a_t),   
            subtitle = paste0("all datasets (n=", nDS, ")"))+
    scale_x_discrete(name="Datasets ranked by decreasing AUC FCC ratio", labels = rep(labsymbol, nDS ))  + 
    scale_y_discrete(name="FCC score")  + 
    labs(fill = "Density")+
    scale_fill_gradient( high="red", low="blue", na.value = "white" )  +
    theme(
      axis.text.x = element_text(colour = plot_dt$ds_col, size=12),
      axis.text.y= element_text(colour = "black", size=12),
      axis.title.x = element_text(colour = "black", size=14, face="bold"),
      axis.title.y = element_text(colour = "black", size=14, face="bold"),
      plot.title = element_text(hjust=0.5, size=16, face="bold"),
      plot.subtitle = element_text(hjust=0.5, size=14, face="italic"),
      panel.background = element_rect(fill = "transparent")
      # legend.background =  element_rect()
    )

  outFile <- file.path(outFolder, paste0("FCC_score_dist_allDS_", a_t, ".", plotType))
  ggsave(density_plot, filename = outFile,  height=myHeightGG, width=myWidthGG)
  cat(paste0("... written: ", outFile, "\n"))
  
}



# 
# v+geom_raster(aes(fill = Value),interpolate = T)+scale_fill_gradientn(colours = terrain.colors(10))                                     
# 
# ggplot(plot_dt, aes(x = dataset, y = density_x, fill = density_y, color=density_y))+ 
#   geom_tile()
# 
# 
# plot_dt <- do.call(rbind, by(all_result_dt, all_result_dt$dataset, function(x) {
#   hicds <- unique(x$hicds)
#   stopifnot(length(hicds) == 1)
#   exprds <- unique(x$exprds)
#   stopifnot(length(exprds) == 1)
#   ds_dt <- density(x$FCC)
#   data.frame(
#     hicds=hicds,
#     exprds=exprds,
#     density_x  = ds_dt$x,
#     density_y  = ds_dt$y,
#     stringsAsFactors = FALSE
#   )
# }))
# 
# plot_dt$dataset <- file.path(plot_dt$hicds, plot_dt$exprds)
# plot_dt$dataset <- factor(plot_dt$dataset, levels=ds_levels)
# stopifnot(!is.na(plot_dt$dataset))
# 
# 
# plot_dt$heatmap_x <- as.numeric(plot_dt$dataset) 
# heatmap_x <- as.numeric(plot_dt$dataset)
# heatmap_y <- plot_dt$density_x
# 
# library(akima)
# library(tidyverse)
# resolution <- 0.1 # you can increase the resolution by decreasing this number (warning: the resulting dataframe size increase very quickly)
# ip_plot_dt <- interp(x=heatmap_x, y=heatmap_y, z=plot_dt$density_y, 
#             xo=seq(min(heatmap_x),max(heatmap_x),by=resolution), 
#             yo=seq(min(heatmap_y),max(heatmap_y),by=resolution), duplicate="mean")
# 
# 
# res <- ip_plot_dt$z %>% 
#   magrittr::set_colnames(ip_plot_dt$y) %>% 
#   as_tibble() %>% 
#   mutate(x=ip_plot_dt$x) %>% 
#   gather(y, z, -x, convert=TRUE)
# 
# ggplot(res, aes(x = x, y = y, fill = z, color=z))+ 
#   geom_tile()
# 
# 
# stop("--ok\n")
# plot_pattern <- "" # then for loop -> heatmap for RANDOMMIDPOS RANDOMMIDPOSDISC
# 
# plot_dt <- plot_dt[order(plot_dt$dataset, plot_dt$density_x),]
# filled.contour(x=as.numeric(plot_dt$dataset), y = plot_dt$density_x, z = plot_dt$density_y)
# 
# myData <- mapply(rnorm, 1000, 200, mean=seq(-50,50,0.5))
# 
# myDensities <- apply(myData, 2, density, from=-500, to=500)
# 
# Ys <- sapply(myDensities, "[", "y")
# 
# img <- do.call(cbind, Ys)
# 
# filled.contour(x=1:ncol(img), y=myDensities[[1]]$x, t(img))
# 
# 
# library(lattice)
# levelplot( density_y ~  heatmap_x+density_x, data=plot_dt)
# 
# 
# image(xtabs(density_y~heatmap_x+density_x, plot_dt))
# 
# 
# 
# 
# 
# 
# all_result_dt$hicds_lab <- gsub(".+RANDOMSHIFT_40kb", "RANDOMSHIFT", 
#                                 gsub(".+RANDOMMIDPOSDISC_40kb", "RANDOMMIDPOSDISC",
#                                      gsub(".+RANDOMMIDPOSSTRICT_40kb", "RANDOMMIDPOSSTRICT",
#                                           gsub(".+RANDOMNBRGENES_40kb", "RANDOMNBRGENES",
#                                                gsub(".+PERMUTG2T_40kb", "PERMUTG2T", 
#                                                     gsub(".+RANDOMMIDPOS_40kb", "RANDOMMIDPOS", all_result_dt$hicds))))))
# all_result_dt$hicds_lab[! all_result_dt$hicds_lab %in% rd_patterns] <- "OBSERVED"
# 
# all_types <- unique(all_result_dt$hicds_lab)
# 
# for(a_t in all_types) {
#   
#   plot_dt <- all_result_dt[all_result_dt$hicds_lab == a_t,]
#   stopifnot(nrow(plot_dt) > 0)
#   
#   nDS <- length(unique(file.path(plot_dt$hicds, plot_dt$exprds)))
#   
#   my_x <- plot_dt$adjCombPval
#   my_y <- plot_dt$FCC
#   
#   x_lab <- "TAD adj. comb. p-val [log10]"
#   y_lab <- "TAD FCC"
#   
#   outFile <- file.path(outFolder, paste0("FCC_vs_adjCombPval_", a_t, ".", plotType))
#   do.call(plotType, list(outFile, height=myHeight, width=myWidth))
#   par(bty="L")
#   densplot(
#     x=my_x,
#     y=my_y,
#     xlab=x_lab,
#     ylab=y_lab,
#     main = paste0("FCC vs. p-val - ", a_t),
#     cex.main=plotCex,
#     cex.lab = plotCex,
#     cex.axis = plotCex,
#     cex=0.7
#     
#   )
#   mtext(side=3, text = paste0("all DS (n=", nDS, ")"))
#   addCorr(x=my_x, y=my_y, bty="n")
#   
# }
# 


