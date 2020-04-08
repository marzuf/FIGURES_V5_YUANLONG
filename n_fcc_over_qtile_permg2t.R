# Rscript n_fcc_over_qtile_permg2t.R

outFolder <- "N_FCC_OVER_QTILE_PERMG2T"
dir.create(outFolder)


plotType <- "svg"
myHeight <- myWidth <- 7
plotCex <- 1.4

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("../FIGURES_V2_YUANLONG/settings.R")

require(foreach)
require(doMC)
registerDoMC(40)
require(ggplot2)

buildData <- FALSE

probQt <- 0.75
fccThresh <- 0.75

keepPermut <- 1000

plotSub <- paste0("probQt=", probQt, "; fccThresh=", fccThresh, "; keepPermut=", keepPermut)


rd_type <- "PERMG2T"

if(buildData){
  hicds = all_hicds[1]
  hicds = all_hicds[2]
  all_fcc_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
    exprds = all_exprds[[paste0(hicds)]][1]
    hicds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      
      cat("... start ", hicds, " - ", exprds, " \n")
      
      rd_file <- file.path("../v2_Yuanlong_Cancer_HiC_data_TAD_DA/PIPELINE/OUTPUT_FOLDER", hicds, exprds, "8cOnlyFCC_runAllDown/prodSignedRatio_permDT.Rdata")
      rd_fcc_dt <- get(load(rd_file))
      rd_fcc_dt <- rd_fcc_dt[,1:keepPermut]
      rd_fcc <- as.numeric(rd_fcc_dt)
      
      fccFile <- file.path("../v2_Yuanlong_Cancer_HiC_data_TAD_DA/PIPELINE/OUTPUT_FOLDER", hicds, exprds, "8cOnlyFCC_runAllDown/all_obs_prodSignedRatio.Rdata")
      obs_fcc <- get(load(fccFile))
      
     dt1 <- data.frame(
        hicds=hicds,
        exprds=exprds,
        fcc_type="random",
        fcc_value = rd_fcc,
        stringsAsFactors = FALSE
      )
     dt2 <- data.frame(
       hicds=hicds,
       exprds=exprds,
       fcc_type="observed",
       fcc_value = as.numeric(obs_fcc),
       stringsAsFactors = FALSE
     )
     rbind(dt1,dt2)
    }
    hicds_dt
  }
  outFile <- file.path(outFolder, "all_fcc_dt.Rdata")
  save(all_fcc_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFolder, "all_fcc_dt.Rdata")
  all_fcc_dt <- get(load(outFile))
}

all_fcc_dt$dataset <- file.path(all_fcc_dt$hicds, all_fcc_dt$exprds)

ds = unique(all_fcc_dt$dataset)[1]
nOverQt_dt <- foreach(ds = unique(all_fcc_dt$dataset), .combine='rbind') %dopar% {
  curr_dt <- all_fcc_dt[all_fcc_dt$dataset == ds,]
  rd_qt <- quantile(curr_dt$fcc_value[curr_dt$fcc_type == "random"], probs=probQt)
  nObs <- sum(curr_dt$fcc_type == "observed")
  nRd <- sum(curr_dt$fcc_type == "random")
  nOverQt <- sum(curr_dt$fcc_value[curr_dt$fcc_type == "observed"] >= rd_qt )
  nOverThresh_obs <-   sum(curr_dt$fcc_value[curr_dt$fcc_type == "observed"] >= fccThresh )
  nOverThresh_rd <-   sum(curr_dt$fcc_value[curr_dt$fcc_type == "random"] >= fccThresh )
  wilcoxPval <- wilcox.test(curr_dt$fcc_value[curr_dt$fcc_type == "observed"], 
                            curr_dt$fcc_value[curr_dt$fcc_type == "random"], alternative="greater")$p.value
  data.frame(
    hicds=dirname(ds),
    exprds=basename(ds),
    rd_qt = rd_qt,
    ratioOverQt_obs = nOverQt/nObs,
    ratioOverThresh_obs = nOverThresh_obs/nObs,
    ratioOverThresh_rd = nOverThresh_rd/nRd,
    wilcoxPval = wilcoxPval,
    stringsAsFactors = FALSE
  )
}
outFile <- file.path(outFolder, "nOverQt_dt.Rdata")
save(nOverQt_dt, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))


dotcols <- all_cols[all_cmps[nOverQt_dt$exprds]]

for(i in 5:(ncol(nOverQt_dt)-1)) {
  
  xvar <- colnames(nOverQt_dt)[i]
  my_x <- nOverQt_dt[,i]
  
  for(j in (i+1):(ncol(nOverQt_dt))) {
    
    yvar <- colnames(nOverQt_dt)[j]
    my_y <- nOverQt_dt[,j]
    
    if(yvar=="wilcoxPval") my_y <- -log10(my_y)

    outFile <- file.path(outFolder, paste0(yvar, "_vs_", xvar, "_cmp_values_dotplot.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    plot(
      x = my_x,
      y = my_y,
      xlab=xvar,
      ylab=yvar,
      pch=16,
      col = dotcols,
      cex=0.7,
      cex.lab=plotCex,
      cex.axis=plotCex,
      cex.main=plotCex
    )
    if(yvar=="wilcoxPval") abline(h=-log10(0.05), lty=2, col="red")
    mtext(side=3, paste0("all DS (n=", nrow(nOverQt_dt), "); ", plotSub))
    
    
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
  }
  
}

p_fcc_boxplot <-  ggplot(all_fcc_dt, aes(x = fcc_type, y = fcc_value)) + 
  labs(title=paste0("random=", rd_type),y ="FCC score", x="")+
  # geom_boxplot(notch = TRUE, outlier.shape=NA)+
  # geom_jitter(aes(colour = cond), position=position_jitterdodge())+
  geom_boxplot(notch = TRUE, outlier.shape=NA)+
  # geom_point(position=position_jitterdodge(),  alpha=0.5) +
  facet_wrap(~dataset,switch="x") + 
  theme( 
    plot.title = element_text(hjust = 0.5, face = "bold", size=16),
    plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14)
  )

outFile <- file.path(outFolder, paste0("allDS_cmp_FCC_random_boxplot.", plotType))
ggsave(plot = p_fcc_boxplot, filename = outFile, height=myHeightGG*1.5, width = myWidthGG*1.5)
cat(paste0("... written: ", outFile, "\n"))


nOverQt_dt$dataset <- file.path(nOverQt_dt$hicds, nOverQt_dt$exprds)

id_vars <- c("hicds", "exprds", "dataset")

plot_vars <- c("ratioOverThresh", "ratioOverQt")

plotSub <- paste0("probQt=", probQt, "; fccThresh=", fccThresh)

for(p_v in plot_vars) {
  
  p_cols <- colnames(nOverQt_dt)[grepl(p_v, colnames(nOverQt_dt))]
  
  plot_dt <- nOverQt_dt[,c(id_vars, p_cols)]
  require(reshape2)
  plot_dt <- melt(plot_dt, id=id_vars)
  plot_dt$variable <- ifelse(plot_dt$variable == paste0(p_v, "_rd"), "random",
                             ifelse(plot_dt$variable == paste0(p_v, "_obs"), "obs.", NA))
  stopifnot(!is.na(plot_dt$variable))
  
  tmp_dt <- plot_dt[plot_dt$variable == "obs.",]
  tmp_dt <- tmp_dt[order(tmp_dt$value, decreasing=TRUE),]
  ds_levels <- tmp_dt$dataset
  
  labcols <- all_cols[all_cmps[basename(ds_levels)]]
  
  plot_dt$dataset <- factor(plot_dt$dataset, levels=ds_levels)
  stopifnot(!is.na(plot_dt$dataset))
  
  
  p_barplot_ratio <- ggplot(plot_dt, aes(x = dataset, y =  value, fill=variable))+
    labs(title=paste0(rd_type, " - ", p_v), fill="", y = paste0(p_v), subtitle = paste0(plotSub))+
    geom_bar(position="dodge", stat = "identity") +
    scale_x_discrete(labels = rep("°", length(ds_levels)))+
    theme( 
      plot.title = element_text(hjust = 0.5, face = "bold", size=16),
      plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14),
      panel.grid = element_blank(),
      panel.grid.major.y = element_line(colour = "grey"),
      panel.grid.minor.y = element_line(colour = "grey"),
      axis.line.x= element_line(size = .2, color = "black"),
      axis.line.y = element_line(size = .2, color = "black"),
      axis.text.y = element_text(color="black", hjust=1,vjust = 0.5, size=12),
      axis.text.x =element_text(color=labcols, hjust=0.5,vjust = 0.5, size=12, face="bold"),
      # axis.ticks.x = element_blank(),
      axis.title.y = element_text(color="black", size=13),
      # axis.title.x = element_text(color="black", size=13),
      axis.title.x = element_blank(),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "transparent"),
      legend.background =  element_rect(),
      legend.text = element_text(size=12),
      legend.key = element_blank(),
      legend.title = element_text(face="bold")
    )
  
  outFile <- file.path(outFolder, paste0("all_DS_", p_v, "_barplot.", plotType))
  ggsave(plot = p_barplot_ratio, filename = outFile, height=myHeightGG, width = myWidthGG*1.2)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  
}


















