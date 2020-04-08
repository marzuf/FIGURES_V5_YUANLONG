plotCex <- 1.4
plotType <- "svg"

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("../FIGURES_V2_YUANLONG/settings.R")

# Rscript cmp_nFCC_aucRatio.R

outFolder <- "CMP_NFCC_AUCRATIO"
dir.create(outFolder)

####################################################################################################################
####################################################################################################################
####################################################################################################################

nOverQt_dt <- get(load("N_FCC_OVER_QTILE_V2_EITHER/nOverQt_dt.Rdata"))
aucRatio_dt <- get(load("../FIGURES_V3_YUANLONG/RANDOM_FCC_AUC_RATIO_MEANCORRPERMUT_V2_EITHER/rd_RorL_auc_fcc_ratio_dt.Rdata"))

merged_dt <- merge(nOverQt_dt, aucRatio_dt, by=c("hicds", "exprds"))

merged_dt$ratioOverThresh_obsOverRd <- merged_dt$ratioOverThresh_obs/merged_dt$ratioOverThresh_rd

merged_dt$dataset <- NULL

dotcols <- all_cols[all_cmps[merged_dt$exprds]]

xvar <- "rd_fcc_auc"
yvar <- "ratioOverThresh_obs"

for(yvar in c("ratioOverThresh_obsOverRd", "ratioOverThresh_obs")) {
  
  outFile <- file.path(outFolder, paste0(yvar, "_", xvar, "_meanCorrPermut_v2_either.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  
  plot(x=merged_dt[,paste0(xvar)],
       y=merged_dt[,paste0(yvar)],
       cex=0.7,
       col=dotcols,
       cex.main=plotCex,
       cex.lab=plotCex,
       cex.axis=plotCex,
       main="MEANCORRPERMUT_V2_EITHER",
       xlab=xvar,
       ylab=yvar,
       pch=16)
  addCorr(x=merged_dt[,paste0(xvar)],
          y=merged_dt[,paste0(yvar)],bty="n", legPos="topleft")  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
}

####################################################################################################################
####################################################################################################################
####################################################################################################################

nOverQt_dt <- get(load("N_FCC_OVER_QTILE_RANDOMMIDPOSSTRICT//nOverQt_dt.Rdata"))
aucRatio_dt <- get(load("../FIGURES_V5_YUANLONG/RANDOM_FCC_AUC_RATIO_RANDOMMIDPOS_V2/all_auc_ratio_dt_RANDOMMIDPOSSTRICT.Rdata"))

merged_dt <- merge(nOverQt_dt, aucRatio_dt, by=c("hicds", "exprds"))

merged_dt$ratioOverThresh_obsOverRd <- merged_dt$ratioOverThresh_obs/merged_dt$ratioOverThresh_rd

merged_dt$dataset <- NULL

dotcols <- all_cols[all_cmps[merged_dt$exprds]]

xvar <- "rd_fcc_auc"
yvar <- "ratioOverThresh_obs"


for(yvar in c("ratioOverThresh_obsOverRd", "ratioOverThresh_obs")) {
  
  outFile <- file.path(outFolder, paste0(yvar, "_", xvar, "_randommidposstrict.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  
  
  plot(x=merged_dt[,paste0(xvar)],
       y=merged_dt[,paste0(yvar)],
       cex=0.7,
       col=dotcols,
       cex.main=plotCex,
       cex.lab=plotCex,
       cex.axis=plotCex,
       main="RANDOMMIDPOSSTRICT",
       xlab=xvar,
       ylab=yvar,
       pch=16)
  
  addCorr(x=merged_dt[,paste0(xvar)],
          y=merged_dt[,paste0(yvar)],bty="n", legPos="topleft")
  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
}

####################################################################################################################
####################################################################################################################
####################################################################################################################


nOverQt_dt <- get(load("N_FCC_OVER_QTILE_PERMG2T//nOverQt_dt.Rdata"))
aucRatio_dt <- get(load("../FIGURES_V3_YUANLONG/BARPLOT_FCC_AUC_RATIO/all_dt.Rdata"))

merged_dt <- merge(nOverQt_dt, aucRatio_dt, by=c("hicds", "exprds"))

merged_dt$ratioOverThresh_obsOverRd <- merged_dt$ratioOverThresh_obs/merged_dt$ratioOverThresh_rd

merged_dt$dataset <- NULL

dotcols <- all_cols[all_cmps[merged_dt$exprds]]

xvar <- "fcc_auc"
yvar <- "ratioOverThresh_obs"

for(yvar in c("ratioOverThresh_obsOverRd", "ratioOverThresh_obs")) {
  
  outFile <- file.path(outFolder, paste0(yvar, "_", xvar, "_permg2t.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  
  plot(x=merged_dt[,paste0(xvar)],
       y=merged_dt[,paste0(yvar)],
       cex=0.7,
       col=dotcols,
       cex.main=plotCex,
       cex.lab=plotCex,
       cex.axis=plotCex,
       main="PERMG2T",
       xlab=xvar,
       ylab=yvar,
       pch=16)
  addCorr(x=merged_dt[,paste0(xvar)],
          y=merged_dt[,paste0(yvar)],bty="n", legPos="topleft")
  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
}
####################################################################################################################
####################################################################################################################
####################################################################################################################