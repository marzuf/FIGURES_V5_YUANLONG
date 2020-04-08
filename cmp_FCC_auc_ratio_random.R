# Rscript cmp_FCC_auc_ratio_random.R

outFolder <- "CMP_FCC_AUC_RATIO_RANDOM"
dir.create(outFolder)


plotType <- "svg"
myHeight <- myWidth <- 7
plotCex <- 1.4

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("../FIGURES_V2_YUANLONG/settings.R")

auc_permg2t_dt <- get(load("../FIGURES_V3_YUANLONG/BARPLOT_FCC_AUC_RATIO/all_dt.Rdata"))
colnames(auc_permg2t_dt)[colnames(auc_permg2t_dt) == "fcc_auc"] <- "PERMG2T_FCC_AUC"

auc_randommidposstrict_dt <- get(load("RANDOM_FCC_AUC_RATIO_RANDOMMIDPOS_V2/all_auc_ratio_dt_RANDOMMIDPOSSTRICT.Rdata"))
colnames(auc_randommidposstrict_dt)[colnames(auc_randommidposstrict_dt) == "rd_fcc_auc"] <- "RANDOMMIDPOSSTRICT_FCC_AUC"
auc_randommidposdisc_dt <- get(load("RANDOM_FCC_AUC_RATIO_RANDOMMIDPOS_V2/all_auc_ratio_dt_RANDOMMIDPOSDISC.Rdata"))
colnames(auc_randommidposdisc_dt)[colnames(auc_randommidposdisc_dt) == "rd_fcc_auc"] <- "RANDOMMIDPOSDISC_FCC_AUC"
auc_randommidpos_dt <- get(load("RANDOM_FCC_AUC_RATIO_RANDOMMIDPOS_V2/all_auc_ratio_dt_RANDOMMIDPOS.Rdata"))
colnames(auc_randommidpos_dt)[colnames(auc_randommidpos_dt) == "rd_fcc_auc"] <- "RANDOMMIDPOS_FCC_AUC"

auc_meanCorrPermut_data <- get(load("../FIGURES_V3_YUANLONG/RANDOM_FCC_AUC_RATIO_MEANCORRPERMUT_V2/all_permut_fcc.Rdata"))
auc_meanCorrPermut_data2 <- unlist(lapply(auc_meanCorrPermut_data, function(sub) lapply(sub, function(x) x[["auc_ratio_rd_meanRL"]])))
auc_meanCorrPermut_dt <- data.frame(
  hicds = gsub("(.+)\\..+", "\\1", names(auc_meanCorrPermut_data2)),
  exprds = gsub(".+\\.(.+)", "\\1", names(auc_meanCorrPermut_data2)),
  MEANCORRPERMUT_meanRL_FCC_AUC = as.numeric(auc_meanCorrPermut_data2),
  stringsAsFactors = FALSE
)

auc_meanCorrPermut_data <- get(load("../FIGURES_V3_YUANLONG/RANDOM_FCC_AUC_RATIO_MEANCORRPERMUT_V3/all_permut_fcc.Rdata"))
auc_meanCorrPermut_data2 <- unlist(lapply(auc_meanCorrPermut_data, function(sub) lapply(sub, function(x) x[["auc_ratio_rd_meanRL"]])))
auc_meanCorrPermut_withTAD_dt <- data.frame(
  hicds = gsub("(.+)\\..+", "\\1", names(auc_meanCorrPermut_data2)),
  exprds = gsub(".+\\.(.+)", "\\1", names(auc_meanCorrPermut_data2)),
  MEANCORRPERMUT_meanRL_withTAD_FCC_AUC = as.numeric(auc_meanCorrPermut_data2),
  stringsAsFactors = FALSE
)


auc_meanCorrPermutRorL_data <- get(load("../FIGURES_V3_YUANLONG/RANDOM_FCC_AUC_RATIO_MEANCORRPERMUT_V2_EITHER/all_permut_fcc.Rdata"))
auc_meanCorrPermutRorL_data2 <- unlist(lapply(auc_meanCorrPermutRorL_data, function(sub) lapply(sub, function(x) x[["auc_ratio_rd_RorL"]])))
auc_meanCorrPermutRorL_dt <- data.frame(
  hicds = gsub("(.+)\\..+", "\\1", names(auc_meanCorrPermutRorL_data2)),
  exprds = gsub(".+\\.(.+)", "\\1", names(auc_meanCorrPermutRorL_data2)),
  MEANCORRPERMUT_RorL_FCC_AUC = as.numeric(auc_meanCorrPermutRorL_data2),
  stringsAsFactors = FALSE
)

auc_meanCorrPermutRorL_data <- get(load("../FIGURES_V3_YUANLONG/RANDOM_FCC_AUC_RATIO_MEANCORRPERMUT_V3_EITHER/all_permut_fcc.Rdata"))
auc_meanCorrPermutRorL_data2 <- unlist(lapply(auc_meanCorrPermutRorL_data, function(sub) lapply(sub, function(x) x[["auc_ratio_rd_RorL"]])))
auc_meanCorrPermutRorL_withTAD_dt <- data.frame(
  hicds = gsub("(.+)\\..+", "\\1", names(auc_meanCorrPermutRorL_data2)),
  exprds = gsub(".+\\.(.+)", "\\1", names(auc_meanCorrPermutRorL_data2)),
  MEANCORRPERMUT_RorL_withTAD_FCC_AUC = as.numeric(auc_meanCorrPermutRorL_data2),
  stringsAsFactors = FALSE
)

auc_meanCorrPermutRorL_data <- get(load("../FIGURES_V3_YUANLONG/RANDOM_FCC_AUC_RATIO_MEANCORRPERMUT_V4_EITHER/all_permut_fcc.Rdata"))
auc_meanCorrPermutRorL_data2 <- unlist(lapply(auc_meanCorrPermutRorL_data, function(sub) lapply(sub, function(x) x[["auc_ratio_rd_RorL"]])))
auc_meanCorrPermutRorL_replace_dt <- data.frame(
  hicds = gsub("(.+)\\..+", "\\1", names(auc_meanCorrPermutRorL_data2)),
  exprds = gsub(".+\\.(.+)", "\\1", names(auc_meanCorrPermutRorL_data2)),
  MEANCORRPERMUT_RorL_replace_FCC_AUC = as.numeric(auc_meanCorrPermutRorL_data2),
  stringsAsFactors = FALSE
)



auc_onePerm_data <- get(load("../FIGURES_V3_YUANLONG/RANDOM_FCC_AUC_RATIO_ONEPERM/all_permut_fcc.Rdata"))
auc_onePerm_data2 <- unlist(lapply(auc_onePerm_data, function(sub) lapply(sub, function(x) x[["auc_ratio_rd_onePermut"]])))
auc_onePerm_dt <- data.frame(
  hicds = gsub("(.+)\\..+", "\\1", names(auc_onePerm_data2)),
  exprds = gsub(".+\\.(.+)", "\\1", names(auc_onePerm_data2)),
  ONEPERMG2T_FCC_AUC = as.numeric(auc_onePerm_data2),
  stringsAsFactors = FALSE
)


plot_dt <- merge(auc_meanCorrPermut_dt, 

			merge(auc_meanCorrPermut_withTAD_dt,
			
			merge(auc_meanCorrPermutRorL_dt, 

			merge(auc_meanCorrPermutRorL_withTAD_dt,

			merge(auc_meanCorrPermutRorL_replace_dt,

merge(auc_permg2t_dt,

                 merge( auc_randommidposstrict_dt, 
						merge(auc_randommidposdisc_dt, 

merge(auc_onePerm_dt,auc_randommidpos_dt, 

								by = c("hicds", "exprds"), all=TRUE),  by = c("hicds", "exprds"), all=TRUE),  
by = c("hicds", "exprds"), all=TRUE) , by = c("hicds", "exprds"), all=TRUE) , by = c("hicds", "exprds"), all=TRUE) , 
								 by = c("hicds", "exprds"), all=TRUE) , by = c("hicds", "exprds"), all=TRUE),by = c("hicds", "exprds"), all=TRUE), by = c("hicds", "exprds"), all=TRUE)

outFile <- file.path(outFolder, "plot_dt.Rdata")
save(plot_dt, file=outFile, version=2)

stopifnot(!is.na(plot_dt))

dotcols <- all_cols[all_cmps[plot_dt$exprds]]

for(r1 in 3:(ncol(plot_dt)-1)) {
 
  var1 <- colnames(plot_dt)[r1]
   
  for(r2 in (r1+1):ncol(plot_dt)) {
  
    var2 <- colnames(plot_dt)[r2]
   
    myx <-  plot_dt[,c(var1)]
    myy <-  plot_dt[,c(var2)]
    
    outFile <- file.path(outFolder, paste0(var2, "_vs_", var1, "_FCC_aucratio.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myHeight))
    par(bty="L")
    plot(
      x = myx,
      y= myy,
      xlab = paste0(var1),
      ylab = paste0(var2),
      main=paste0("CMP FCC AUC ratio"),
      pch=16,
      cex=0.7,
col = dotcols,
      cex.lab=plotCex,
      cex.main=plotCex,
      cex.axis=plotCex
    )
    mtext(side=3, text=paste0("all DS (n=", nrow(plot_dt), ")"))
    
    addCorr(
      x = myx,
      y= myy,
      legPos="topleft",
      bty="n"
    )
legend("topright", 
       # legend = paste0(labsymbol, " ", names(all_cols)),
       legend = paste0(names(all_cols)),
       col=all_cols,
       pch=16,
#       cex = plotCex,
 bty="n"
)
    
    
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    
    
  }
}

load("CMP_FCC_AUC_RATIO_RANDOM/plot_dt.Rdata")
xdt <- plot_dt[,c("hicds", "exprds", "RANDOMMIDPOSSTRICT_FCC_AUC", "PERMG2T_FCC_AUC")]
xdt[xdt$RANDOMMIDPOSSTRICT_FCC_AUC >= 1.45,]
# ENCSR504OTV_transverse_colon_40kb TCGAcoad_msi_mss
xdt[xdt$PERMG2T_FCC_AUC >= 1.5 & xdt$RANDOMMIDPOSSTRICT_FCC_AUC <= 1.2,]
# ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_wt_mutCTNNB1



