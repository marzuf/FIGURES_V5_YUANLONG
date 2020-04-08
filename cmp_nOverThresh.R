source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("../FIGURES_V2_YUANLONG/settings.R")

outFolder <- "CMP_NOVERQT"
dir.create(outFolder)

plotCex <- 1.4

v2_either_dt <- get(load("N_FCC_OVER_QTILE_V2_EITHER/nOverQt_dt.Rdata"))
colnames(v2_either_dt)[3:ncol(v2_either_dt)] <- paste0(colnames(v2_either_dt)[3:ncol(v2_either_dt)], "_v2_either")

randommidpos_dt <- get(load("N_FCC_OVER_QTILE_RANDOMMIDPOSSTRICT//nOverQt_dt.Rdata"))
colnames(randommidpos_dt)[3:ncol(randommidpos_dt)] <- paste0(colnames(randommidpos_dt)[3:ncol(randommidpos_dt)], "_randommidposstrict")

merged_dt <- merge(v2_either_dt, randommidpos_dt, by=c("hicds", "exprds"))

plot_var <- "ratioOverQt_obs"

dotcols <- all_cols[all_cmps[merged_dt$exprds]]

outFile <- file.path(outFolder, paste0(plot_var, "randommidposstrict_vs_v2_either"))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(
  x = merged_dt[,paste0(plot_var, "_v2_either")],
  y = merged_dt[,paste0(plot_var, "_randommidposstrict")],
  xlab="v2_either",
  ylab="randommidposstrict",
  main=plot_var,
  cex.main=plotCex,
  cex.lab=plotCex,
  cex.axis=plotCex,
  pch=16,
  col=dotcols,
  cex=0.7
)
addCorr(
  x = merged_dt[,paste0(plot_var, "_v2_either")],
  y = merged_dt[,paste0(plot_var, "_randommidposstrict")],
  bty="n",
  legPos="bottomright"
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))
####################################################################################################################

v2_either_dt <- get(load("N_FCC_OVER_QTILE_V2_EITHER/nOverQt_dt.Rdata"))
colnames(v2_either_dt)[3:ncol(v2_either_dt)] <- paste0(colnames(v2_either_dt)[3:ncol(v2_either_dt)], "_v2_either")

permg2t_dt <- get(load("N_FCC_OVER_QTILE_PERMG2T//nOverQt_dt.Rdata"))
colnames(permg2t_dt)[3:ncol(permg2t_dt)] <- paste0(colnames(permg2t_dt)[3:ncol(permg2t_dt)], "_permg2t")

merged_dt <- merge(v2_either_dt, permg2t_dt, by=c("hicds", "exprds"))

plot_var <- "ratioOverQt_obs"

dotcols <- all_cols[all_cmps[merged_dt$exprds]]

outFile <- file.path(outFolder, paste0(plot_var, "permg2t_vs_v2_either"))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(
  x = merged_dt[,paste0(plot_var, "_v2_either")],
  y = merged_dt[,paste0(plot_var, "_permg2t")],
  xlab="v2_either",
  ylab="permg2t",
  main=plot_var,
  cex.main=plotCex,
  cex.lab=plotCex,
  cex.axis=plotCex,
  pch=16,
  col=dotcols,
  cex=0.7
)
addCorr(
  x = merged_dt[,paste0(plot_var, "_v2_either")],
  y = merged_dt[,paste0(plot_var, "_permg2t")],
  bty="n",
  legPos="bottomright"
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

####################################################################################################################

randommidpos_dt <- get(load("N_FCC_OVER_QTILE_RANDOMMIDPOSSTRICT//nOverQt_dt.Rdata"))
colnames(randommidpos_dt)[3:ncol(randommidpos_dt)] <- paste0(colnames(randommidpos_dt)[3:ncol(randommidpos_dt)], "_randommidpos")

permg2t_dt <- get(load("N_FCC_OVER_QTILE_PERMG2T//nOverQt_dt.Rdata"))
colnames(permg2t_dt)[3:ncol(permg2t_dt)] <- paste0(colnames(permg2t_dt)[3:ncol(permg2t_dt)], "_permg2t")

merged_dt <- merge(randommidpos_dt, permg2t_dt, by=c("hicds", "exprds"))

plot_var <- "ratioOverQt_obs"

dotcols <- all_cols[all_cmps[merged_dt$exprds]]

outFile <- file.path(outFolder, paste0(plot_var, "permg2t_vs_randommidpos"))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(
  x = merged_dt[,paste0(plot_var, "_randommidpos")],
  y = merged_dt[,paste0(plot_var, "_permg2t")],
  xlab="randommidpos",
  ylab="permg2t",
  main=plot_var,
  cex.main=plotCex,
  cex.lab=plotCex,
  cex.axis=plotCex,
  pch=16,
  col=dotcols,
  cex=0.7
)
addCorr(
  x = merged_dt[,paste0(plot_var, "_randommidpos")],
  y = merged_dt[,paste0(plot_var, "_permg2t")],
  bty="n",
  legPos="bottomright"
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))



