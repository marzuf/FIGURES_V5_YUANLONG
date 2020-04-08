# NEW VERSION: RANK BY MINIMAL # OF SAMPLES OR BY SUM # OF SAMPLES


plotType <- "svg"


source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("../FIGURES_V2_YUANLONG/settings.R")
require(ggsci)

# Rscript nSignifGenes_nSignifTADs_sumSamp_byCond.R

outFolder <- file.path("NSIGNIFGENES_NSIGNIFTADS_SUMSAMP_BYCOND")
dir.create(outFolder, recursive = TRUE)

outFile_model <- file.path(outFolder, "outFile_model.txt")
file.remove(outFile_model)

dotpch <- 19
# segcol <-  "#BEBEBE19"
segcol <- "grey"
dotCex <- 1.1

setDir <- "/media/electron"
setDir <- ""
mainFolder <- file.path("../v2_Yuanlong_Cancer_HiC_data_TAD_DA")
settingFolder <- file.path(mainFolder, "PIPELINE", "INPUT_FILES")


inDT <- get(load("../v2_Yuanlong_Cancer_HiC_data_TAD_DA/GENE_RANK_TAD_RANK/all_gene_tad_signif_dt.Rdata"))
geneDT <- inDT[,c("hicds", "exprds", "entrezID", "adj.P.Val")]
geneDT <- unique(geneDT)
nSignifGenes_dt <- aggregate(adj.P.Val~hicds + exprds, data = geneDT, FUN=function(x) sum(x<=geneSignifThresh))
colnames(nSignifGenes_dt)[colnames(nSignifGenes_dt) == "adj.P.Val"] <- "nSignifGenes"

tadDT <- inDT[,c("hicds", "exprds", "region", "tad_adjCombPval")]
tadDT <- unique(tadDT)
nSignifTADs_dt <- aggregate(tad_adjCombPval~hicds + exprds, data = tadDT, FUN=function(x) sum(x<=tadSignifThresh))
colnames(nSignifTADs_dt)[colnames(nSignifTADs_dt) == "tad_adjCombPval"] <- "nSignifTADs"

nSignif_dt <- merge(nSignifGenes_dt, nSignifTADs_dt, by=c("hicds", "exprds"), all=TRUE)
stopifnot(!is.na(nSignif_dt))

nSignif_dt <- nSignif_dt[order(nSignif_dt$nSignifTADs, decreasing = TRUE),]

# retrieve the number of samples

nSignif_dt$dataset <- file.path(nSignif_dt$hicds, nSignif_dt$exprds)

nSignif_dt$sumSample <- sapply(nSignif_dt$dataset, function(x) {
  settingFile <- file.path(settingFolder, dirname(x), paste0("run_settings_", basename(x), ".R"))
  stopifnot(file.exists(settingFile))
  source(settingFile)
  samp1 <- get(load(file.path(setDir, sample1_file)))
  samp2 <- get(load(file.path(setDir, sample2_file)))
  sum(c(length(samp1), length(samp2)))
})

nSignif_dt$dataset <- paste0(nSignif_dt$hicds, "\n", nSignif_dt$exprds)


nSignif_dt <- nSignif_dt[order(nSignif_dt$sumSample),]


nSignif_dt$ds_rank <- 1:nrow(nSignif_dt)

nSignif_dt$cond <- all_cmps[nSignif_dt$exprds]

nSignif_dt$labcols <- all_cols[all_cmps[nSignif_dt$exprds]]


for(x_var in c("sumSample", "ds_rank")) {
  for(y_var in c("nSignifGenes", "nSignifTADs")) {
    
    m_wt_mut <- lm(as.formula(paste0(y_var, "~", x_var)), data = nSignif_dt[nSignif_dt$cond == "wt_vs_mut",])
    m_subtypes <- lm(as.formula(paste0(y_var, "~", x_var)), data = nSignif_dt[nSignif_dt$cond == "subtypes",])
    m_norm_tum <- lm(as.formula(paste0(y_var, "~", x_var)), data = nSignif_dt[nSignif_dt$cond == "norm_vs_tumor",])
    m_all <- lm(as.formula(paste0(y_var, "~", x_var)), data = nSignif_dt)
    
    outFile <- file.path(outFolder, paste0(y_var, "_vs_", x_var, "_geneSignif",geneSignifThresh, "_tadSignif", tadSignifThresh, "_withFit.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    
    my_x <- nSignif_dt[,paste0(x_var)]
    my_y <- nSignif_dt[,paste0(y_var)]
    par(bty="L")
    plot(
      x = my_x,
      y = my_y,
      main  = paste0(y_var, "_vs_", x_var),
      xlab = x_var,
      ylab=y_var,
      cex=0.7,
      pch=16,
      col=nSignif_dt$labcols,
      cex.axis=plotCex,
      cex.lab=plotCex,
      cex.main=plotCex
    )
    if(y_var == "nSignifGenes") {
      mtext(side=3, text = paste0("all DS (n=", nrow(nSignif_dt), "); gene adj. p-val<=", geneSignifThresh ))  
    } else if(y_var == "nSignifTADs") {
      mtext(side=3, text = paste0("all DS (n=", nrow(nSignif_dt), "); TAD adj. p-val_tadSignif", tadSignifThresh))  
      
    }
    
    
    
    abline(m_wt_mut, col = all_cols["wt_vs_mut"], lty=2)
    abline(m_subtypes, col = all_cols["subtypes"], lty=2)
    abline(m_norm_tum, col = all_cols["norm_vs_tumor"], lty=2)
    abline(m_all, col ="black", lty=2)
    
    legend(
      "topleft",
      legend = c(names(all_cols), "all"),
      col = c(all_cols, "black"),
      bty="n",
      lty=1
    )
    
  }
}



m1 <- lm(nSignifTADs ~ ds_rank, data = nSignif_dt)
m2 <- lm(nSignifGenes ~ ds_rank, data = nSignif_dt)
m3 <- lm(nSignifTADs ~ ds_rank+cond, data = nSignif_dt)
m4 <- lm(nSignifGenes ~ ds_rank+cond, data = nSignif_dt)

sink(outFile_model, append=TRUE)
print("************************** ALL DATASETS\n")
print(summary(m1))
print("\n----------------------------------------------------\n")
print(summary(m2))
print("\n----------------------------------------------------\n")
print(summary(m3))
print("\n----------------------------------------------------\n")
print(summary(m4))
sink()
    


init_nSignif_dt <- nSignif_dt

all_conds <- unique(nSignif_dt$cond)
curr_cond = all_conds[1]
for(curr_cond in all_conds) {

	
	nSignif_dt <- init_nSignif_dt[init_nSignif_dt$cond == curr_cond,]


	labcols <- nSignif_dt$labcols
	
m1 <- lm(nSignifTADs ~ ds_rank, data = nSignif_dt)
m2 <- lm(nSignifGenes ~ ds_rank, data = nSignif_dt)

sink(outFile_model, append=TRUE)
print(paste0("************************** ", curr_cond, "\n"))
print(summary(m1))
print("\n----------------------------------------------------\n")
print(summary(m2))
sink()


labcols <- all_cols[all_cmps[nSignif_dt$exprds]]

maxTADs <- max(ceiling(nSignif_dt$nSignifTADs/10)*10)
maxGenes <- max(ceiling(nSignif_dt$nSignifGenes/1000)*1000)

nSignif_dt$nSignifTADs_rescaled <- nSignif_dt$nSignifTADs/maxTADs * maxGenes
nSignif_dt$nSignifGenes_rescaled <- nSignif_dt$nSignifGenes/maxGenes * maxTADs

outFile <- file.path(outFolder, paste0("nSignifGenes_nSignifTADs_all_ds_withSymb_geneSignif",geneSignifThresh, "_tadSignif", tadSignifThresh, "_", curr_cond, ".", plotType))
do.call(plotType, list(outFile, height=myHeight*1.2, width=myWidth*2))
dev.control(displaylist="enable")
par(bty="U", family=fontFamily)
par(mar = c(5,5,2,5))
plot(
  NULL,
  # x = 1:nrow(nSignif_dt),
  # y = nSignif_dt$nSignifTADs, 
  # col = tad_signif_col,
  xlim = c(1, nrow(nSignif_dt)),
  ylim=c(0, maxTADs),
  ylab = "# signif. TADs",
  xlab = "",
  pch=dotpch,
  cex = dotCex,
  main = paste0("# signif. features"),
  cex.axis=plotCex,
  cex.main = plotCex,
  cex.lab = plotCex,
  col.lab = tad_signif_col,
  axes = FALSE
)


segments(
  x0= 1:nrow(nSignif_dt),
  y0 = nSignif_dt$nSignifGenes_rescaled,
  x1= 1:nrow(nSignif_dt),
  y1=nSignif_dt$nSignifTADs,
  col = segcol
)

points(
  x = 1:nrow(nSignif_dt),
  y = nSignif_dt$nSignifTADs, 
  col = tad_signif_col,
  cex = dotCex,
  pch=dotpch
)


abline(lm(nSignif_dt$nSignifTADs~c(1:nrow(nSignif_dt))), col = tad_signif_col, lty=2)

legend("topleft",
       legend=c(paste0("gene signif.: adj. p-val <= ", geneSignifThresh),
                paste0("TAD signif.: adj. p-val <= ", tadSignifThresh)),
       bty="n")



mtext(side=3, line=-1, text=paste0("all datasets - n= ", length(unique(nSignif_dt$dataset))))
# mtext(side=1, col = labcols, text = nSignif_dt$dataset, at= 1:nrow(nSignif_dt), las=2, cex =0.5)
axis(2, col = tad_signif_col, col.ticks = tad_signif_col, col.axis=tad_signif_col, at=seq(from=0, to = maxTADs, by=10))
axis(1, labels=F, lwd.ticks = -1)
# axis(1, labels=F, at=1:nrow(nSignif_dt))
par(new = T, family=fontFamily)
plot(
  x = 1:nrow(nSignif_dt),
  y = nSignif_dt$nSignifGenes,
  col = gene_signif_col,
  ylab = NA,
  ylim=c(0, maxGenes),
  xlab = NA,
  pch=dotpch,
  cex.axis=plotCex,
  cex.main = plotCex,
  cex.lab = plotCex,
  axes = FALSE
)
axis(side=4, col = gene_signif_col, col.ticks = gene_signif_col, col.axis=gene_signif_col, at = seq(from=0, to=maxGenes, by=1000))
mtext(side = 4, line = 3, '# signif. genes', col=gene_signif_col,  cex=plotCex)

abline(lm(nSignif_dt$nSignifGenes~c(1:nrow(nSignif_dt))), col = gene_signif_col, lty=2)


legend("bottom", 
       # legend = paste0(labsymbol, " ", names(all_cols)),
       legend = paste0(names(all_cols)),
       col=all_cols,
       pch=15,
       # lty=c(1,2),
       horiz=TRUE,
       inset=c(0,-0.12),
       cex = plotCex,
       xpd=TRUE, bty="n"
)
signifPlot <- recordPlot()



invisible(dev.off())
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("nSignifGenes_nSignifTADs_all_ds_withLeg_geneSignif",geneSignifThresh, "_tadSignif", tadSignifThresh, "_", curr_cond, ".", plotType))
do.call(plotType, list(outFile, height=myHeight*1.2, width=myWidth*2.5))
par(mar = c(5,5,2,5))
replayPlot(signifPlot) 
mtext(side=1, col = labcols, text = nSignif_dt$dataset, at= 1:nrow(nSignif_dt), las=2, cex = 0.6)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))




outFile <- file.path(outFolder, paste0("nSignifGenes_nSignifTADs_all_ds_withNbrSamp_geneSignif",geneSignifThresh, "_tadSignif", tadSignifThresh, "_", curr_cond, ".", plotType))
do.call(plotType, list(outFile, height=myHeight*1.2, width=myWidth*2.5))
par(mar = c(5,5,2,5))
replayPlot(signifPlot) 
mtext(side=1, col = labcols, text = paste0(nSignif_dt$sumSample, " ", labsymbol), at= 1:nrow(nSignif_dt),  las=2,cex = 0.9, line = 0)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


outFile <- file.path(outFolder, paste0("nSignifGenes_nSignifTADs_all_ds_withNbrSamp_geneSignif",geneSignifThresh, "_tadSignif", tadSignifThresh, "_", curr_cond, ".", plotType))
do.call(plotType, list(outFile, height=myHeight*1.2, width=myWidth*2.5))
par(mar = c(5,5,2,5))
replayPlot(signifPlot) 
mtext(side=1, col = labcols, text = paste0(nSignif_dt$sumSample, " ", labsymbol), at= 1:nrow(nSignif_dt),  las=2,cex = 0.9, line = 0)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

dotcols <- all_cols[all_cmps[nSignif_dt$exprds]]

outFile <- file.path(outFolder, paste0("nSignifGenes_vs_nSignifTADs_all_ds_",geneSignifThresh, "_tadSignif", tadSignifThresh, "_", curr_cond, ".", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
par(bty="L")
plot(
  x = nSignif_dt$nSignifTADs,
  y = nSignif_dt$nSignifGenes,
  col = dotcols,
  ylab = "# signif. genes",
  xlab = "# signif. TADs",
  pch=16,
cex=0.7,
	main = paste0("# signif. features"),
  cex.axis=plotCex,
  cex.main = plotCex,
  cex.lab = plotCex
)
mtext(side=3, text = paste0(curr_cond, " - n=", nrow(nSignif_dt), "; gene adj. p-val <= ", geneSignifThresh, " - TAD adj. p-val <= ", tadSignifThresh))
addCorr(x=nSignif_dt$nSignifTADs,y=nSignif_dt$nSignifGenes,bty="n")
#legend("bottomleft", 
#       # legend = paste0(labsymbol, " ", names(all_cols)),
#       legend = paste0(names(all_cols)),
#       col=all_cols,
#       pch=16,
#       cex = plotCex,
# bty="n"
#)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


outFile <- file.path(outFolder, paste0("nSignifGenes_vs_sumNbrSample_all_ds_",geneSignifThresh, "_tadSignif", tadSignifThresh, "_", curr_cond, ".", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
par(bty="L")
plot(
  x = nSignif_dt$sumSample,
  y = nSignif_dt$nSignifGenes,
  col = dotcols,
  ylab = "# signif. genes",
  xlab = "Sum # of samples",
  pch=16,
cex=0.7,
	main = paste0("# signif. genes vs. sample size"),
  cex.axis=plotCex,
  cex.main = plotCex,
  cex.lab = plotCex
)
mtext(side=3, text = paste0(curr_cond, " - n=", nrow(nSignif_dt), "; gene adj. p-val <= ", geneSignifThresh))
addCorr(x=nSignif_dt$sumSample,y=nSignif_dt$nSignifGenes,bty="n")
#legend("bottomleft", 
#       # legend = paste0(labsymbol, " ", names(all_cols)),
#       legend = paste0(names(all_cols)),
#       col=all_cols,
#       pch=16,
#       cex = plotCex,
# bty="n"
#)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


outFile <- file.path(outFolder, paste0("nSignifTADs_vs_sumNbrSample_all_ds_",geneSignifThresh, "_tadSignif", tadSignifThresh, "_", curr_cond, ".", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
par(bty="L")
plot(
  x = nSignif_dt$sumSample,
  y = nSignif_dt$nSignifTADs,
  col = dotcols,
  ylab = "# signif. TADs",
  xlab = "Sum # of samples",
  pch=16,
cex=0.7,
	main = paste0("# signif. TADs vs. sample size"),
  cex.axis=plotCex,
  cex.main = plotCex,
  cex.lab = plotCex
)
mtext(side=3, text = paste0(curr_cond, " - n=", nrow(nSignif_dt), "; TAD adj. p-val <= ", tadSignifThresh))
addCorr(x=nSignif_dt$sumSample,y=nSignif_dt$nSignifGenes,bty="n")
#legend("bottomleft", 
#       # legend = paste0(labsymbol, " ", names(all_cols)),
#       legend = paste0(names(all_cols)),
#       col=all_cols,
#       pch=16,
#       cex = plotCex,
# bty="n"
#)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))




}
cat(paste0("... written: ", outFile_model, "\n"))



