
options(scipen=100)

# Rscript random_FCC_AUC_ratio_randommidpos_v2.R

script_name <- "random_FCC_AUC_ratio_randommidpos_v2.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

buildData <- FALSE

require(flux)
require(foreach)
require(doMC)
require(reshape2)
require(ggplot2)
require(ggpubr)
require(ggsci)
registerDoMC(40)

outFolder <- file.path("RANDOM_FCC_AUC_RATIO_RANDOMMIDPOS_V2")
dir.create(outFolder, recursive = TRUE)

runFolder <- file.path("..", "v2_Yuanlong_Cancer_HiC_data_TAD_DA")

permutFolder <- file.path(runFolder, "PIPELINE", "OUTPUT_FOLDER")

plotMargin <- c(1,2,1,1)

ggsci_pal <- "lancet"
ggsci_subpal <- ""


source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

plotType <- "svg"
myHeight <- 7
myWidth <- 10
source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("../FIGURES_V2_YUANLONG/settings.R")

cumsumCurve_col <- "darkgrey"



fcc_auc_dt <- get(load("../FIGURES_V3_YUANLONG/BARPLOT_FCC_AUC_RATIO/all_dt.Rdata"))


script8_name <- "8cOnlyFCC_runAllDown"


permutFolder <- file.path("../v2_Yuanlong_Cancer_HiC_data_TAD_DA/PIPELINE/OUTPUT_FOLDER")
all_hicds <- list.files(file.path(permutFolder))
all_hicds <- all_hicds[! (grepl("RANDOM", all_hicds) | grepl("PERMUT", all_hicds))]
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(permutFolder, x)))

all_rd <- c("RANDOMMIDPOS", "RANDOMMIDPOSDISC", "RANDOMMIDPOSSTRICT")
# all_rd <- c("RANDOMMIDPOSSTRICT")

rd = all_rd [1]
for(rd in all_rd){
  hicds = all_hicds[1]
    all_auc_ratio_dt <- foreach(hicds = all_hicds, .combine='rbind') %do%{
      exprds = all_exprds[[paste0(hicds)]][1]
      hicds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
        
        cat(paste0("... start ", rd, " - ", hicds, " - ", exprds, "\n"))
        
        rd_hicds <- gsub("_40kb", paste0("_", rd, "_40kb"), hicds)
        
        fcc_file <- file.path(permutFolder, hicds, exprds, script8_name, "all_obs_prodSignedRatio.Rdata")
        stopifnot(file.exists(fcc_file))
        
        rd_fcc_file <- file.path(permutFolder, rd_hicds, exprds, script8_name, "all_obs_prodSignedRatio.Rdata")
        stopifnot(file.exists(rd_fcc_file))
        
        obs_fcc <- get(load(fcc_file))
        rd_fcc <- get(load(rd_fcc_file))
        
        
        obs_fcc_sorted <- sort(obs_fcc, decreasing = TRUE)
        rd_fcc_sorted <- sort(rd_fcc, decreasing = TRUE)
        
        nTot  <- min(c(length(obs_fcc_sorted), length(rd_fcc_sorted)))
        x_val <- 1:nTot
        
        
        cumsum_obs <- cumsum(obs_fcc_sorted[1:nTot])
        cumsum_rd <- cumsum(rd_fcc_sorted[1:nTot])
        
        auc_obs <- auc(x = x_val, y = cumsum_obs)
        auc_rd <- auc(x = x_val, y = cumsum_rd)
        
        auc_ratio_rd <- auc_obs/auc_rd
        
        
        outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_cumsumFCC_obs_", rd,".", plotType))
        do.call(plotType, list(outFile, height=myHeight, width=myWidth))
        par(bty="L")
        plot(
          x = 1:nTot,
          y= cumsum_obs,
          xlab = "TADs ranked by decreasing FCC",
          ylab = "Cumsum FCC",
          ylim = range(c(cumsum_obs, cumsum_rd)),
          main=paste0("Cumsum FCC scores obs-", rd),
          type="l",
          cex.lab=plotCex,
          cex.main=plotCex,
          cex.axis=plotCex
        )
        lines(
          x=1:nTot,
          y=cumsum_rd,
          col = cumsumCurve_col
        )
        mtext(side=3, text=paste0(hicds, " - ", exprds))
        
        # what I could do: merge all curve and take 95-qtile from all permuts
        
        legend(
          "topleft",
          paste0("# obs. TADs = ",  length(obs_fcc_sorted), "\n# ", rd, " TADs  = ", length(rd_fcc_sorted)),
          bty="n"
        )
        legend(
          "bottomright",
          lty=c(-1, 1, 1),
          legend=c(paste0("AUC ratio = ", round(auc_ratio_rd, 2)), "obs.", rd),
          col = c("black", "black", cumsumCurve_col),
          bty="n"
        )
        foo <- dev.off()
        cat(paste0("... written: ", outFile, "\n"))
        
        data.frame(
          hicds =hicds,
          exprds=exprds,
          rd_fcc_auc=auc_ratio_rd,
          stringsAsFactors = FALSE
        )
        
        
      }
      hicds_dt
    }

  save(all_auc_ratio_dt,file = file.path(outFolder, paste0("all_auc_ratio_dt_", rd, ".Rdata")), version=2)
  cat(paste0("... written all_auc_ratio_dt.Rdata\n"))
  
  plot_dt <- merge(fcc_auc_dt, all_auc_ratio_dt, by =c("hicds", "exprds"), all.x=TRUE, all.y=TRUE)
  stopifnot(!is.na(plot_dt))
  
  dotcols <- all_cols[all_cmps[plot_dt$exprds]]

  outFile <- file.path(outFolder, paste0("allDS_FCC_aucratio_obs_vs_", rd,".", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myHeight))
  par(bty="L")
  plot(
    x = plot_dt$fcc_auc,
    y= plot_dt$rd_fcc_auc,
    xlab = "FCC AUC ratio (observed)",
    ylab = paste0("FCC AUC ratio (", rd, ")"),
    main=paste0("FCC AUC ratio obs-", rd),
	col = dotcols,
    pch=16,
    cex=0.7,
    cex.lab=plotCex,
    cex.main=plotCex,
    cex.axis=plotCex
  )
  mtext(side=3, text=paste0("all DS (n=", nrow(plot_dt), ")"))
  
  addCorr(
    x = plot_dt$fcc_auc,
    y= plot_dt$rd_fcc_auc,
    bty="n"
  )
legend("topleft", 
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

