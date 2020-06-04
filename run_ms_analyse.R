### ------------------------------------------------------------------------ ###
### analyse "multi-species" GA runs ####
### ------------------------------------------------------------------------ ###

library(doParallel)
library(doRNG)
library(GA)
library(ggplot2)
library(cowplot)
library(Cairo)
library(tidyr)
library(dplyr)
library(FLCore)
library(FLash)
library(FLBRP)
library(mseDL)

### ------------------------------------------------------------------------ ###
### collate results from HPC ####
### ------------------------------------------------------------------------ ###
stocks <- read.csv("input/stocks.csv", stringsAsFactors = FALSE)
stocks <- stocks[-29, ]
brps <- readRDS("input/OMs/brps.rds")

#path_scn <- "output/500_50/ms/SSB_idx_rfb_exp_error/whg_bll_lem_ane_jnd/"
path_scn <- "output/500_50/ms/full/whg_bll_lem_ane_jnd/"

### GA results
GA_res <- readRDS(paste0(path_scn, "res.rds"))
### individual MSE runs
GA_runs <- readRDS(paste0(path_scn, "/runs.rds"))
which.max(sapply(GA_runs, "[[", "obj"))
GA_runs[[1]]
### optimised parameters
par_opt <- as.list(GA_runs[[1]]$pars)

### ------------------------------------------------------------------------ ###
### plot GA results ####
### ------------------------------------------------------------------------ ###

png(paste0(path_scn, "res.png"), 
    width = 15, height = 10, units = "cm", 
    res = 300, type = "cairo")
plot(GA_res)
dev.off()

### ------------------------------------------------------------------------ ###
### recreate MSE runs ####
### ------------------------------------------------------------------------ ###
### default, optimised

### run optimised solution
library(doParallel)
cl <- makeCluster(10)
registerDoParallel(cl)
cl_length <- length(cl)
req_pckgs <- c("FLCore", "FLash", "mseDL", "GA", "doParallel", "doRNG", "FLBRP")
for (i in req_pckgs) library(package = i, character.only = TRUE)
source("funs.R"); source("GA_funs.R")
. <- foreach(i = seq(cl_length)) %dopar% {
  for (i in req_pckgs) library(package = i, character.only = TRUE)
  source("funs.R"); source("GA_funs.R")
}
### load stocks
stocks <- read.csv("input/stocks.csv", stringsAsFactors = FALSE)[-29, ]
stock <- stocks$stock[c(22,23,24,25,26)]
names(stock) <- stock
input <- lapply(stock, function(x) {
  readRDS(paste0("input/500_50/OM_2_mp_input/", x,
                 ".rds"))
})

# scenario <- "SSB_idx_rfb_exp_error"
scenario <- "full"
### prepare input object(s) for MSE
input <- lapply(input, function(x) {
  ### OEM: activate uncertainty
  x$oem@args$idx_dev <- TRUE
  x$oem@args$ssb <- TRUE
  x$oem@args$lngth <- TRUE
  x$oem@args$lngth_dev <- TRUE
  ### IEM: do not activate uncertainty
  x$iem@args$use_dev <- FALSE
  ### catch rule components
  x$ctrl.mp$ctrl.est@args$comp_r <- TRUE
  x$ctrl.mp$ctrl.est@args$comp_f <- TRUE
  x$ctrl.mp$ctrl.est@args$comp_b <- TRUE
  ### catch lag fixed
  x$ctrl.mp$ctrl.est@args$catch_lag <- 1
  ### parallelise
  x$genArgs$nblocks <- 10
  x$cut_hist <- FALSE ### retain full history
  return(x)
})
### run MSE with default parameters
res_mp_def <- lapply(input, function(x) {
  do.call(mpDL, x)
})
saveRDS(res_mp_def, file = paste0(path_scn, "1_2_3_1_1_1_1.rds"))
# res_mp_def <- readRDS(paste0(path_scn, "1_2_3_1_1_1_1.rds"))
### optimised parameters
input_opt <- lapply(input, function(x) {
  ### insert optimised parameters
  x$ctrl.mp$ctrl.est@args$idxB_lag     <- par_opt$lag_idx
  x$ctrl.mp$ctrl.est@args$idxB_range_1 <- par_opt$range_idx_1
  x$ctrl.mp$ctrl.est@args$idxB_range_2 <- par_opt$range_idx_2
  x$ctrl.mp$ctrl.est@args$catch_range  <- par_opt$range_catch
  x$ctrl.mp$ctrl.phcr@args$exp_r <- par_opt$exp_r
  x$ctrl.mp$ctrl.phcr@args$exp_f <- par_opt$exp_f
  x$ctrl.mp$ctrl.phcr@args$exp_b <- par_opt$exp_b
  return(x)
})
### run MSE with optimised parameters
res_mp_opt <- lapply(input_opt, function(x) {
  do.call(mpDL, x)
})
saveRDS(res_mp_opt, file = paste0(path_scn, paste0(par_opt, collapse = "_"), 
                                  ".rds"))
# res_mp_opt <- readRDS(paste0(path_scn, paste0(par_opt, collapse = "_"), 
#                                   ".rds"))

### ------------------------------------------------------------------------ ###
### plot time series ####
### ------------------------------------------------------------------------ ###
quantiles <- c(0.05, 0.25, 0.5, 0.75, 0.95)
dfs <- lapply(names(res_mp_def), function(x) {#browser()
  df_opt <- metrics(res_mp_opt[[x]]@stock, 
                    metrics = list(SSB = ssb, F = fbar, Catch = catch))
  df_opt$SSB <- df_opt$SSB/c(input_opt[[x]]$refpts["msy", "ssb"])
  df_opt$F <- df_opt$F/c(input_opt[[x]]$refpts["msy", "harvest"])
  df_opt$Catch <- df_opt$Catch/c(input_opt[[x]]$refpts["msy", "yield"])
  df_opt_p <- lapply(df_opt, function(x) {
    quantile(x, quantiles)
  })
  df_opt_p <- as.data.frame(df_opt_p)
  df_opt_p <- df_opt_p %>% tidyr::spread(key = iter, value = data)
  df_opt_iter <- lapply(df_opt, function(x) {
    iter(x, c(1, 3, 4, 10))
  })
  df_opt_iter <- as.data.frame(df_opt_iter)
  df_opt_iter$iter <- an(ac(df_opt_iter$iter))
  #df_opt <- full_join(df_opt_p, df_opt_iter)
  
  df_def <- metrics(res_mp_def[[x]]@stock, 
                    metrics = list(SSB = ssb, F = fbar, Catch = catch))
  df_def$SSB <- df_def$SSB/c(input_opt[[x]]$refpts["msy", "ssb"])
  df_def$F <- df_def$F/c(input_opt[[x]]$refpts["msy", "harvest"])
  df_def$Catch <- df_def$Catch/c(input_opt[[x]]$refpts["msy", "yield"])
  df_def <- lapply(df_def, iterMedians)
  df_def <- as.data.frame(df_def)
  names(df_def)[7] <- "default"
  
  df_comb <- full_join(df_opt_p, df_def)
  names(df_comb)[9] <- "optimised"
  df_comb <- df_comb %>% gather(key = "parameter", value = `50%`, 
                                "optimised", "default")
  df_comb$stock <- x
  df_opt_iter$stock <- x
  return(list(quantiles = df_comb, iter = df_opt_iter))
})
dfs_quantiles <- do.call(rbind, lapply(dfs, "[[", "quantiles"))
dfs_iter <- do.call(rbind, lapply(dfs, "[[", "iter"))
levels(dfs_quantiles$qname) <- c("SSB/Bmsy", "F/Fmsy", "Catch/MSY")
levels(dfs_iter$qname) <- c("SSB/Bmsy", "F/Fmsy", "Catch/MSY")

ggplot(data = dfs_quantiles,
       aes(x = year - 100, y = `50%`)) +
  geom_ribbon(aes(ymin = `5%`, ymax = `95%`, alpha = "90%")) +
  geom_ribbon(aes(ymin = `25%`, ymax = `75%`, alpha = "50%")) +
  scale_alpha_manual(name = "interval", values = c(0.50, 0.30)) +
  geom_line(aes(linetype = parameter)) +
  scale_linetype_manual(name = "parameters", 
                        values = c(default = "dotted", optimised = "solid")) +
  geom_line(data = dfs_iter,
            aes(x = year - 100, y = data, colour = as.factor(iter)),
            alpha = 0.4, show.legend = FALSE) + 
  facet_grid(qname ~ stock, scales = "free_y") +
  geom_vline(xintercept = 0.5, alpha = 0.5) +
  #geom_hline(yintercept = 1, colour = "red", alpha = 0.3) +
  theme_bw() +
  labs(x = "year", y = "")
ggsave(filename = paste0(path_scn, "timeseries_iter.png"),
       width = 30, height = 17, units = "cm", dpi = 300, type = "cairo")


### ------------------------------------------------------------------------ ###
### plot stats ####
### ------------------------------------------------------------------------ ###

# pos_def <- which(sapply(lapply(lapply(GA_runs, "[[", "pars"), function(x) { 
#   x == c(1, 2, 3, 1, 1, 1, 1)}), sum) == 7)
pos_def <- which(sapply(lapply(lapply(GA_runs, "[[", "pars"), function(x) { 
  x == c(1, 2, 3, 1, 1, 1, 1, 1)}), sum) == 8)
stats <- lapply(GA_runs[c(1, pos_def)], function(x) {
  tmp <- as.data.frame(apply(as.data.frame(x$stats), 2, unlist))
  tmp[2:6, ] <- -tmp[2:6, ]
  tmp$stat <- row.names(tmp)
  tmp <- tmp[, c(6, 1:5)]
  tmp <- gather(tmp, key = "stock", value = "data", 2:6)
})
stats[[1]]$parameter <- "optimised"
stats[[2]]$parameter <- "default"
stats <- do.call(rbind, stats)
### factor levels
stats$stat <- factor(stats$stat,
                     levels = c("obj", "obj_SSB", "obj_F", 
                                "obj_C", "obj_risk", 
                                "obj_icv", "SSB", "SSB_rel",
                                "Fbar", "Fbar_rel",
                                "Catch", "Catch_rel",
                                "ICV",
                                "risk_collapse", "risk_Blim", "risk_halfBmsy",
                                "risk_Bmsy"),
                     labels = c("objective", "SSB objective", "F objective",
                                "Catch objective", "risk objective", 
                                "ICV objective", "SSB", "SSB/Bmsy",
                                "F", "F/Fmsy", 
                                "Catch", "Catch/MSY",
                                "ICV",
                                "collapse risk", "Blim risk", "Bmsy/2 risk",
                                "Bmsy risk"))

### plot
ggplot(data = stats, 
       aes(x = stock, y = data, fill = parameter)) +
  geom_col(position = position_dodge()) +
  facet_wrap(~ stat, scales = "free") +
  theme_bw() +
  labs(y = "")
ggsave(filename = paste0(path_scn, "stats.png"),
       width = 30, height = 17, units = "cm", dpi = 300, type = "cairo")



### ------------------------------------------------------------------------ ###
### multiplier ####
### ------------------------------------------------------------------------ ###
runs <- readRDS("C:/Users/sf02/OneDrive - CEFAS/WKLIFEVII/wklife9/component_r/output/1000_50/ms/trial/gut_whg_bll_lem_ane_jnd_sar_her_san/multiplier_runs.rds")
res <- readRDS("C:/Users/sf02/OneDrive - CEFAS/WKLIFEVII/wklife9/component_r/output/1000_50/ms/trial/gut_whg_bll_lem_ane_jnd_sar_her_san/multiplier_res.rds")

runs[[1]]

plot(sapply(runs, function(x) x$obj) ~ sapply(runs, function(x) x$pars["multiplier"]))
