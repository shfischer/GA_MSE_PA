library(doParallel)
library(GA)
library(ggplot2)
library(tidyr)
library(dplyr)
library(FLCore)

### ------------------------------------------------------------------------ ###
### collate results from HPC ####
### ------------------------------------------------------------------------ ###
stocks <- read.csv("input/stocks.csv", stringsAsFactors = FALSE)
stocks <- stocks[-29, ]
brps <- readRDS("input/OMs/brps.rds")

path_scn <- "output/500_50/SSB_idx_r_only/"
path_scn <- "output/500_50/SSB_idx_r_only_error/"
path_scn <- "output/500_50/SSB_idx_rfb_exp/"
path_scn <- "output/500_50/SSB_idx_rfb_exp_error/"
path_scn <- "output/500_50/SSB_idx_rfb_exp_error_Itrigger/"
### GA results
GA_res <- foreach(stock = stocks$stock, .packages = c("GA"), 
                  .errorhandling = "remove", .final = unlist) %dopar% {
  ret <- list(readRDS(paste0(path_scn, stock, "/res.rds")))
  names(ret) <- stock
  ret
}
saveRDS(GA_res, file = paste0(path_scn, "res_combined.rds"))
GA_res <- c(readRDS(paste0(path_scn, "res_combined_IC.rds")),
            readRDS(paste0(path_scn, "res_combined_UEA.rds")))
GA_res <- GA_res[stocks$stock]
saveRDS(GA_res, file = paste0(path_scn, "res_combined.rds"))

### MSE runs
GA_runs <- foreach(stock = stocks$stock, .packages = c("GA"), 
                  .errorhandling = "remove", 
                  .final = function(x){unlist(x, recursive = FALSE)}) %dopar% {
  ret <- list(readRDS(paste0(path_scn, stock, "/runs.rds")))
  names(ret) <- stock
  ret
}
saveRDS(GA_runs, file = paste0(path_scn, "runs_combined.rds"))
GA_runs <- c(readRDS(paste0(path_scn, "runs_combined_IC.rds")),
            readRDS(paste0(path_scn, "runs_combined_UEA.rds")))
GA_runs <- GA_runs[stocks$stock]
GA_runs <- do.call(rbind, GA_runs)
saveRDS(GA_runs, file = paste0(path_scn, "runs_combined.rds"))

### ------------------------------------------------------------------------ ###
### plot GA results ####
### ------------------------------------------------------------------------ ###
path_scn <- "output/500_50/SSB_idx_rfb_exp/"
path_scn <- "output/500_50/SSB_idx_rfb_exp_error/"
GA_res <- readRDS(file = paste0(path_scn, "res_combined.rds"))
. <- foreach(res_i = GA_res[!sapply(GA_res, is.null)], 
             stock = names(GA_res)) %do% {
  dir.create(paste0(path_scn, "plots/GA_res"), recursive = TRUE)
  png(paste0(path_scn, "plots/GA_res/", stock, ".png"), 
      width = 15, height = 10, units = "cm", 
      res = 300, type = "cairo")
  plot(res_i)
  dev.off()
}

### get some stats about GA runs
GA_stats <- foreach(res_i = GA_res, stock = names(GA_res), 
                    .combine = rbind) %do% {
  cbind(stock = stock, iter = seq(nrow(res_i@summary)), 
        data.frame(res_i@summary))
}
GA_stats %>% 
  gather(key = "key", value = "value", max:min) %>%
  filter(key %in% c("max", "mean", "min")) %>%
  ggplot(aes(x = iter, y = value, linetype = key)) +
  geom_line() + 
  facet_wrap(~ stock) +
  theme_bw()

### ------------------------------------------------------------------------ ###
### optimised: without error ####
### ------------------------------------------------------------------------ ###
path_scn <- "output/500_50/SSB_idx_rfb_exp/"

### run optimised solution
library(doParallel)
cl <- makeCluster(10)
registerDoParallel(cl)
cl_length <- length(cl)
req_pckgs <- c("FLCore", "FLash", "mseDL", "GA", "doParallel", "doRNG")
for (i in req_pckgs) library(package = i, character.only = TRUE)
source("funs.R"); source("GA_funs.R")
. <- foreach(i = seq(cl_length)) %dopar% {
  for (i in req_pckgs) library(package = i, character.only = TRUE)
  source("funs.R"); source("GA_funs.R")
}
stocks <- read.csv("input/stocks.csv", stringsAsFactors = FALSE)
stocks <- stocks[-29, ]
smry <- read.csv(paste0(path_scn, "summary.csv"), stringsAsFactors = FALSE)

. <- foreach(stock = stocks$stock, .errorhandling = "remove") %do% {
  cat(stock, "\n"); flush.console()
  input <- readRDS(paste0("input/500_50/OM_2_mp_input/", stock, ".rds"))
  input$oem@args$idx_dev <- FALSE
  input$oem@args$ssb <- TRUE
  input$oem@args$lngth <- TRUE
  input$oem@args$lngth_dev <- FALSE
  ### catch rule components
  input$ctrl.mp$ctrl.est@args$comp_r <- TRUE
  input$ctrl.mp$ctrl.est@args$comp_f <- TRUE
  input$ctrl.mp$ctrl.est@args$comp_b <- TRUE
  ### catch lag fixed
  input$ctrl.mp$ctrl.est@args$catch_lag <- 1
  ### optimised parameters
  params <- smry[smry$stock == stock, ]
  input$ctrl.mp$ctrl.est@args$idxB_lag     <- params$lag_idx
  input$ctrl.mp$ctrl.est@args$idxB_range_1 <- params$range_idx_1
  input$ctrl.mp$ctrl.est@args$idxB_range_2 <- params$range_idx_2
  input$ctrl.mp$ctrl.est@args$catch_range  <- params$range_catch
  input$ctrl.mp$ctrl.phcr@args$exp_r <- params$exp_r
  input$ctrl.mp$ctrl.phcr@args$exp_f <- params$exp_f
  input$ctrl.mp$ctrl.phcr@args$exp_b <- params$exp_b
  ### parallelise
  input$genArgs$nblocks <- 10
  input$cut_hist <- FALSE
  invisible(gc()); invisible(foreach(i = seq(cl_length)) %dopar% gc())
  ### run
  res_mp <- do.call(mpDL, input)
  ### save
  dir.create(paste0(path_scn, stock), recursive = TRUE)
  saveRDS(res_mp, paste0(path_scn, stock, "/res_", params$file))
  ### plot
  stk <- res_mp@stock
  stk <- window(stk, start = -25)
  stk[, ac(-25:50)] <- stk[, ac(75:150)]
  stk <- window(stk, end = 50)
  p <- plot(stk, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), iter = 11:15) +
    ylim(0, NA) + labs(x = "year") +
    geom_hline(data = data.frame(qname = "SSB", 
                                 data = c(input$refpts["msy", "ssb"])),
               aes(yintercept = data), linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = 0, alpha = 0.5) +
    scale_x_continuous(expand = c(0, 0))
  ggsave(filename = paste0(path_scn, stock, "/res_", 
                           gsub(x = params$file, pattern = "rds", 
                                replacement = "png")),
         width = 30, height = 20, units = "cm", dpi = 300, type = "cairo", 
         plot = p)
  ### run with default parameters
  input$ctrl.mp$ctrl.est@args$idxB_lag     <- 1
  input$ctrl.mp$ctrl.est@args$idxB_range_1 <- 2
  input$ctrl.mp$ctrl.est@args$idxB_range_2 <- 3
  input$ctrl.mp$ctrl.est@args$catch_range  <- 1
  input$ctrl.mp$ctrl.phcr@args$exp_r <- 1
  input$ctrl.mp$ctrl.phcr@args$exp_f <- 1
  input$ctrl.mp$ctrl.phcr@args$exp_b <- 1
  res_mp <- do.call(mpDL, input)
  saveRDS(res_mp, paste0(path_scn, stock, "/res_1_2_3_1_1_1_1.rds"))
  stk <- res_mp@stock
  stk <- window(stk, start = -25)
  stk[, ac(-25:50)] <- stk[, ac(75:150)]
  stk <- window(stk, end = 50)
  p <- plot(stk, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), iter = 11:15) +
    ylim(0, NA) + labs(x = "year") +
    geom_hline(data = data.frame(qname = "SSB", 
                                 data = c(input$refpts["msy", "ssb"])),
               aes(yintercept = data), linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = 0, alpha = 0.5) +
    scale_x_continuous(expand = c(0, 0))
  ggsave(filename = paste0(path_scn, stock, "/res_1_2_3_1_1_1_1.png"),
         width = 30, height = 20, units = "cm", dpi = 300, type = "cairo",
         plot = p)
  rm(input, res_mp)
  invisible(gc()); invisible(foreach(i = seq(cl_length)) %dopar% gc())
  return(NULL)
}

### get stats for default scenario
runs_combined <- readRDS(paste0(path_scn, "runs_combined.rds"))
runs <- runs_combined[runs_combined$file == "1_2_3_1_1_1_1.rds", ]
stats_default <- merge(stocks[, c("k", "stock")], runs)
stats_default <- stats_default[match(stocks$stock, stats_default$stock), ]
write.csv(stats_default, file = paste0(path_scn, "stats_default.csv"),
          row.names = FALSE)
### format optimised solutions
smry <- read.csv(paste0(path_scn, "summary.csv"), stringsAsFactors = FALSE)
stats_opt <- merge(stocks[, c("k", "stock")], smry)
stats_opt <- stats_opt[match(stocks$stock, stats_opt$stock), ]
write.csv(stats_opt, file = paste0(path_scn, "stats_opt.csv"),
          row.names = FALSE)

stats_cf <- stats_opt[, c("k", "stock", "lag_idx", "range_idx_1", "range_idx_2",
                          "range_catch", "exp_r", "exp_f", "exp_b")]
stats_cf$fitness_default <- stats_default$fitness
stats_cf$fitness_opt <- stats_opt$fitness
stats_cf$fitness_improvement <- (1 - stats_cf$fitness_opt/stats_cf$fitness_default)*100
stats_cf$dev_default <- stats_default$Bmsy_dev
stats_cf$dev_opt <- stats_opt$Bmsy_dev
stats_cf$dev_improvement <- (1 - stats_cf$dev_opt/stats_cf$dev_default)*100
write.csv(stats_cf, file = paste0(path_scn, "stats_cf.csv"),
          row.names = FALSE)

### ------------------------------------------------------------------------ ###
### optimised: with uncertainty ####
### ------------------------------------------------------------------------ ###
path_scn <- "output/500_50/SSB_idx_rfb_exp_error_Itrigger/"

### run optimised solution
library(doParallel)
cl <- makeCluster(20)
registerDoParallel(cl)
cl_length <- length(cl)
req_pckgs <- c("FLCore", "FLash", "mseDL", "GA", "doParallel", "doRNG")
for (i in req_pckgs) library(package = i, character.only = TRUE)
source("funs.R"); source("GA_funs.R")
. <- foreach(i = seq(cl_length)) %dopar% {
  for (i in req_pckgs) library(package = i, character.only = TRUE)
  source("funs.R"); source("GA_funs.R")
}
stocks <- read.csv("input/stocks.csv", stringsAsFactors = FALSE)
stocks <- stocks[-29, ]
smry <- read.csv(paste0(path_scn, "summary.csv"), stringsAsFactors = FALSE)

. <- foreach(stock = stocks$stock, .errorhandling = "remove") %do% {
  cat(stock, "\n"); flush.console()
  input <- readRDS(paste0("input/500_50/OM_2_mp_input/", stock, ".rds"))
  input$oem@args$idx_dev <- TRUE
  input$oem@args$ssb <- TRUE
  input$oem@args$lngth <- TRUE
  input$oem@args$lngth_dev <- TRUE
  ### catch rule components
  input$ctrl.mp$ctrl.est@args$comp_r <- TRUE
  input$ctrl.mp$ctrl.est@args$comp_f <- TRUE
  input$ctrl.mp$ctrl.est@args$comp_b <- TRUE
  ### catch lag fixed
  input$ctrl.mp$ctrl.est@args$catch_lag <- 1
  ### optimised parameters
  params <- smry[smry$stock == stock, ]
  input$ctrl.mp$ctrl.est@args$idxB_lag     <- params$lag_idx
  input$ctrl.mp$ctrl.est@args$idxB_range_1 <- params$range_idx_1
  input$ctrl.mp$ctrl.est@args$idxB_range_2 <- params$range_idx_2
  input$ctrl.mp$ctrl.est@args$catch_range  <- params$range_catch
  input$ctrl.mp$ctrl.phcr@args$exp_r <- params$exp_r
  input$ctrl.mp$ctrl.phcr@args$exp_f <- params$exp_f
  input$ctrl.mp$ctrl.phcr@args$exp_b <- params$exp_b
  ### parallelise
  input$genArgs$nblocks <- cl_length
  input$cut_hist <- FALSE
  ### use Itrigger #############################################################
  input$ctrl.mp$ctrl.est@args$I_trigger[] <- c(input$I_loss$SSB_idx_dev) * 1.4
  ###
  invisible(gc()); invisible(foreach(i = seq(cl_length)) %dopar% gc())
  ### run
  res_mp <- do.call(mpDL, input)
  ### save
  dir.create(paste0(path_scn, stock), recursive = TRUE)
  saveRDS(res_mp, paste0(path_scn, stock, "/res_", params$file))
  ### plot
  stk <- res_mp@stock
  stk <- window(stk, start = -25)
  stk[, ac(-25:50)] <- stk[, ac(75:150)]
  stk <- window(stk, end = 50)
  p <- plot(stk, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), iter = 11:15) +
    ylim(0, NA) + labs(x = "year") +
    geom_hline(data = data.frame(qname = "SSB", 
                                 data = c(input$refpts["msy", "ssb"])),
               aes(yintercept = data), linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = 0, alpha = 0.5) +
    scale_x_continuous(expand = c(0, 0))
  ggsave(filename = paste0(path_scn, stock, "/res_", 
                           gsub(x = params$file, pattern = "rds", 
                                replacement = "png")),
         width = 30, height = 20, units = "cm", dpi = 300, type = "cairo", 
         plot = p)
  ### run with default parameters
  input$ctrl.mp$ctrl.est@args$idxB_lag     <- 1
  input$ctrl.mp$ctrl.est@args$idxB_range_1 <- 2
  input$ctrl.mp$ctrl.est@args$idxB_range_2 <- 3
  input$ctrl.mp$ctrl.est@args$catch_range  <- 1
  input$ctrl.mp$ctrl.phcr@args$exp_r <- 1
  input$ctrl.mp$ctrl.phcr@args$exp_f <- 1
  input$ctrl.mp$ctrl.phcr@args$exp_b <- 1
  res_mp <- do.call(mpDL, input)
  saveRDS(res_mp, paste0(path_scn, stock, "/res_1_2_3_1_1_1_1.rds"))
  stk <- res_mp@stock
  stk <- window(stk, start = -25)
  stk[, ac(-25:50)] <- stk[, ac(75:150)]
  stk <- window(stk, end = 50)
  p <- plot(stk, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), iter = 11:15) +
    ylim(0, NA) + labs(x = "year") +
    geom_hline(data = data.frame(qname = "SSB", 
                                 data = c(input$refpts["msy", "ssb"])),
               aes(yintercept = data), linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = 0, alpha = 0.5) +
    scale_x_continuous(expand = c(0, 0))
  ggsave(filename = paste0(path_scn, stock, "/res_1_2_3_1_1_1_1.png"),
         width = 30, height = 20, units = "cm", dpi = 300, type = "cairo",
         plot = p)
  rm(input, res_mp)
  invisible(gc()); invisible(foreach(i = seq(cl_length)) %dopar% gc())
  return(NULL)
}

### get stats for default scenario
runs_combined <- readRDS(paste0(path_scn, "runs_combined.rds"))
runs <- runs_combined[runs_combined$file == "1_2_3_1_1_1_1.rds", ]
stats_default <- merge(stocks[, c("k", "stock")], runs)
stats_default <- stats_default[match(stocks$stock, stats_default$stock), ]
write.csv(stats_default, file = paste0(path_scn, "stats_default.csv"),
          row.names = FALSE)
### format optimised solutions
smry <- read.csv(paste0(path_scn, "summary.csv"), stringsAsFactors = FALSE)
stats_opt <- merge(stocks[, c("k", "stock")], smry)
stats_opt <- stats_opt[match(stocks$stock, stats_opt$stock), ]
write.csv(stats_opt, file = paste0(path_scn, "stats_opt.csv"),
          row.names = FALSE)

stats_cf <- stats_opt[, c("k", "stock", "lag_idx", "range_idx_1", "range_idx_2",
                          "range_catch", "exp_r", "exp_f", "exp_b")]
stats_cf$fitness_default <- stats_default$fitness
stats_cf$fitness_opt <- stats_opt$fitness
stats_cf$fitness_improvement <- (1 - stats_cf$fitness_opt/stats_cf$fitness_default)*100
stats_cf$dev_default <- stats_default$Bmsy_dev
stats_cf$dev_opt <- stats_opt$Bmsy_dev
stats_cf$dev_improvement <- (1 - stats_cf$dev_opt/stats_cf$dev_default)*100
write.csv(stats_cf, file = paste0(path_scn, "stats_cf.csv"),
          row.names = FALSE)


### ------------------------------------------------------------------------ ###
### plot SSBs ####
### ------------------------------------------------------------------------ ###

# path_scn <- "output/500_50/SSB_idx_r_only/"
# path_scn <- "output/500_50/SSB_idx_r_only_error/"
# path_scn <- "output/500_50/SSB_idx_rfb_exp/"
# path_scn <- "output/500_50/SSB_idx_rfb_exp_error/"
# path_scn <- "output/500_50/SSB_idx_rfb_exp_error_Itrigger/"
smry <- read.csv(paste0(path_scn, "summary.csv"), stringsAsFactors = FALSE)

# stocks_subset <- stocks[-29, c("stock", "k")]
# stocks_subset$name <- with(stocks_subset, paste0("k=", k, ", ", stock))
# ssbs <- foreach(stock = stocks_subset$stock, .errorhandling = "pass") %do% {
#   refpts_i <- refpts(brps[[stock]])
#   file_def <- "res_1_2_3_1_1_1_1.rds"
#   file_opt <- paste0("res_", smry$file[smry$stock == stock])
#   res_mp_def <- readRDS(paste0(path_scn, stock, "/", file_def))
#   res_mp_opt <- readRDS(paste0(path_scn, stock, "/", file_opt))
#   ssb_rel_def <- ssb(stock(res_mp_def))/c(refpts_i["msy","ssb"])
#   ssb_rel_opt <- ssb(stock(res_mp_opt))/c(refpts_i["msy","ssb"])
#   return(list(def = ssb_rel_def, opt = ssb_rel_opt))
# }
# ssbs_def <- lapply(ssbs, "[[", "def")
# ssbs_opt <- lapply(ssbs, "[[", "opt")
# names(ssbs_def) <- names(ssbs_opt) <- stocks_subset$name
# ssbs_def <- ssbs_def[sapply(ssbs_def, is, "FLQuant")]
# ssbs_opt <- ssbs_opt[sapply(ssbs_opt, is, "FLQuant")]
# saveRDS(ssbs_def, file = paste0(path_scn, "SSBs_default.rds"))
# saveRDS(ssbs_opt, file = paste0(path_scn, "SSBs_opt.rds"))
ssbs_def <- readRDS(paste0(path_scn, "SSBs_default.rds"))
ssbs_opt<- readRDS(paste0(path_scn, "SSBs_opt.rds"))

df_quant_def <- as.data.frame(FLQuants(ssbs_def)) %>%
  mutate(year = year - 100) %>%
  group_by(qname, year) %>%
  summarise(q5 = quantile(data, 0.05), q25 = quantile(data, 0.25),
            q50 = quantile(data, 0.5), q75 = quantile(data, 0.75),
            q95 = quantile(data, 0.95))
df_def <- as.data.frame(FLQuants(ssbs_def)) %>%
  mutate(year = year - 100) %>%
  group_by(qname, year) %>%
  filter(iter %in% 1:500)
df_quant_opt <- as.data.frame(FLQuants(ssbs_opt)) %>%
  mutate(year = year - 100) %>%
  group_by(qname, year) %>%
  summarise(q5 = quantile(data, 0.05), q25 = quantile(data, 0.25),
            q50 = quantile(data, 0.5), q75 = quantile(data, 0.75),
            q95 = quantile(data, 0.95))
df_opt <- as.data.frame(FLQuants(ssbs_opt)) %>%
  mutate(year = year - 100) %>%
  group_by(qname, year) %>%
  filter(iter %in% 1:500)
ggplot(data = df_quant_def,
       aes(x = year, y = q50)) +
  geom_ribbon(aes(ymin = q5, ymax = q95), colour = "grey", alpha = 0.3,
              linetype = 0) +
  geom_ribbon(aes(ymin = q25, ymax = q75), colour = "grey", alpha = 0.5,
              linetype = 0) +
  geom_line() +
  geom_line(data = df_def %>%
              filter(iter %in% c(1, 4, 11, 22, 64)), 
            aes(x = year, y = data, colour = iter), 
            show.legend = FALSE) +
  theme_bw() +
  labs(x = "year", y = expression(SSB/B[MSY])) +
  facet_wrap(~ qname) +
  ylim(0, NA) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_vline(xintercept = 0.5, alpha = 0.5)
ggsave(filename = paste0(path_scn, "SSB_default_iter.png"),
       width = 30, height = 20, units = "cm", dpi = 200, type = "cairo")
ggplot(data = df_quant_opt,
       aes(x = year, y = q50)) +
  geom_ribbon(aes(ymin = q5, ymax = q95), colour = "grey", alpha = 0.3,
              linetype = 0) +
  geom_ribbon(aes(ymin = q25, ymax = q75), colour = "grey", alpha = 0.5,
              linetype = 0) +
  geom_line() +
  geom_line(data = df_opt %>%
              filter(iter %in% c(1, 4, 11, 22, 64)),
            aes(x = year, y = data, colour = iter), 
            show.legend = FALSE) +
  theme_bw() +
  labs(x = "year", y = expression(SSB/B[MSY])) +
  facet_wrap(~ qname) +
  ylim(0, NA) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_vline(xintercept = 0.5, alpha = 0.5)
ggsave(filename = paste0(path_scn, "SSB_optimised_iter.png"),
       width = 30, height = 20, units = "cm", dpi = 200, type = "cairo")

### ------------------------------------------------------------------------ ###
### update stock plots ####
### ------------------------------------------------------------------------ ###
# path_scn <- "output/500_50/SSB_idx_r_only/"
# path_scn <- "output/500_50/SSB_idx_r_only_error/"
# path_scn <- "output/500_50/SSB_idx_rfb_exp/"
# path_scn <- "output/500_50/SSB_idx_rfb_exp_error/"
# path_scn <- "output/500_50/SSB_idx_rfb_exp_error_Itrigger/"
smry <- read.csv(paste0(path_scn, "summary.csv"), stringsAsFactors = FALSE)
. <- foreach(stock = stocks$stock, refpts = lapply(brps[stocks$stock], refpts),
             .errorhandling = "remove") %dopar% {
  cat(stock, "\n"); flush.console()
  file_opt <- smry[smry$stock == stock, "file"]
  ### optimised solution
  res_mp <- readRDS(paste0(path_scn, stock, "/res_", file_opt))
  stk <- res_mp@stock
  stk <- window(stk, start = -25)
  stk[, ac(-25:50)] <- stk[, ac(75:150)]
  stk <- window(stk, end = 50)
  p <- plot(stk, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), 
            iter = c(1, 4, 11, 22, 64)) +
    ylim(0, NA) + labs(x = "year") +
    geom_hline(data = data.frame(qname = c("SSB", "Catch", "F", "Rec"), 
                                 data = c(refpts["msy", "ssb"],
                                          refpts["msy", "yield"],
                                          refpts["msy", "harvest"],
                                          refpts["msy", "rec"])),
               aes(yintercept = data), linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = 0, alpha = 0.5) +
    scale_x_continuous(expand = c(0, 0))
  ggsave(filename = paste0(path_scn, stock, "/res_", 
                           gsub(x = params$file, pattern = ".rds", 
                                replacement = "_new_iter.png")),
         width = 30, height = 20, units = "cm", dpi = 300, type = "cairo", 
         plot = p)
  ### default
  res_mp <- readRDS(paste0(path_scn, stock, "/res_1_2_3_1_1_1_1.rds"))
  stk <- res_mp@stock
  stk <- window(stk, start = -25)
  stk[, ac(-25:50)] <- stk[, ac(75:150)]
  stk <- window(stk, end = 50)
  p <- plot(stk, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), 
            iter = c(1, 4, 11, 22, 64)) +
    ylim(0, NA) + labs(x = "year") +
    geom_hline(data = data.frame(qname = c("SSB", "Catch", "F", "Rec"), 
                                 data = c(refpts["msy", "ssb"],
                                          refpts["msy", "yield"],
                                          refpts["msy", "harvest"],
                                          refpts["msy", "rec"])),
               aes(yintercept = data), linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = 0, alpha = 0.5) +
    scale_x_continuous(expand = c(0, 0))
  ggsave(filename = paste0(path_scn, stock, "/res_1_2_3_1_1_1_1_new_iter.png"),
         width = 30, height = 20, units = "cm", dpi = 300, type = "cairo",
         plot = p)
}


### ------------------------------------------------------------------------ ###
### r only: optimised: without/with error ####
### ------------------------------------------------------------------------ ###
path_scn <- "output/500_50/SSB_idx_r_only_error/"
smry <- read.csv(paste0(path_scn, "summary.csv"), stringsAsFactors = FALSE)

stocks_subset <- stocks$stock
. <- foreach(stock = stocks_subset, .errorhandling = "remove") %do% {
  cat(stock, "\n"); flush.console()
  input <- readRDS(paste0("input/500_50/OM_2_mp_input/", stock, ".rds"))
  input$oem@args$idx_dev <- TRUE
  input$oem@args$ssb <- TRUE
  input$oem@args$lngth <- FALSE
  input$oem@args$lngth_dev <- FALSE
  ### catch rule components
  input$ctrl.mp$ctrl.est@args$comp_r <- TRUE
  input$ctrl.mp$ctrl.est@args$comp_f <- FALSE
  input$ctrl.mp$ctrl.est@args$comp_b <- FALSE
  ### catch lag fixed
  input$ctrl.mp$ctrl.est@args$catch_lag <- 1
  ### optimised parameters
  params <- smry[smry$stock == stock, ]
  input$ctrl.mp$ctrl.est@args$idxB_lag     <- params$lag_idx
  input$ctrl.mp$ctrl.est@args$idxB_range_1 <- params$range_idx_1
  input$ctrl.mp$ctrl.est@args$idxB_range_2 <- params$range_idx_2
  input$ctrl.mp$ctrl.est@args$catch_range  <- params$range_catch
  ### parallelise
  input$genArgs$nblocks <- 10
  input$cut_hist <- FALSE
  invisible(gc()); invisible(foreach(i = seq(cl_length)) %dopar% gc())
  ### run
  res_mp <- do.call(mpDL, input)
  ### save
  dir.create(paste0(path_scn, stock), recursive = TRUE)
  saveRDS(res_mp, paste0(path_scn, stock, "/res_", params$file))
  ### plot
  stk <- res_mp@stock
  stk <- window(stk, start = -25)
  stk[, ac(-25:50)] <- stk[, ac(75:150)]
  stk <- window(stk, end = 50)
  p <- plot(stk, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), iter = 11:15) +
    ylim(0, NA) + labs(x = "year") +
    geom_hline(data = data.frame(qname = "SSB", 
                                 data = c(input$refpts["msy", "ssb"])),
               aes(yintercept = data), linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = 0, alpha = 0.5) +
    scale_x_continuous(expand = c(0, 0))
  ggsave(filename = paste0(path_scn, stock, "/res_", 
                           gsub(x = params$file, pattern = "rds", 
                                replacement = "png")),
         width = 30, height = 20, units = "cm", dpi = 300, type = "cairo", 
         plot = p)
  ### run with default parameters
  input$ctrl.mp$ctrl.est@args$idxB_lag     <- 1
  input$ctrl.mp$ctrl.est@args$idxB_range_1 <- 2
  input$ctrl.mp$ctrl.est@args$idxB_range_2 <- 3
  input$ctrl.mp$ctrl.est@args$catch_range  <- 1
  res_mp <- do.call(mpDL, input)
  saveRDS(res_mp, paste0(path_scn, stock, "/res_1_2_3_1_1_1_1.rds"))
  stk <- res_mp@stock
  stk <- window(stk, start = -25)
  stk[, ac(-25:50)] <- stk[, ac(75:150)]
  stk <- window(stk, end = 50)
  p <- plot(stk, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), iter = 11:15) +
    ylim(0, NA) + labs(x = "year") +
    geom_hline(data = data.frame(qname = "SSB", 
                                 data = c(input$refpts["msy", "ssb"])),
               aes(yintercept = data), linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = 0, alpha = 0.5) +
    scale_x_continuous(expand = c(0, 0))
  ggsave(filename = paste0(path_scn, stock, "/res_1_2_3_1_1_1_1.png"),
         width = 30, height = 20, units = "cm", dpi = 300, type = "cairo",
         plot = p)
  rm(input, res_mp)
  invisible(gc()); invisible(foreach(i = seq(cl_length)) %dopar% gc())
  return(NULL)
}

### get stats for default scenario
runs_combined <- readRDS(paste0(path_scn, "runs_combined.rds"))
runs_combined <- do.call(rbind, runs_combined)
runs <- runs_combined[runs_combined$file == "1_2_3_1.rds", ]
stats_default <- merge(stocks[, c("k", "stock")], runs)
stats_default <- stats_default[match(stocks$stock, stats_default$stock), ]
write.csv(stats_default, file = paste0(path_scn, "stats_default.csv"),
          row.names = FALSE)
### format optimised solutions
smry <- read.csv(paste0(path_scn, "summary.csv"), stringsAsFactors = FALSE)
stats_opt <- merge(stocks[, c("k", "stock")], smry)
stats_opt <- stats_opt[match(stocks$stock, stats_opt$stock), ]
write.csv(stats_opt, file = paste0(path_scn, "stats_opt.csv"),
          row.names = FALSE)

stats_cf <- stats_opt[, c("k", "stock", "lag_idx", "range_idx_1", "range_idx_2",
                          "range_catch")]
stats_cf$fitness_default <- stats_default$fitness
stats_cf$fitness_opt <- stats_opt$fitness
stats_cf$fitness_improvement <- (1 - stats_cf$fitness_opt/stats_cf$fitness_default)*100
stats_cf$dev_default <- stats_default$Bmsy_dev
stats_cf$dev_opt <- stats_opt$Bmsy_dev
stats_cf$dev_improvement <- (1 - stats_cf$dev_opt/stats_cf$dev_default)*100
write.csv(stats_cf, file = paste0(path_scn, "stats_cf.csv"),
          row.names = FALSE)


