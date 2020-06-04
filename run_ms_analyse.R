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
### new runs, new objective function, 04/20 ####
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### collate results - optimised parameters ####
### ------------------------------------------------------------------------ ###
scenario <- "trial"
fhist <- "one-way"#"random"
n_iter <- 500
n_yrs <- 50

stocks <- read.csv("input/stocks.csv", stringsAsFactors = FALSE)

stocks_subset <- stocks$stock#[21:27]#"pol"
names(stocks_subset) <- stocks_subset

### load GA results
res_lst <- lapply(stocks_subset, function(x) {
  file <- paste0("output/", n_iter, "_", n_yrs, "/ms/", scenario, "/", 
                        fhist, "/", x, "/lag_idx-range_idx_1-range_idx_2",
                        "-exp_r-exp_f-exp_b-interval-multiplier--obj_SSB_C_",
                        "risk_ICV_res.rds")
  if (file.exists(file))
    tmp <- readRDS(file)
  else
    NULL
})
res_lst <- res_lst[!sapply(res_lst, is.null)]

res_par <- lapply(res_lst, function(x) {
  tmp <- x@solution[1,]
  tmp[c(1:4, 8)] <- round(tmp[c(1:4, 8)])
  tmp[5:7] <- round(tmp[5:7], 1)
  tmp[9] <- round(tmp[9], 2)
  return(tmp)
})
saveRDS(res_par, paste0("output/", n_iter, "_", n_yrs, "/ms/trial/all_stocks_",
                        fhist, "_opt_pars.rds"))



### objective function trials for pollack
fhist <- "random"#"one-way"
fs <- list.files(path = paste0("output/500_50/ms/trial/", fhist, "/pol/"),
                 pattern = "*_res.rds")
fs <- fs[grep(x = fs, pattern = "lag_idx-range_idx_1-range_idx_2-exp_r-exp_f-exp_b-interval-multiplier--obj_")]
trials <- data.frame(file = fs)
trials$obj_fun <- gsub(x = trials$file, pattern = "lag_idx-range_idx_1-range_idx_2-exp_r-exp_f-exp_b-interval-multiplier--obj_|_res\\.rds", replacement = "")
trials$obj_fun
res_lst <- lapply(trials$file, function(x) {
  readRDS(paste0("output/500_50/ms/trial/", fhist, "/pol/", x))
})
res_par <- lapply(res_lst, function(x) {
  tmp <- x@solution[1,]
  tmp[c(1:4, 8)] <- round(tmp[c(1:4, 8)])
  tmp[5:7] <- round(tmp[5:7], 1)
  tmp[9] <- round(tmp[9], 2)
  return(tmp)
})
names(res_par) <- trials$obj_fun
saveRDS(res_par, paste0("output/500_50/ms/trial/pol_", fhist, 
                        "_obj_funs_res.rds"))

### "multi-species" runs
ms_groups <- list(low = 1:12, medium = 13:20, high = 21:29)
ms_groups <- sapply(ms_groups, function(x) {paste0(stocks$stock[x], collapse = "_")})
names(ms_groups) <- unlist(ms_groups)

res_lst <- lapply(ms_groups, function(x) {
  readRDS(paste0("output/500_50/ms/trial/one-way/", x, "/",
                 "lag_idx-range_idx_1-range_idx_2-exp_r-exp_f-exp_b-interval-", 
                 "multiplier--obj_SSB_C_risk_ICV_res.rds"))
})
res_par <- lapply(res_lst, function(x) {
  tmp <- x@solution[1,]
  tmp[c(1:4, 8)] <- round(tmp[c(1:4, 8)])
  tmp[5:7] <- round(tmp[5:7], 1)
  tmp[9] <- round(tmp[9], 2)
  return(tmp)
})
saveRDS(res_par, paste0("output/500_50/ms/trial/ms_res.rds"))

### ------------------------------------------------------------------------ ###
### recreate MSE runs ####
### ------------------------------------------------------------------------ ###
### default, optimised

### run optimised solution
library(doParallel)
n_cores <- ifelse(parallel::detectCores() > 20, 20, 10)
cl <- makeCluster(n_cores)
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
# stocks <- read.csv("input/stocks.csv", stringsAsFactors = FALSE)[-29, ]
# stocks <- names(stocks_subset)
# names(stocks) <- stocks
stocks <- names(res_lst)
names(stocks) <- stocks
#stocks <- rev(stocks)


fhist <- "one-way"#"random"
scenario <- "trial"
n_iter <- 500
n_yrs <- 50


### by stock
for (stock in stocks) {
  
  rm(res_mp_def, res_mp_zero, res_mp_opt)
  
  input <- lapply(stock, function(x) {
    readRDS(paste0("input/", n_iter, "_", n_yrs, "/OM_2_mp_input/", fhist, "/", x,
                  ".rds"))
  })

  ### prepare input object(s) for MSE
  input <- lapply(input, function(x) {
    ### OEM: activate uncertainty
    x$oem@args$idx_dev <- TRUE
    x$oem@args$ssb <- FALSE
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
    x$genArgs$nblocks <- n_cores
    x$cut_hist <- FALSE ### retain full history
    return(x)
  })
  ### run MSE with default parameters
  res_mp_def <- lapply(input, function(x) {
    do.call(mpDL, x)
  })
  path_out <- paste0("output/", n_iter, "_", n_yrs, "/ms/", scenario, "/", fhist, 
                     "/")
  saveRDS(res_mp_def, file = paste0(path_out, stock,
                                    "/mp_1_2_3_1_1_1_1_2_1.rds"))
  
    ### zero catch
    input_zero <- lapply(input, function(x) {
      x$ctrl.mp$ctrl.phcr@args$multiplier <- 0
      return(x)
    })
    res_mp_zero <- lapply(input_zero, function(x) {
      do.call(mpDL, x)
    })
    saveRDS(res_mp_zero, file = paste0(path_out, stock,
                                      "/mp_1_2_3_1_1_1_1_2_0.rds"))
    
    ### optimised parameters
    params <- res_par[[stock]]
    
    input_opt <- lapply(input, function(x) {
      ### insert optimised parameters
      x$ctrl.mp$ctrl.est@args$idxB_lag     <- params[1]
      x$ctrl.mp$ctrl.est@args$idxB_range_1 <- params[2]
      x$ctrl.mp$ctrl.est@args$idxB_range_2 <- params[3]
      x$ctrl.mp$ctrl.est@args$catch_range  <- params[4]
      x$ctrl.mp$ctrl.phcr@args$exp_r <- params[5]
      x$ctrl.mp$ctrl.phcr@args$exp_f <- params[6]
      x$ctrl.mp$ctrl.phcr@args$exp_b <- params[7]
      x$ctrl.mp$ctrl.hcr@args$interval <- params[8]
      x$ctrl.mp$ctrl.phcr@args$multiplier <- params[9]
      return(x)
    })
    ### run MSE with optimised parameters
    res_mp_opt <- lapply(input_opt, function(x) {
      do.call(mpDL, x)
    })
    saveRDS(res_mp_opt, file = paste0(path_out, stock, "/mp_", 
                                      paste0(params, collapse = "_"), 
                                      ".rds"))
  
}
### pollack objective function trials

for (obj_fun in names(res_par)) {
  
  rm(res_mp_opt)
  
  input <- lapply("pol", function(x) {
    readRDS(paste0("input/", n_iter, "_", n_yrs, "/OM_2_mp_input/", fhist, "/",
                   x, ".rds"))
  })
  names(input) <- "pol"

  ### prepare input object(s) for MSE
  input <- lapply(input, function(x) {
    ### OEM: activate uncertainty
    x$oem@args$idx_dev <- TRUE
    x$oem@args$ssb <- FALSE
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
    x$genArgs$nblocks <- 20
    x$cut_hist <- FALSE ### retain full history
    return(x)
  })

  ### optimised parameters
  params <- res_par[[obj_fun]]
  
  input_opt <- lapply(input, function(x) {
    ### insert optimised parameters
    x$ctrl.mp$ctrl.est@args$idxB_lag     <- params[1]
    x$ctrl.mp$ctrl.est@args$idxB_range_1 <- params[2]
    x$ctrl.mp$ctrl.est@args$idxB_range_2 <- params[3]
    x$ctrl.mp$ctrl.est@args$catch_range  <- params[4]
    x$ctrl.mp$ctrl.phcr@args$exp_r <- params[5]
    x$ctrl.mp$ctrl.phcr@args$exp_f <- params[6]
    x$ctrl.mp$ctrl.phcr@args$exp_b <- params[7]
    x$ctrl.mp$ctrl.hcr@args$interval <- params[8]
    x$ctrl.mp$ctrl.phcr@args$multiplier <- params[9]
    return(x)
  })
  ### run MSE with optimised parameters
  res_mp_opt <- lapply(input_opt, function(x) {
    do.call(mpDL, x)
  })
  path_out <- paste0("output/", n_iter, "_", n_yrs, "/ms/", scenario, "/", 
                     fhist, "/")
  saveRDS(res_mp_opt, file = paste0(path_out, "pol", "/mp_", 
                                    paste0(params, collapse = "_"), 
                                    ".rds"))
  
}

### ------------------------------------------------------------------------ ###
### pollack: objective functions and fishing histories - plots ####
### ------------------------------------------------------------------------ ###

### optimised parameters for objective function trials (one-way)
pol_par_ow_obj <- readRDS("output/500_50/ms/trial/pol_one-way_obj_funs_res.rds")
pol_par_rnd_obj <- readRDS("output/500_50/ms/trial/pol_random_obj_funs_res.rds")

### format 
pol_pars <- rbind(t(as.data.frame(pol_par_ow_obj)),
                  t(as.data.frame(pol_par_rnd_obj)))
pol_pars <- as.data.frame(pol_pars)
pol_pars$obj_def <- row.names(pol_pars)
pol_pars$fhist <- c(rep("one-way", length(pol_par_ow_obj)),
                    rep("random", length(pol_par_rnd_obj)))
pol_pars$optimised <- TRUE
pol_pars <- do.call(rbind, list(pol_pars, 
                                c(1, 2, 3, 1, 1, 1, 1, 2, 1, "", "one-way", 
                                  FALSE),
                                c(1, 2, 3, 1, 1, 1, 1, 2, 1, "", "random", 
                                  FALSE)))
pol_pars$seq <- seq(nrow(pol_pars))
pol_pars <- pol_pars %>% 
  group_by(seq) %>%
  mutate(obj_def = gsub(x = obj_def, pattern = "\\.[0-9]{1,}", 
                        replacement = ""),
         file = paste0(c(lag_idx, range_idx_1, range_idx_2, range_catch, exp_r, 
                       exp_f, exp_b, interval, multiplier), collapse = "_"))



pol_pars %>% print(n = Inf, width = Inf)

### plot MSE time series ####
res <- lapply(pol_pars$seq, function(x) {
  res <- readRDS(paste0("output/500_50/ms/trial/", pol_pars$fhist[x], "/pol/",
                        "mp_", pol_pars$file[x], 
                        ".rds"))[[1]]
  ### stock metrics
  SSBs <- FLCore::ssb(res@stock)
  Fs <- FLCore::fbar(res@stock)
  Cs <- FLCore::catch(res@stock)
  yrs <- dim(SSBs)[2]
  its <- dim(SSBs)[6]
  ### collapse correction
  if (isTRUE(TRUE)) {
    ### find collapses
    cd <- sapply(seq(its), function(x) {
      min_yr <- min(which(SSBs[,,,,, x] < 1))
      if (is.finite(min_yr)) {
        all_yrs <- min_yr:yrs
      } else {
        all_yrs <- NA
      }
      all_yrs + (x - 1)*yrs
    })
    cd <- unlist(cd)
    cd <- cd[which(!is.na(cd))]
    ### remove values
    SSBs@.Data[cd] <- 0
    Cs@.Data[cd] <- 0
    Fs@.Data[cd] <- 0
  }
  input <- readRDS(paste0("input/500_50/OM_2_mp_input/", pol_pars$fhist[x], 
                          "/pol.rds"))
  Bmsy <- c(input$refpts["msy", "ssb"])
  Fmsy <- c(input$refpts["msy", "harvest"])
  Cmsy <- c(input$refpts["msy", "yield"])
  Blim <- input$Blim
  ### summarise
  qnts <- FLQuants(SSB = SSBs/Bmsy, F = Fs/Fmsy, Catch = Cs/Cmsy)
  qnts <- lapply(qnts, quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
  df_qnts <- as.data.frame(qnts)[, c("year", "iter", "data", "qname")]
  df_qnts <- bind_cols(df_qnts, pol_pars[rep(x, nrow(df_qnts)), ])
  df_qnts <- df_qnts %>% tidyr::spread(key = iter, value = data)
  
  return(df_qnts)
  
})
res_df <- do.call(rbind, res)
res_df <- res_df %>%
  mutate(year = year - 100,
         group = ifelse(year < 1, "history", "projection"))
res_df <- res_df %>% bind_rows(res_df %>% filter(year == 0) %>%
                                 mutate(group = "projection"))
res_df <- res_df %>% 
  mutate(group = factor(group, levels = c("history", "projection")))
res_df <- res_df %>% 
  filter(year >= -25)
res_df$obj_def[res_df$obj_def == ""] <- "not optimised"
res_df$obj_def[res_df$obj_def == "C"] <- "Catch"
res_df$obj_def[res_df$obj_def == "SSB"] <- "SSB"
res_df$obj_def[res_df$obj_def == "SSB_C_risk_ICV"] <- "SSB+Catch+\nrisk+ICV"
res_df$obj_def[res_df$obj_def == "SSB_F_C_risk_ICV"] <- "SSB+F+Catch+\nrisk+ICV"
res_df$obj_def[res_df$obj_def == "SSB_risk_ICV"] <- "SSB+risk+ICV"
res_df$obj_def <- as.factor(res_df$obj_def)
res_df$obj_def <- factor(res_df$obj_def, levels = levels(res_df$obj_def)[c(2, 1, 3, 6, 4, 5)])
res_df$`50%` <- ifelse(res_df$group == "history" & res_df$optimised == TRUE, 
                       NA, res_df$`50%`)


plot_ssb_rel <- res_df %>% filter(qname == "SSB") %>%
  ggplot(aes(x = year, y = `50%`, colour = obj_def, linetype = obj_def)) +
  geom_ribbon(data = res_df %>% filter(qname == "SSB" & optimised == FALSE),
              aes(x = year, ymin = `5%`, ymax = `95%`, fill = "90%"), 
              linetype = 0, colour = 0, show.legend = FALSE, alpha = 0.5) +
  geom_ribbon(data = res_df %>% filter(qname == "SSB" & optimised == FALSE),
              aes(x = year, ymin = `25%`, ymax = `75%`, fill = "50%"),
              linetype = 0, colour = 0, show.legend = FALSE, alpha = 0.5) +
  scale_fill_manual("confidence\ninterval",
                    values = c("50%" = "grey50", "90%" = "grey80")) +
  geom_line(size = 0.3, show.legend = FALSE) +
  scale_colour_manual("fitness function", 
    values = setNames(c("black", "black", scales::hue_pal()(4)),
                      c("not optimised", "SSB+Catch+\nrisk+ICV", "Catch", "SSB",
                        "SSB+risk+ICV", "SSB+F+Catch+\nrisk+ICV"))) +
  scale_linetype_manual("fitness function", 
    values = setNames(c("solid", 11, 22, 1121, 31, 33),#c(1, 2, 3:6),
                      c("not optimised", "SSB+Catch+\nrisk+ICV", "Catch", "SSB",
                        "SSB+risk+ICV", "SSB+F+Catch+\nrisk+ICV"))) +
  facet_grid(fhist ~ group, scales = "free_x", space = "free_x") +
  theme_bw(base_size = 8, base_family = "sans") +
  scale_x_continuous(breaks = seq(-20, 50, 10), expand = c(0, 0)) +
  ylim(c(0, 4.1)) +
  labs(y = expression(SSB/B[MSY]), x = "") +
  theme(panel.spacing.x = unit(0, units = "cm"),
        panel.spacing.y = unit(1, units = "pt"),
        axis.title.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = unit(x = c(1, 3, 1, 3), units = "pt"))

plot_fbar_rel <- res_df %>% filter(qname == "F") %>%
  ggplot(aes(x = year, y = `50%`, colour = obj_def, linetype = obj_def)) +
  geom_ribbon(data = res_df %>% filter(qname == "F" & optimised == FALSE),
              aes(x = year, ymin = `5%`, ymax = `95%`, fill = "90%"), 
              linetype = 0, colour = 0, show.legend = FALSE, alpha = 0.5) +
  geom_ribbon(data = res_df %>% filter(qname == "F" & optimised == FALSE),
              aes(x = year, ymin = `25%`, ymax = `75%`, fill = "50%"),
              linetype = 0, colour = 0, show.legend = FALSE, alpha = 0.5) +
  scale_fill_manual("confidence\ninterval",
                    values = c("50%" = "grey50", "90%" = "grey80")) +
  geom_line(size = 0.3, show.legend = FALSE) +
  scale_colour_manual("fitness function", 
    values = setNames(c("black", "black", scales::hue_pal()(4)),
                      c("not optimised", "SSB+Catch+\nrisk+ICV", "Catch", "SSB",
                        "SSB+risk+ICV", "SSB+F+Catch+\nrisk+ICV"))) +
  scale_linetype_manual("fitness function", 
    values = setNames(c("solid", 11, 22, 1121, 31, 33),#c(1, 2, 3:6),
                      c("not optimised", "SSB+Catch+\nrisk+ICV", "Catch", "SSB",
                        "SSB+risk+ICV", "SSB+F+Catch+\nrisk+ICV"))) +
  facet_grid(fhist ~ group, scales = "free_x", space = "free_x") +
  theme_bw(base_size = 8, base_family = "sans") +
  scale_x_continuous(breaks = seq(-20, 50, 10), expand = c(0, 0)) +
  labs(y = expression(F/F[MSY]), x = "year") +
  theme(panel.spacing.x = unit(0, units = "cm"),
        panel.spacing.y = unit(1, units = "pt"),
        axis.title.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text.x = element_blank(),
        plot.margin = unit(x = c(1, 3, 1, 3), units = "pt"))

plot_catch_rel <- res_df %>% filter(qname == "Catch") %>%
  ggplot(aes(x = year, y = `50%`, colour = obj_def, linetype = obj_def)) +
  geom_ribbon(data = res_df %>% filter(qname == "Catch" & optimised == FALSE),
              aes(x = year, ymin = `5%`, ymax = `95%`, fill = "90%"), 
              linetype = 0, colour = 0, alpha = 0.5) +
  geom_ribbon(data = res_df %>% filter(qname == "Catch" & optimised == FALSE),
              aes(x = year, ymin = `25%`, ymax = `75%`, fill = "50%"),
              linetype = 0, colour = 0, alpha = 0.5) +
  scale_fill_manual("confidence\nintervals",
                    values = c("50%" = "grey50", "90%" = "grey80")) +
  geom_line(size = 0.3) +
  scale_colour_manual("fitness\nfunction", 
    values = setNames(c("black", "black", scales::hue_pal()(4)),
                      c("not optimised", "SSB+Catch+\nrisk+ICV", "Catch", "SSB",
                        "SSB+risk+ICV", "SSB+F+Catch+\nrisk+ICV"))) +
  scale_linetype_manual("fitness\nfunction", 
    values = setNames(c("solid", 11, 22, 1121, 31, 33),#c(1, 2, 3:6),
                      c("not optimised", "SSB+Catch+\nrisk+ICV", "Catch", "SSB",
                        "SSB+risk+ICV", "SSB+F+Catch+\nrisk+ICV"))) +
  facet_grid(fhist ~ group, scales = "free_x", space = "free_x") +
  theme_bw(base_size = 8, base_family = "sans") +
  scale_x_continuous(breaks = seq(-20, 50, 10), expand = c(0, 0)) +
  ylim(0, 3.1) +
  labs(y = expression(Catch/MSY), x = "year") +
  theme(panel.spacing.x = unit(0, units = "cm"),
        panel.spacing.y = unit(1, units = "pt"), 
        legend.position = "right",
        legend.key.height = unit(1, "lines"),
        legend.key.width = unit(1, "lines"),
        strip.text.x = element_blank(),
        plot.margin = unit(x = c(1, 3, 1, 3), units = "pt"))

p_pol_trajectories <- plot_grid(
  plot_grid(plot_ssb_rel, plot_fbar_rel, 
            plot_catch_rel + theme(legend.position = "none"),
            ncol = 1, align = "v", rel_heights = c(1.1, 1, 1.2)),
  get_legend(plot_catch_rel), rel_widths = c(1, 0.4), ncol = 2)

ggsave(filename = "output/plots/pol_trials.png", plot = p_pol_trajectories,
       width = 8.5, height = 10, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/pol_trials.pdf", plot = p_pol_trajectories,
       width = 8.5, height = 8.5, units = "cm", dpi = 600)


### add GA progress to same plot
pol_ga <- readRDS("output/500_50/ms/trial/one-way/pol/lag_idx-range_idx_1-range_idx_2-exp_r-exp_f-exp_b-interval-multiplier--obj_SSB_C_risk_ICV_res.rds")
ga_df <- as.data.frame(pol_ga@summary)
ga_df$generation <- seq(nrow(ga_df))
# ga_df <- rbind(ga_df, c(rep(NA, 6), 0))
# ga_df <- rbind(ga_df, c(rep(NA, 6), 30))
p_ga <- ga_df %>% 
  mutate(best = max) %>%
  pivot_longer(c(best, median), names_to = "key", values_to = "value") %>%
  ggplot(aes(x = generation, y = value, linetype = key, shape = key,
             ymin = min, ymax = max)) +
  geom_ribbon(alpha = 0.5, show.legend = FALSE, fill = "grey80") +
  geom_line(size = 0.3) +
  geom_point(size = 0.7) + 
  facet_grid("one-way" ~ NA) +
  scale_linetype_discrete("") + scale_shape_discrete("") +
  theme_bw(base_size = 8, base_family = "sans") +
  ylim(c(NA, 0)) +
  scale_x_continuous(expand = expansion(add = c(1, 1))) +
  scale_y_continuous(expand = expansion(add = c(0.4, abs(max(ga_df$max))))) +
  labs(x = "generations", y = "fitness value") +
  theme(legend.position = "right", legend.background = element_blank(),
        legend.key = element_blank(), 
        legend.key.size = unit(0.5, "lines"), legend.key.width = unit(1, "lines"),
        strip.text.x = element_blank(), legend.title = element_blank())

### combine
p_pol <- plot_grid(p_pol_trajectories, 
          plot_grid(plot_grid(NULL, p_ga + theme(legend.position = "none"),
                              rel_widths = c(0.05, 1), nrow = 1), 
                    get_legend(p_ga), nrow = 1, rel_widths = c(1, 0.4)), 
          ncol = 1, rel_heights = c(4, 1), 
          labels = c("(a)", "(b)"), label_size = 10, hjust = -0.1, vjust = 1.1)
ggsave(filename = "output/plots/pol_trials_combined.png", plot = p_pol,
       width = 8.5, height = 13, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/pol_trials.pdf", plot = p_pol,
       width = 8.5, height = 13, units = "cm", dpi = 600)


### plot stats ####

pol_pars %>% print(n = Inf, width = Inf)

stats_pol <- lapply(pol_pars$seq, function(x) {#browser()
  if (isTRUE(as.logical(pol_pars$optimised[x]))) {
    ### stats from optimised solution
    tmp <- readRDS(paste0("output/500_50/ms/trial/", pol_pars$fhist[x], "/pol/",
                          "lag_idx-range_idx_1-range_idx_2-exp_r-exp_f-exp_b",
                          "-interval-multiplier--obj_", pol_pars$obj_def[x],
                          "_runs.rds"))
    ### add optimised solution
    tmp <- as.list(as.data.frame(t(tmp[[pol_pars$file[x]]]$stats)))
    tmp <- lapply(tmp, "[[", 1)
    tmp <- c(as.list(pol_pars[x, ]), tmp)
    ### summary of genetic algorithm
    ga <- readRDS(paste0("output/500_50/ms/trial/", pol_pars$fhist[x], "/pol/",
                          "lag_idx-range_idx_1-range_idx_2-exp_r-exp_f-exp_b",
                          "-interval-multiplier--obj_", pol_pars$obj_def[x],
                          "_res.rds"))
    tmp$ga_popSize <- ga@popSize
    tmp$ga_iter <- ga@iter
    tmp$ga_fitnessValue <- ga@fitnessValue
  } else {
    ### extract stats from default objective function run
    tmp <- readRDS(paste0("output/500_50/ms/trial/", pol_pars$fhist[x], "/pol/",
                          "lag_idx-range_idx_1-range_idx_2-exp_r-exp_f-exp_b",
                          "-interval-multiplier--obj_", "SSB_C_risk_ICV",
                          "_runs.rds"))
    tmp <- as.list(as.data.frame(t(tmp[[pol_pars$file[x]]]$stats)))
    tmp <- lapply(tmp, "[[", 1)
    tmp <- c(as.list(pol_pars[x, ]), tmp)
    tmp$ga_popSize <- NA
    tmp$ga_iter <- NA
    tmp$ga_fitnessValue <- NA
  }
  return(tmp)
})
stats_pol <- as.data.frame(do.call(rbind, stats_pol))
stats_pol <- as.data.frame(lapply(stats_pol, unlist))
### format labels for plotting
stats_pol$obj_label <- as.character(stats_pol$obj_def)
stats_pol$obj_label[stats_pol$obj_label == ""] <- "not optimised"
stats_pol$obj_label[stats_pol$obj_label == "C"] <- "Catch"
stats_pol$obj_label[stats_pol$obj_label == "SSB"] <- "SSB"
stats_pol$obj_label[stats_pol$obj_label == "SSB_C_risk_ICV"] <- "SSB+Catch+\nrisk+ICV"
stats_pol$obj_label[stats_pol$obj_label == "SSB_F_C_risk_ICV"] <- "SSB+F+Catch+\nrisk+ICV"
stats_pol$obj_label[stats_pol$obj_label == "SSB_risk_ICV"] <- "SSB+risk+ICV"
stats_pol$obj_label <- as.factor(stats_pol$obj_label)
stats_pol$obj_label <- factor(stats_pol$obj_label, 
                              levels = levels(stats_pol$obj_label)[c(2, 1, 3, 6, 4, 5)])

stats_pol <- stats_pol %>%
  pivot_longer(c(SSB_rel, Fbar_rel, Catch_rel, risk_Blim, ICV, ga_fitnessValue), 
               names_to = "key", values_to = "value") %>%
  mutate(stat = factor(key, levels = c("SSB_rel", "Fbar_rel", "Catch_rel", "risk_Blim", 
                                       "ICV", "ga_fitnessValue"), 
                       labels = c("SSB/B[MSY]", "F/F[MSY]", "Catch/MSY", 
                                  "B[lim]~risk", "ICV", "fitness~value")))

p_pol_stats <- stats_pol %>% 
  ggplot(aes(x = obj_label, y = value, fill = obj_label,
             colour = obj_label)) +
  geom_col(position = "dodge", show.legend = FALSE, width = 0.8) +
  facet_grid(stat ~ fhist, scales = "free", space = "free_x", switch = "y",
             labeller = "label_parsed") +
  labs(y = "", x = "fitness function") +
  theme_bw(base_size = 8, base_family = "sans") +
  theme(panel.spacing.x = unit(0, units = "cm"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.placement.y = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        plot.margin = unit(x = c(1, 3, 3, 3), units = "pt"))

ggsave(filename = "output/plots/pol_trials_stats.png",
       width = 8.5, height = 13, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/pol_trials_stats.pdf",
       width = 8.5, height = 13, units = "cm", dpi = 600)


### combine all plots ####
p <- plot_grid(p_pol, p_pol_stats, ncol = 2, rel_widths = c(1.2, 1),
               labels = c("", "(c)"), label_size = 10, hjust = -0.1, vjust = 1.1)
ggsave(filename = "output/plots/pol_combined.png", plot = p,
       width = 17, height = 17, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/pol_combined.pdf", plot = p,
       width = 17, height = 13, units = "cm", dpi = 600)

### ------------------------------------------------------------------------ ###
### all stocks stats: default vs. optimised ####
### ------------------------------------------------------------------------ ###

### optimised parameters
res_par <- readRDS(paste0("output/", n_iter, "_", n_yrs, 
                          "/ms/trial/all_stocks_", fhist, "_opt_pars.rds"))

res_par <- t(as.data.frame(res_par))
res_par <- as.data.frame(res_par)
res_par$stock <- row.names(res_par)
res_par$fhist <- "one-way"
res_par$optimised <- TRUE
tmp <- expand.grid(1, 2, 3, 1, 1, 1, 1, 2, 1, stocks$stock, "one-way", FALSE)
names(tmp) <- names(res_par)
res_par <- rbind(res_par, tmp)
res_par$seq <- seq(nrow(res_par))
res_par <- res_par %>% 
  group_by(seq) %>%
  mutate(file = paste0(c(lag_idx, range_idx_1, range_idx_2, range_catch, exp_r, 
                         exp_f, exp_b, interval, multiplier), collapse = "_"))
res_par %>% print(n = Inf, width = Inf)

stats <- lapply(res_par$seq, function(x) {#browser()
  if (isTRUE(as.logical(res_par$optimised[x]))) {
    ### stats from optimised solution
    tmp <- readRDS(paste0("output/500_50/ms/trial/", res_par$fhist[x], "/", 
                          res_par$stock[x], "/",
                          "lag_idx-range_idx_1-range_idx_2-exp_r-exp_f-exp_b",
                          "-interval-multiplier--obj_", "SSB_C_risk_ICV",
                          "_runs.rds"))
    ### add optimised solution
    tmp <- as.list(as.data.frame(t(tmp[[res_par$file[x]]]$stats)))
    tmp <- lapply(tmp, "[[", 1)
    tmp <- c(as.list(res_par[x, ]), tmp)
    ### summary of genetic algorithm
    ga <- readRDS(paste0("output/500_50/ms/trial/", res_par$fhist[x], "/", 
                         res_par$stock[x], "/",
                         "lag_idx-range_idx_1-range_idx_2-exp_r-exp_f-exp_b",
                         "-interval-multiplier--obj_", "SSB_C_risk_ICV",
                         "_res.rds"))
    tmp$ga_popSize <- ga@popSize
    tmp$ga_iter <- ga@iter
    tmp$ga_fitnessValue <- ga@fitnessValue
  } else {
    ### extract stats from default objective function run
    tmp <- readRDS(paste0("output/500_50/ms/trial/", res_par$fhist[x], "/", 
                          res_par$stock[x], "/",
                          "lag_idx-range_idx_1-range_idx_2-exp_r-exp_f-exp_b",
                          "-interval-multiplier--obj_", "SSB_C_risk_ICV",
                          "_runs.rds"))
    tmp <- as.list(as.data.frame(t(tmp[[res_par$file[x]]]$stats)))
    tmp <- lapply(tmp, "[[", 1)
    tmp <- c(as.list(res_par[x, ]), tmp)
    tmp$ga_popSize <- NA
    tmp$ga_iter <- NA
    tmp$ga_fitnessValue <- -sum(c(abs(tmp$SSB_rel - 1), abs(tmp$Catch_rel - 1),
                                  tmp$risk_Blim, tmp$ICV))#NA
  }
  return(tmp)
})
stats <- as.data.frame(do.call(rbind, stats))
stats <- as.data.frame(lapply(stats, unlist))
stats <- full_join(stats, stocks[, c("stock", "k")])
stats$stock <- factor(stats$stock, levels = stocks$stock[order(stocks$k)])
saveRDS(stats, file = "output/500_50/ms/trial/all_stocks_one-way_stats.rds")
stats <- readRDS("output/500_50/ms/trial/all_stocks_one-way_stats.rds")

stats_plot <- stats %>% 
  pivot_longer(c(SSB_rel, Fbar_rel, Catch_rel, risk_Blim, ICV, 
                 ga_fitnessValue)) %>%
  mutate(name = ifelse(name == "SSB_rel", "SSB/B[MSY]", name),
         name = ifelse(name == "Fbar_rel", "F/F[MSY]", name),
         name = ifelse(name == "Catch_rel", "Catch/MSY", name),
         name = ifelse(name == "risk_Blim", "B[lim]~risk", name),
         name = ifelse(name == "ICV", "ICV", name),
         name = ifelse(name == "ga_fitnessValue", "fitness~value", 
                       name)) %>%
  mutate(name = factor(name, levels = unique(name)[c(1, 2, 3, 4, 5, 6)]),
         rule = ifelse(optimised == TRUE, "optimised", "default"))

p_stats <- stats_plot %>%
  ggplot(aes(x = stock, y = value, fill = rule, colour = rule)) +
  geom_col(position = position_dodge2(padding = 0.4), width = 0.8) +
  scale_colour_discrete("catch\nrule") +
  scale_fill_discrete("catch\nrule") + 
  facet_grid(name ~ "stock~specific~optimisation", 
             scales = "free_y", switch = "y", 
             labeller = "label_parsed") +
  labs(x = "stock", y = "") +
  theme_bw(base_size = 8, base_family = "sans") +
  theme(strip.placement.y = "outside",
        strip.background.y = element_blank(), 
        strip.text.y = element_text(size = 8),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.title.y = element_blank(), 
        legend.key.size = unit(0.3, "lines"))

ggsave(filename = "output/plots/all_stocks_stats.png", plot = p_stats,
       width = 8.5, height = 12, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/all_stocks_stats.pdf", plot = p_stats,
       width = 8.5, height = 12, units = "cm", dpi = 600)

### ------------------------------------------------------------------------ ###
### multi-species stats ####
### ------------------------------------------------------------------------ ###

### optimised parameters
res_par <- readRDS("output/500_50/ms/trial/ms_res.rds")

res_par <- t(as.data.frame(res_par))
res_par <- as.data.frame(res_par)
res_par$stocks <- row.names(res_par)
res_par$fhist <- "one-way"
res_par$optimised <- TRUE
res_par$group <- c("low", "medium", "high")
tmp <- data.frame(1, 2, 3, 1, 1, 1, 1, 2, 1, res_par$stocks, "one-way", FALSE,
                  res_par$group)
names(tmp) <- names(res_par)
res_par <- rbind(res_par, tmp)
res_par$seq <- seq(nrow(res_par))
res_par <- res_par %>% 
  group_by(seq) %>%
  mutate(file = paste0(c(lag_idx, range_idx_1, range_idx_2, range_catch, exp_r, 
                         exp_f, exp_b, interval, multiplier), collapse = "_"))
res_par %>% print(n = Inf, width = Inf)

stats_ms <- lapply(res_par$seq, function(x) {#browser()
  ### stats by run
  runs <- readRDS(paste0("output/500_50/ms/trial/one-way/", 
                         res_par$stocks[x], "/",
                         "lag_idx-range_idx_1-range_idx_2-exp_r-exp_f-exp_b",
                         "-interval-multiplier--obj_", "SSB_C_risk_ICV",
                         "_runs.rds"))
  tmp <- as.data.frame(t(runs[[res_par$file[x]]]$stats))
  tmp$stock <- row.names(tmp)
  if (isTRUE(as.logical(res_par$optimised[x]))) {
    ### summary of genetic algorithm
    ga <- readRDS(paste0("output/500_50/ms/trial/one-way/", 
                         res_par$stocks[x], "/",
                         "lag_idx-range_idx_1-range_idx_2-exp_r-exp_f-exp_b",
                         "-interval-multiplier--obj_", "SSB_C_risk_ICV",
                         "_res.rds"))
    tmp$ga_popSize <- ga@popSize
    tmp$ga_iter <- ga@iter
    tmp$ga_fitnessValue <- ga@fitnessValue
  } else {
    tmp$ga_popSize <- NA
    tmp$ga_iter <- NA
    tmp$ga_fitnessValue <- -sum(c(abs(unlist(tmp$SSB_rel) - 1), 
                                  abs(unlist(tmp$Catch_rel) - 1),
                                  unlist(tmp$risk_Blim), unlist(tmp$ICV)))#NA
  }
  tmp <- bind_cols(res_par[rep(x, nrow(tmp)), ], tmp)
  return(tmp)
})
stats_ms <- as.data.frame(do.call(rbind, stats_ms))
stats_ms <- as.data.frame(lapply(stats_ms, unlist))
stats_ms <- full_join(stats_ms, stocks[, c("stock", "k")])
stats_ms$stock <- factor(stats_ms$stock, levels = stocks$stock[order(stocks$k)])
saveRDS(stats_ms, file = "output/500_50/ms/trial/ms_stats.rds")
stats_ms <- readRDS("output/500_50/ms/trial/ms_stats.rds")

stats_ms_plot <- stats_ms %>% 
  select(stocks, stock, group, optimised, ga_fitnessValue) %>%
  group_by(group, optimised) %>%
  summarise(ga_fitnessValue = first(ga_fitnessValue)) %>%
  mutate(stock = "group", k = Inf) %>%
  ungroup() %>%
  full_join(stats_ms)
stats_ms_plot$stock <- factor(stats_ms_plot$stock, 
                              levels = c(stocks$stock[order(stocks$k)], "group"))
stats_ms_plot <- stats_ms_plot %>%
  mutate(ga_fitnessValue = ifelse(stock == "group", ga_fitnessValue, NA))

stats_ms_plot <- stats_ms_plot %>% 
  pivot_longer(c(SSB_rel, Fbar_rel, Catch_rel, risk_Blim, ICV, 
                 ga_fitnessValue)) %>%
  mutate(name = ifelse(name == "SSB_rel", "SSB/B[MSY]", name),
         name = ifelse(name == "Fbar_rel", "F/F[MSY]", name),
         name = ifelse(name == "Catch_rel", "Catch/MSY", name),
         name = ifelse(name == "risk_Blim", "B[lim]~risk", name),
         name = ifelse(name == "ICV", "ICV", name),
         name = ifelse(name == "ga_fitnessValue", "fitness~value", 
                       name)) %>%
  mutate(name = factor(name, levels = unique(name)[c(1, 2, 3, 4, 5, 6)]),
         rule = ifelse(optimised == TRUE, "optimised", "default"),
         k_group = paste0(group, "-italic(k)~group"),
         k_group = factor(k_group, levels = unique(k_group)[c(2, 3, 1)]))

p_stats_ms <- stats_ms_plot %>%
  filter(name %in% c("SSB/B[MSY]", "Catch/MSY", "ICV", "F/F[MSY]", 
                     "B[lim]~risk", "fitness~value")) %>%
  ggplot(aes(x = stock, y = value, fill = rule, colour = rule)) +
  geom_col(position = position_dodge2(padding = 0.4), width = 0.8) +
  scale_colour_discrete("catch rule") +
  scale_fill_discrete("catch rule") + 
  facet_grid(name ~ k_group, space = "free_x", switch = "y",
             scales = "free",# strip.position = "left", 
             labeller = "label_parsed") +
  labs(x = "stock", y = "") +
  theme_bw(base_size = 8, base_family = "sans") +
  theme(strip.placement.y = "outside",
        strip.background.y = element_blank(), 
        panel.spacing.x = unit(0, units = "cm"),
        strip.text.y = element_text(size = 8),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.title.y = element_blank(), 
        legend.key.size = unit(0.3, "lines"),
        legend.position = c(0.15, 0.1),
        legend.background = element_blank())

ggsave(filename = "output/plots/ms_stats.png", plot = p_stats_ms,
       width = 8.5, height = 12, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/ms_stats.pdf", plot = p_stats_ms,
       width = 8.5, height = 12, units = "cm", dpi = 600)


### combine single and multi-species stats plots ####
plot_grid(plot_grid(NULL, p_stats + theme(legend.position = "none"),
                    ncol = 1, rel_heights = c(0, 1)),#c(0.0395, 1)), 
          plot_grid(p_stats_ms,# + theme(legend.position = "none"),
                    NULL,#get_legend(p_stats_ms),
                    ncol = 1, rel_heights = c(1, 0)),#0.17)),
          labels = c("(a)", "(b)"), label_size = 10, rel_widths = c(0.95, 1.05))
ggsave(filename = "output/plots/stats_ms_and_single.png",
       width = 17, height = 12, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/stats_ms_and_single.pdf",
       width = 17, height = 12, units = "cm", dpi = 600)

