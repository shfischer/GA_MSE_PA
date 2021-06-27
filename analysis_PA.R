### ------------------------------------------------------------------------ ###
### analyse MSY GA runs with PA fitness function ####
### ------------------------------------------------------------------------ ###

req_pckgs <- c("doParallel", "doRNG", "mse", "GA", "ggplot2", "cowplot", 
               "scales", "Cairo", "tidyr", "dplyr", "FLCore", "FLash", "FLBRP")
for (i in req_pckgs) library(package = i, character.only = TRUE)
library(RColorBrewer)

### function for transforming y-axis origin in plots
trans_from <- function(from = 1) {
  trans <- function(x) x - from
  inv <- function(x) x + from
  trans_new("from", trans, inv, 
            domain = c(from, Inf))
}

### load stock specifications
stocks <- read.csv("input/stocks.csv", stringsAsFactors = FALSE)

source("funs_GA.R")
source("funs.R")

### ------------------------------------------------------------------------ ###
### plot penalty function ####
### ------------------------------------------------------------------------ ###

penalty <- function(x, negative = FALSE, max = 5,
                    inflection = 0.06, steepness = 0.5e+3) {
  y <- max / (1 + exp(-(x - inflection)*steepness))
  if (isTRUE(negative)) y <- -y
  return(y)
}

p <- ggplot() +
  geom_function(fun = penalty, n = 1000) +
  geom_vline(xintercept = 0.05, colour = "red") +
  theme_bw(base_size = 8, base_family = "sans") +
  #scale_x_continuous(expand = c(0, 0)) +
  xlim(c(0, 1)) + 
  coord_cartesian(xlim = c(0, 0.2)) +
  labs(x = expression(B[lim]~risk),
       y = expression(penalty))
p
ggsave(filename = "output/plots/PA/Blim_penalty_curve.png",
       width = 8.5, height = 6, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/PA/Blim_penalty_curve.pdf",
      width = 8.5, height = 6, units = "cm", dpi = 600)


### ------------------------------------------------------------------------ ###
### collate results - pollack rfb-rule component explorations ####
### ------------------------------------------------------------------------ ###
n_yrs <- 50
n_iter <- 500

### get optimised parameterisation
pol_PA <- foreach(stock = "pol", .combine = rbind) %:%
  foreach(optimised = c("default", "mult", "cap", "mult_cap", "all", "all_cap"),
          .combine = rbind) %:%
  foreach(scenario = "PA", .combine = rbind) %:%
  foreach(stat_yrs = "more", .combine = rbind) %:%
  foreach(fhist = c("one-way", "random"), .combine = rbind) %do% {#browser()
    #browser()
    ### find files
    file_name <- switch(optimised,
      "default" = "multiplier--obj_ICES_MSYPA",
      "mult" = "multiplier--obj_ICES_MSYPA",
      "cap" = "upper_constraint-lower_constraint--obj_ICES_MSYPA",
      "mult_cap" = paste0("multiplier-upper_constraint-lower_constraint--",
                          "obj_ICES_MSYPA"),
      "all" = paste0("lag_idx-range_idx_1-range_idx_2-exp_r-exp_f-exp_b-",
                     "interval-multiplier--obj_ICES_MSYPA"),
      "all_cap" = paste0("lag_idx-range_idx_1-range_idx_2-exp_r-exp_f-exp_b-",
                         "interval-multiplier-upper_constraint-lower_constraint",
                         "--obj_ICES_MSYPA"))
    fs <- list.files(path = paste0("output/", n_iter, "_", n_yrs, "/", scenario, 
                                   "/", fhist, "/", stock, "/"), 
                     pattern = paste0("^", file_name, "_res_", stat_yrs, ".rds"))
    if (length(fs) < 1) return(NULL)
    trials <- data.frame(file = fs)
    trials$obj_fun <- "MSYPA"
    trials$fhist <- fhist
    trials$scenario <- scenario
    trials$optimised <- optimised
    trials$stat_yrs <- stat_yrs
    trials$stock <- stock
    ### load GA results
    res_lst <- lapply(trials$file, function(x) {
      readRDS(paste0("output/", n_iter, "_", n_yrs, "/", scenario, "/", fhist, 
                     "/", stock, "/", x))
    })
    res_par <- lapply(res_lst, function(x) {
      tmp <- x@solution[1,]
      tmp[c(1:4, 8)] <- round(tmp[c(1:4, 8)])
      tmp[5:7] <- round(tmp[5:7], 1)
      tmp[9] <- round(tmp[9], 2)
      tmp[10] <- ifelse(is.nan(tmp[10]), Inf, tmp[10])
      tmp[10:11] <- round(tmp[10:11], 2)
      return(tmp)
    })
    names(res_par) <- trials$obj_fun
    ### default (non-optimised?)
    if (identical(optimised, "default")) {
      trials <- trials[1, ]
      res_lst <- res_lst[1]
      res_par <- res_par[1]
      res_par[[1]][] <- c(1, 2, 3, 1, 1, 1, 1, 2, 1, Inf, 0)
    }
    res <- as.data.frame(do.call(rbind, res_par))
    res$file <- trials$file
    res$solution <- sapply(res_par, paste0, collapse = "_")
    res$fitness <- sapply(res_lst, slot, "fitnessValue")
    res$iter <- sapply(res_lst, slot, "iter")
    ### combine
    res <- merge(trials, res)
    ### load stats of solution
    res_stats <- lapply(split(res, seq(nrow(res))), function(x) {
      res_runs <- readRDS(paste0("output/", n_iter, "_", n_yrs, "/", scenario, 
                                 "/", fhist, "/", stock, "/",
                                 gsub(x = x$file, 
                                      pattern = paste0("res_", stat_yrs), 
                                      replacement = "runs")))
      stats_i <- t(res_runs[[x$solution]]$stats)
      as.data.frame(lapply(as.data.frame(stats_i), unlist))
    })
    res_stats <- do.call(rbind, res_stats)
    res <- cbind(res, res_stats)
    ### calculate fitness value for non-optimised rule
    if (identical(optimised, "default")) {
      res$fitness <- sapply(split(res, seq(nrow(res))), function(x) {
        -sum(abs(x$SSB_rel - 1), abs(x$Catch_rel - 1), 
             penalty(x = unlist(x$risk_Blim), negative = FALSE, max = 5, 
                 inflection = 0.06, steepness = 0.5e+3),
             x$ICV)
      })
    }
    return(res)
}
saveRDS(pol_PA, file = "output/pol_PA_components_stats.rds")
write.csv(pol_PA, file = "output/pol_PA_components_stats.csv", row.names = FALSE)

### ------------------------------------------------------------------------ ###
### collate results - all stocks with multiplier ####
### ------------------------------------------------------------------------ ###

### get results
all_mult <- foreach(stock = stocks$stock, .combine = rbind) %:%
  foreach(fhist = c("one-way", "random"), .combine = rbind) %do% {
    
    runs <- readRDS(paste0("output/500_50/PA/", fhist, "/", stock,
                           "/multiplier--obj_ICES_MSYPA_runs.rds"))
    runs <- lapply(runs, function(x) {
      data.frame(c(as.list(x$pars), x$stats[, 1], 
                   stock = stock, fhist = fhist), 
                 stringsAsFactors = FALSE)
    })
    runs <- do.call(rbind, runs)
    row.names(runs) <- NULL
    return(runs)
}
names(all_mult)[12:22] <- paste0(names(all_mult)[12:22], "_all")
names(all_mult)[12:89] <- sapply(names(all_mult)[12:89], function(x) {
  x <- gsub(x = x, pattern = "_all", replacement = "-all")
  x <- gsub(x = x, pattern = "_first10", replacement = "-first10")
  x <- gsub(x = x, pattern = "_41to50", replacement = "-41to50")
  x <- gsub(x = x, pattern = "_last10", replacement = "-last10")
  x <- gsub(x = x, pattern = "_firsthalf", replacement = "-firsthalf")
  x <- gsub(x = x, pattern = "_lastfhalf", replacement = "-lastfhalf")
  x <- gsub(x = x, pattern = "_11to50", replacement = "-11to50")
  x
})
all_mult <- all_mult %>%
  pivot_longer(cols = "risk_Blim-all":"ICV-11to50",
               names_to = c("stat", "period"),
               names_pattern = "(.*)-(.*)",
               values_to = "value") %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  ### correct ICV when multiplier = 0
  mutate(ICV = ifelse(multiplier == 0, 0, ICV))

saveRDS(all_mult, file = "output/all_stocks_PA_multiplier_stats.rds")
write.csv(all_mult, file = "output/all_stocks_PA_multiplier_stats.csv", 
          row.names = FALSE)

### ------------------------------------------------------------------------ ###
### collate results - all stocks - all scenarios ####
### ------------------------------------------------------------------------ ###
n_yrs <- 50
n_iter <- 500

### get optimised parameterisation
all_GA <- foreach(stock = stocks$stock, .combine = bind_rows) %:%
  foreach(optimised = c("zero-fishing", "default", "mult", "all", "all_cap"),
          .combine = bind_rows) %:%
  foreach(capped = c(TRUE, FALSE), .combine = bind_rows) %:%
  foreach(scenario = c("PA", "MSY", "PA_capped"), .combine = bind_rows) %:%
  foreach(stat_yrs = "more", .combine = bind_rows) %:%
  foreach(fhist = c("one-way", "random"), .combine = bind_rows) %do% {#browser()
    
    ### objective function name
    name_obj <- switch(scenario,
                       "PA" = "ICES_MSYPA",
                       "PA_capped" = "ICES_MSYPA",
                       "MSY" = "SSB_C_risk_ICV")
    ### parameters used in optimisation
    name_pars <- switch(optimised,
      "zero-fishing" = "multiplier",
      "default" = "multiplier",
      "mult" = "multiplier",
      "cap" = "upper_constraint-lower_constraint",
      "mult_cap" = paste0("multiplier-upper_constraint-lower_constraint--"),
      "all" = paste0("lag_idx-range_idx_1-range_idx_2-exp_r-exp_f-exp_b-",
                     "interval-multiplier"),
      "all_cap" = paste0("lag_idx-range_idx_1-range_idx_2-exp_r-exp_f-exp_b-",
                         "interval-multiplier-upper_constraint-lower_constraint"))
    ### uncertainty cap fixed?
    if (isTRUE(capped)) {
      upper_cap <- 1.2
      lower_cap <- 0.7
      ### add caps if not part of file name
      if (isFALSE(grepl(x = name_pars, 
                       pattern = "upper_constraint-lower_constraint"))) {
        name_pars <- paste0(name_pars, 
                            "-upper_constraint", upper_cap, 
                            "-lower_constraint", lower_cap)
      ### if caps part of filename, add cap value
      } else {
        name_pars <- gsub(x = name_pars, pattern = "upper_constraint",
                          replacement = paste0("upper_constraint", upper_cap))
        name_pars <- gsub(x = name_pars, pattern = "lower_constraint",
                          replacement = paste0("lower_constraint", lower_cap))
      }
    }

    ### stats period
    name_stats <- switch(scenario,
                         "PA" = "_more",
                         "PA_capped" = "_more",
                         "MSY" = switch(optimised,
                                        "mult" = "_more",
                                        "zero-fishing" = "_more",
                                        "default" = "_more", 
                                        ""))
    ### assemble file name
    file_path <- paste0("output/", n_iter, "_", n_yrs, "/", scenario, 
                        "/", fhist, "/", stock, "/")
    file_name <- paste0(name_pars, "--obj_", name_obj)
    file_name_res <- paste0(file_name, "_res", name_stats, ".rds")
    file_name_runs <- paste0(file_name, "_runs", ".rds")
    if (isFALSE(file.exists(paste0(file_path, file_name_res)))) return(NULL)
    res_df <- data.frame(file = file_name_res)
    res_df$obj_fun <- scenario
    res_df$fhist <- fhist
    res_df$scenario <- scenario
    res_df$optimised <- optimised
    res_df$capped <- capped
    res_df$stat_yrs <- stat_yrs
    res_df$stock <- stock
    ### load GA results
    res <- readRDS(paste0(file_path, file_name_res))
    ### get parameters
    res_par <- res@solution[1, ]
    res_par["lag_idx"] <- round(res_par["lag_idx"])
    res_par["range_idx_1"] <- round(res_par["range_idx_1"])
    res_par["range_idx_2"] <- round(res_par["range_idx_2"])
    res_par["range_catch"] <- round(res_par["range_catch"])
    res_par["exp_r"] <- round(res_par["exp_r"], 1)
    res_par["exp_f"] <- round(res_par["exp_f"], 1)
    res_par["exp_b"] <- round(res_par["exp_b"], 1)
    res_par["interval"] <- round(res_par["interval"])
    res_par["multiplier"] <- round(res_par["multiplier"], 2)
    if ("upper_constraint" %in% names(res_par)) {
      res_par["lower_constraint"] <- round(res_par["lower_constraint"], 2)
      res_par["upper_constraint"] <- ifelse(!is.nan(res_par["upper_constraint"]),
                                            round(res_par["upper_constraint"], 2),
                                            Inf)
    }
    ### default or zero-fishing?
    if (identical(optimised, "zero-fishing")) res_par["multiplier"] <- 0
    if (identical(optimised, "default")) res_par["multiplier"] <- 1
    ### combine
    res_df <- cbind(as.data.frame(do.call(rbind, list(res_par))),
                    as.data.frame(do.call(rbind, list(res_df))))
    res_df$solution <- paste0(res_par, collapse = "_")
    res_df$fitness <- res@fitnessValue
    res_df$iter <- res@iter
    
    ### load stats of solution
    res_runs <- readRDS(paste0(file_path, file_name_runs))
    res_stats <- as.data.frame(
      lapply(as.data.frame(t(res_runs[[res_df$solution]]$stats)), unlist))
    res_df <- cbind(res_df, res_stats)
    
    ### calculate fitness value for non-optimised rule
    if (isTRUE(optimised %in% c("default", "zero-fishing"))) {
      if (isTRUE(scenario %in% c("PA", "PA_capped"))) {
        res_df$fitness <- -sum(abs(res_df$SSB_rel - 1),
                               abs(res_df$Catch_rel - 1),
                               res_df$ICV,
                               penalty(x = res_df$risk_Blim, negative = FALSE, 
                                       max = 5, inflection = 0.06, 
                                       steepness = 0.5e+3))
      } else if (isTRUE(scenario == "MSY")) {
        res_df$fitness <- -sum(abs(res_df$SSB_rel - 1),
                               abs(res_df$Catch_rel - 1),
                               res_df$ICV,
                               res_df$risk_Blim)
      }
    }
    return(res_df)
}
saveRDS(all_GA, file = "output/all_stocks_GA_optimised_stats.rds")
write.csv(all_GA, file = "output/all_stocks_GA_optimised_stats.csv",
          row.names = FALSE)

### ------------------------------------------------------------------------ ###
### plot - all stocks multipliers ####
### ------------------------------------------------------------------------ ###

all_mult <- readRDS("output/all_stocks_PA_multiplier_stats.rds")

### format for plotting
all_mult <- all_mult %>%
  filter(period == "all") %>%
  select(multiplier, risk_Blim, SSB_rel, Fbar_rel, Catch_rel, ICV, stock, 
         fhist) %>%
  mutate(ICV = ifelse(multiplier == 0, NA, ICV)) %>%
  pivot_longer(c(SSB_rel, Fbar_rel, Catch_rel, risk_Blim, ICV), 
               names_to = "key", values_to = "value") %>%
  mutate(stat = factor(key, 
                       levels = c("SSB_rel", "Fbar_rel", "Catch_rel", 
                                  "risk_Blim", "ICV"), 
                       labels = c("SSB/B[MSY]", "F/F[MSY]", "Catch/MSY", 
                                  "B[lim]~risk", "ICV")))
all_mult <- all_mult %>% 
  full_join(stocks[, c("stock", "k")]) %>%
  mutate(stock_label = paste0(stock, "~(k==", k, ")")) %>%
  mutate(stock_label = factor(stock_label, levels = unique(stock_label)))
### risk reference values
all_mult_refs <- bind_rows(
  all_mult %>%
    group_by(fhist, stock) %>%
    filter(key == "risk_Blim" & value < 0.055) %>%
    filter(multiplier == max(multiplier)) %>%
    mutate(point = "5%", key = NULL, value = NULL, stat = NULL) %>%
    left_join(all_mult),
  all_mult %>%
    group_by(fhist, stock) %>%
    filter(key == "risk_Blim") %>%
    filter(multiplier == 0) %>%
    mutate(F0_risk = value,
           value = NULL, multiplier = NULL) %>%
    full_join(all_mult) %>%
    filter(key == "risk_Blim") %>%
    filter(value <= (F0_risk + 0.05)) %>%
    filter(multiplier == max(multiplier)) %>%
    mutate(point = "F0+5%", F0_risk = NULL, key = NULL, 
           value = NULL, stat = NULL) %>%
    left_join(all_mult),
  all_mult %>%
    group_by(fhist, stock) %>%
    filter(key == "risk_Blim" & value < 0.105) %>%
    filter(multiplier == max(multiplier)) %>%
    mutate(point = "10%", key = NULL, value = NULL, stat = NULL) %>%
    left_join(all_mult)
) %>%
  mutate(point = factor(point, levels = unique(point)))

saveRDS(all_mult, file = "output/plots/PA/data_all_stocks_multiplier_stats.rds")
all_mult <- readRDS("output/plots/PA/data_all_stocks_multiplier_stats.rds")
saveRDS(all_mult_refs, 
        file = "output/plots/PA/data_all_stocks_multiplier_stats_refs.rds")
all_mult_refs <- readRDS("output/plots/PA/data_all_stocks_multiplier_stats_refs.rds")

### plot some example stocks
p_all_mult <- all_mult %>%
  filter(stock %in% c("ang3", "pol", "tur", "san") &
           stat != "F/F[MSY]") %>%
  ggplot(aes(x = multiplier, y = value, linetype = fhist,
             colour = fhist)) +
  geom_line(size = 0.4) +
  geom_point(data = all_mult_refs %>%
               filter(stock %in% c("ang3", "pol", "tur", "san") &
                        stat != "F/F[MSY]"),
             aes(shape = point), size = 1, stroke = 1) +
  scale_shape_manual("risk limit", values = 2:4, 
                     guide = guide_legend(order = 2)) + 
  scale_colour_brewer("fishing\nhistory", guide = guide_legend(order = 1),
                      palette = "Set1") +
  #scale_colour_discrete("fishing\nhistory", guide = guide_legend(order = 1)) +
  scale_linetype("fishing\nhistory", guide = guide_legend(order = 1)) +
  facet_grid(stat ~ stock_label, 
             scales = "free", space = "free_x", switch = "y",
             labeller = "label_parsed") +
  labs(y = "", x = "multiplier") +
  theme_bw(base_size = 8, base_family = "sans") +
  xlim(0, 1.3) +
  ylim(0, NA) +
  geom_blank(data = data.frame(multiplier = 1, value = 1, fhist = NA)) + ### for axis limits
  theme(strip.placement.y = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        legend.key.height = unit(1, "line"))
ggsave(filename = "output/plots/PA/all_stocks_mult_stats.png", 
       plot = p_all_mult,
       width = 17, height = 13, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/PA/all_stocks_mult_stats.pdf", 
       plot = p_all_mult,
       width = 17, height = 13, units = "cm", dpi = 600)

### plot all stocks for supplementary document
. <- foreach(subset = split(stocks$stock, rep(seq(6), 
                                              each = length(stocks$stock)/5)), 
             p = seq(6)) %do% {
  p_tmp <- all_mult %>%
    filter(stock %in% subset & stat != "F/F[MSY]") %>%
    ggplot(aes(x = multiplier, y = value, linetype = fhist,
               colour = fhist)) +
    geom_line(size = 0.4) +
    geom_point(data = all_mult_refs %>%
                 filter(stock %in% subset & stat != "F/F[MSY]"),
               aes(shape = point), size = 1, stroke = 1) +
    scale_shape_manual("risk limit", values = 2:4, 
                       guide = guide_legend(order = 2)) + 
    scale_colour_brewer("fishing\nhistory", guide = guide_legend(order = 1),
                        palette = "Set1") +
    scale_linetype("fishing\nhistory", guide = guide_legend(order = 1)) +
    facet_grid(stat ~ stock_label, 
               scales = "free", space = "free_x", switch = "y",
               labeller = "label_parsed") +
    labs(y = "", x = "multiplier") +
    theme_bw(base_size = 8, base_family = "sans") +
    xlim(0, 1.3) +
    ylim(0, NA) +
    geom_blank(data = data.frame(multiplier = 1, value = 1, fhist = NA)) + 
    theme(strip.placement.y = "outside",
          strip.background.y = element_blank(),
          strip.text.y = element_text(size = 8),
          legend.key.height = unit(1, "line"))
  ggsave(filename = paste0("output/plots/PA/all_stocks_mult_stats_", p ,".png"),
         plot = p_tmp,
         width = 17, height = 13, units = "cm", dpi = 600, type = "cairo")
  ggsave(filename = paste0("output/plots/PA/all_stocks_mult_stats_", p ,".pdf"),
         plot = p_tmp,
         width = 17, height = 13, units = "cm", dpi = 600)
}

### ------------------------------------------------------------------------ ###
### plot - all stocks multipliers - time periods ####
### ------------------------------------------------------------------------ ###

mult_periods <- readRDS("output/all_stocks_PA_multiplier_stats.rds")
brps <- readRDS("input/brps.rds")

### format for plotting
mult_periods <- mult_periods %>%
  filter(period %in% c("all", "first10", "last10")) %>%
  select(stock, fhist, period, multiplier, risk_Blim, Catch_rel) %>%
  pivot_longer(c(Catch_rel, risk_Blim), 
               names_to = "key", values_to = "value") %>%
  mutate(stat = factor(key, 
                       levels = c("risk_Blim", "Catch_rel"), 
                       labels = c("B[lim]~risk", "Catch/MSY")),
         period = factor(period,
                         levels = c("first10", "last10", "all"),
                         labels = c("first 10 years", 
                                    "last 10 years", 
                                    "all years")))
mult_periods <- mult_periods %>% 
  full_join(stocks[, c("stock", "k")]) %>%
  mutate(k = round(k, 2)) %>%
  mutate(stock_label = paste0(stock, "~(k==", k, ")")) %>%
  mutate(stock_label = factor(stock_label, levels = unique(stock_label)))

### risk reference values
mult_periods_5 <- mult_periods %>%
    group_by(fhist, stock, period) %>%
    filter(key == "risk_Blim" & value < 0.055) %>%
    filter(multiplier == max(multiplier)) %>%
    mutate(point = "5%", key = NULL, value = NULL, stat = NULL)

### plot example: pollack & herring
p_period <- mult_periods %>%
  filter(stock %in% c("pol", "her") & 
           fhist == "one-way") %>%
  ggplot(aes(x = multiplier, y = value,
             colour = period, linetype = period)) +
  geom_hline(data = data.frame(stat = "B[lim]~risk", y = 0.055),
             aes(yintercept = y), colour = "red", size = 0.4) +
  geom_line(size = 0.4) +
  scale_colour_manual("", 
                      values = brewer.pal(4, name = "Set1")[-1],
                      guide = guide_legend(order = 1)) +
  scale_linetype("", guide = guide_legend(order = 1)) +
  labs(y = "", 
       x = expression(multiplier~"("*italic(x)*")")) +
  facet_grid(stat ~ stock_label, scales = "free_y", labeller = label_parsed,
             switch = "y") +
  theme_bw(base_size = 8, base_family = "sans") +
  xlim(0, 2) + ylim(0, 1) +
  theme(legend.position = c(0.7, 0.4),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.key.height = unit(0.5, "lines"),
        strip.background.y = element_blank(),
        strip.placement.y = "outside",
        strip.text = element_text(size = 8),
        axis.title.y = element_blank())
ggsave(filename = "output/plots/PA/period_example.png",
       plot = p_period, width = 8.5, height = 7, units = "cm", dpi = 600, 
       type = "cairo")
ggsave(filename = "output/plots/PA/period_example.pdf",
       plot = p_period, width = 8.5, height = 7, units = "cm", dpi = 600)

### plot for all stocks
p_period_all <- foreach(stock = stocks$stock) %do% {
  p_tmp <- mult_periods %>%
    filter(stock == !!stock) %>%
    ggplot(aes(x = multiplier, y = value,
               colour = period, linetype = period)) +
    geom_hline(data = data.frame(stat = "B[lim]~risk", y = 0.055),
               aes(yintercept = y), colour = "red", size = 0.4) +
    geom_blank(data = data.frame(value = 0:1, multiplier = c(0, 2),
                                 period = NA, fhist = NA)) +
    geom_line(size = 0.4) +
    scale_colour_manual("", 
                        values = brewer.pal(4, name = "Set1")[-1],
                        guide = guide_legend(order = 1)) +
    scale_linetype("", guide = guide_legend(order = 1)) +
    labs(y = "", 
         x = expression(multiplier~"("*italic(x)*")")) +
    # facet_grid(stat ~ stock_label + fhist, scales = "free_y", 
    #            labeller = label_parsed, switch = "y") +
    facet_grid(stat ~ paste0(stock_label, "~'|'~", fhist), 
               scales = "free_y", labeller = label_parsed, switch = "y") +
    theme_bw(base_size = 8, base_family = "sans") +
    xlim(0, 2) + #ylim(0, 1) +
    theme(legend.position = c(0.7, 0.4),
          legend.background = element_blank(),
          legend.key = element_blank(),
          legend.key.height = unit(0.5, "lines"),
          strip.background.y = element_blank(),
          strip.placement.y = "outside",
          strip.text = element_text(size = 8),
          axis.title.y = element_blank())
  p_tmp
}
plot_grid(p_period_all[[1]], p_period_all[[2]], p_period_all[[3]],
          p_period_all[[4]], p_period_all[[5]], p_period_all[[6]],
          ncol = 2)
ggsave(filename = "output/plots/PA/period_example_1.png",
       width = 17, height = 21, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/PA/period_example_1.pdf",
       width = 17, height = 21, units = "cm", dpi = 600)
plot_grid(p_period_all[[7]], p_period_all[[8]], p_period_all[[9]],
          p_period_all[[10]], p_period_all[[11]], p_period_all[[12]],
          ncol = 2)
ggsave(filename = "output/plots/PA/period_example_2.png",
       width = 17, height = 21, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/PA/period_example_2.pdf",
       width = 17, height = 21, units = "cm", dpi = 600)
plot_grid(p_period_all[[13]], p_period_all[[14]], p_period_all[[15]],
          p_period_all[[16]], p_period_all[[17]], p_period_all[[18]],
          ncol = 2)
ggsave(filename = "output/plots/PA/period_example_3.png",
       width = 17, height = 21, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/PA/period_example_3.pdf",
       width = 17, height = 21, units = "cm", dpi = 600)
plot_grid(p_period_all[[19]], p_period_all[[20]], p_period_all[[21]],
          p_period_all[[22]], p_period_all[[23]], p_period_all[[24]],
          ncol = 2)
ggsave(filename = "output/plots/PA/period_example_4.png",
       width = 17, height = 21, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/PA/period_example_4.pdf",
       width = 17, height = 21, units = "cm", dpi = 600)
plot_grid(p_period_all[[25]], p_period_all[[26]], p_period_all[[28]],
          p_period_all[[29]],
          ncol = 2)
ggsave(filename = "output/plots/PA/period_example_5.png",
       width = 17, height = 14, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/PA/period_example_5.pdf",
       width = 17, height = 14, units = "cm", dpi = 600)


### run MSE projections to get time series
# req_pckgs <- c("FLCore", "FLash", "mse", "GA", "doParallel", "doRNG", "FLBRP")
# for (i in req_pckgs) library(package = i, character.only = TRUE)
# source("funs.R")
# cl1 <- makeCluster(5)
# registerDoParallel(cl1)
# cl_length_1 <- length(cl1)
# . <- foreach(i = seq(cl_length_1)) %dopar% {
#   for (i in req_pckgs) library(package = i, character.only = TRUE,
#                                warn.conflicts = FALSE, verbose = FALSE,
#                                quietly = TRUE)
#   source("funs.R", echo = FALSE)
# }
# . <- foreach(scn = split(mult_periods_5, seq(nrow(mult_periods_5)))) %dopar% {
#   input <- readRDS(paste0("input/500_50/OM_2_mp_input/", scn$fhist, "/",
#                           scn$stock, ".rds"))
#   input$ctrl$est@args$comp_m <- scn$multiplier
#   res <- do.call(mp, input)
#   saveRDS(res, file = paste0("output/500_50/PA/", scn$fhist, "/", scn$stock, 
#                              "/mp_1_2_3_1_1_1_1_2_", scn$multiplier, 
#                              "_Inf_0.rds" ))
#   qnts <- collapse_correction(res@stock, yrs = 100:150)
#   saveRDS(qnts, file = paste0("output/500_50/PA/", scn$fhist, "/", scn$stock, 
#                               "/qnts_1_2_3_1_1_1_1_2_", scn$multiplier, 
#                               "_Inf_0.rds" ))
#   return(NULL)
# }

### plot SSB trajectories for multipliers where risk = 5%
# fhist <- "one-way"
# stock <- "pol"
# 
# qnts_def <- mult_periods_5 %>%
#   filter(stock == !!stock & fhist == !!fhist)
# qnts <- foreach(i = split(qnts_def, seq(nrow(qnts_def))),
#                 .combine = bind_rows) %do% {
#   tmp <- readRDS(paste0("output/500_50/PA/", i$fhist, "/",
#                         i$stock, "/qnts_1_2_3_1_1_1_1_2_", i$multiplier,
#                         "_Inf_0.rds"))$ssb
#   brp <- brps[[i$stock]]
#   tmp <- tmp/c(refpts(brp)["msy", "ssb"])
#   tmp <- quantile(tmp, probs = c(0.05, 0.5, 0.95))
#   tmp <- as.data.frame(tmp) %>%
#     select(year, iter, data) %>%
#     mutate(stock = i$stock, fhist = i$fhist, multiplier = i$multiplier)
#   tmp <- tmp %>% full_join(qnts_def)
#   tmp <- tmp %>%
#     pivot_wider(names_from = iter, values_from = data)
#   return(tmp)
# }
# BlimBmsy <- c(refpts(brp)["virgin", "ssb"]/refpts(brp)["msy", "ssb"]) * 0.1628
# 
# p_pol_ssb <- qnts %>%
#   ggplot(aes(x = year - 100)) +
#   geom_line(aes(y = `50%`, colour = period, linetype = period), 
#             size = 0.4, show.legend = FALSE) +
#   geom_line(aes(y = `5%`, colour = period, linetype = period), 
#             size = 0.05, show.legend = FALSE) +
#   geom_line(aes(y = `95%`, colour = period, linetype = period), 
#             size = 0.05, show.legend = FALSE) +
#   geom_ribbon(aes(ymin = `5%`, ymax = `95%`, fill = period), alpha = 0.1,
#               show.legend = FALSE) +
#   geom_vline(xintercept = 0, size = 0.3) +
#   geom_vline(xintercept = 10, size = 0.3) +
#   geom_vline(xintercept = 40, size = 0.3) +
#   geom_vline(xintercept = 50, size = 0.3) +
#   geom_hline(yintercept = BlimBmsy, colour = "red", size = 0.4) +
#   scale_colour_manual("",
#     values = c("first 10 years" = brewer.pal(4, name = "Set1")[2],
#                "last 10 years" = brewer.pal(4, name = "Set1")[3],
#                "all years" = brewer.pal(4, name = "Set1")[4]),
#     guide = guide_legend(order = 1)) +
#   scale_fill_manual("",
#     values = c("first 10 years" = brewer.pal(4, name = "Set1")[2],
#                "last 10 years" = brewer.pal(4, name = "Set1")[3],
#                "all years" = brewer.pal(4, name = "Set1")[4]),
#     guide = guide_legend(order = 1)) +
#   scale_linetype_manual("", 
#     values = c("first 10 years" = "solid",
#                "last 10 years" = "dotted",
#                "all years" = "dashed"),
#     guide = guide_legend(order = 1)) +
#   facet_wrap(~ stock_label, strip.position = "right", labeller = label_parsed) +
#   labs(y = expression(SSB/italic(B)[MSY]), 
#        x = "year") +
#   xlim(0, 50) + #ylim(0, 1)
#   theme_bw(base_size = 8, base_family = "sans") +
#   theme()
# p_pol_ssb

### ------------------------------------------------------------------------ ###
### plot - pollack rfb-rule component explorations ####
### ------------------------------------------------------------------------ ###

pol <- readRDS("output/pol_PA_components_stats.rds")
### create labels for figure
pol_plot_stats <- pol %>%
  mutate(label = factor(pol$optimised,
                        levels = c("default", "mult", "cap", "mult_cap", "all",
                                   "all_cap"),
                        labels = c("default\n(not optimised)", 
                                   "multiplier", 
                                   "uncertainty\ncap", 
                                   "multiplier\nand cap",
                                   "all without\ncap",
                                   "all"))) %>%
  select(fhist, label, SSB_rel, Catch_rel, risk_Blim, ICV, fitness)
### repeat, but add annotations to labels
pol_plot_fitness <- pol %>%
  mutate(label = factor(pol$optimised,
                        levels = c("default", "mult", "cap", "mult_cap", "all",
                                   "all_cap"),
                        labels = c("default\n(not optimised)*", 
                                   "multiplier*", 
                                   "uncertainty\ncap", 
                                   "multiplier\nand cap",
                                   "all without\ncap",
                                   "all*"))) %>%
  select(fhist, label, SSB_rel, Catch_rel, risk_Blim, ICV, fitness)

### recreate fitness elements
pol_fitness <- pol_plot_fitness %>%
  mutate(comp_Catch = Catch_rel - 1,
         comp_SSB = SSB_rel - 1,
         comp_ICV = ICV, 
         comp_risk_penalty = penalty(x = risk_Blim, 
                                     negative = FALSE, max = 5, 
                                     inflection = 0.05 + 0.01, 
                                     steepness = 0.5e+3)) %>%
  #mutate(Catch_rel = NULL, SSB_rel = NULL, ICV = NULL, risk_Blim = NULL) %>%
  pivot_longer(c(comp_Catch, comp_SSB, comp_ICV, comp_risk_penalty),
               names_prefix = "comp_") %>%
  mutate(name = factor(name,
                       levels = rev(c("SSB", "Catch", "ICV", 
                                      "risk_penalty")),
                       labels = rev(c("SSB", "Catch", "ICV", 
                                      "risk-PA\n(penalty)"))))
pol_fitness_dev <- pol_fitness %>%
  filter(name %in% c("SSB", "Catch")) %>%
  mutate(value_sign = ifelse(name == "SSB",
                             -abs(SSB_rel - 1)/2,
                             -abs(SSB_rel - 1) - abs(Catch_rel - 1)/2)) %>%
  mutate(sign = ifelse(name == "SSB",
                       ifelse(SSB_rel > 1, "+", "-"),
                       ifelse(Catch_rel > 1, "+", "-"))) %>%
  filter(abs(value) >= 0.02)

### plot fitness function, split into elements
p_pol_fitness <- pol_fitness %>% 
  mutate(value = -abs(value)) %>%
  ggplot(aes(x = label, y = value, fill = name)) +
  geom_hline(yintercept = 0, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = "stack", width = 0.8, 
           colour = "black", size = 0.1) +
  ### for reversing order in legend (not in plot)
  scale_fill_discrete("fitness elements",
                      breaks = rev(levels(pol_fitness$name))) +
  geom_text(data = pol_fitness_dev,
            aes(x = label, y = value_sign, label = sign), 
            vjust = 0.5, colour = "grey20") +
  facet_grid("fitness" ~ fhist, scales = "free", space = "free_x", switch = "y",
             labeller = "label_parsed") +
  labs(y = "", x = "") +
  theme_bw(base_size = 8, base_family = "sans") +
  theme(panel.spacing.x = unit(0, units = "cm"),
        strip.placement.y = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5,
                                   lineheight = 0.7),
        plot.margin = unit(x = c(1, 3, 3, 3), units = "pt"),
        axis.title.y = element_blank(),
        legend.key.width = unit(1, "lines"),
        legend.key.height = unit(1, "lines"))

### plot stats individually
p_pol_stats_SSB <- pol_plot_stats %>% 
  ggplot(aes(x = label, y = SSB_rel, fill = label)) +
  geom_hline(yintercept = 1, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = position_dodge2(preserve = "single"), width = 0.8, 
           show.legend = FALSE, colour = "black", size = 0.1) +
  scale_fill_grey() +
  facet_grid("SSB/B[MSY]" ~ fhist, scales = "free", space = "free_x", switch = "y",
             labeller = "label_parsed") +
  labs(y = "", x = "fitness function") +
  theme_bw(base_size = 8, base_family = "sans") +
  theme(panel.spacing.x = unit(0, units = "cm"),
        strip.placement.y = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        plot.margin = unit(x = c(1, 3, 0, 3), units = "pt"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_y_continuous(trans = trans_from(), limits = c(0, 3.25))
p_pol_stats_C <- pol_plot_stats %>% 
  ggplot(aes(x = label, y = Catch_rel, fill = label)) +
  geom_hline(yintercept = 1, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = "dodge", show.legend = FALSE, width = 0.8, 
           colour = "black", size = 0.1) +
  scale_fill_grey() +
  facet_grid("Catch/MSY" ~ fhist, scales = "free", space = "free_x", switch = "y",
             labeller = "label_parsed") +
  labs(y = "", x = "fitness function") +
  theme_bw(base_size = 8, base_family = "sans") +
  theme(panel.spacing.x = unit(0, units = "cm"),
        strip.text.x = element_blank(),
        strip.placement.y = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        plot.margin = unit(x = c(0, 3, 0, 3), units = "pt"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_y_continuous(trans = trans_from(), limits = c(0, 3.25))
p_pol_stats_risk <- pol_plot_stats %>% 
  ggplot(aes(x = label, y = risk_Blim, fill = label)) +
  geom_hline(yintercept = 0, linetype = "solid", size = 0.5, colour = "grey") +
  geom_hline(yintercept = 0.05, linetype = "solid", size = 0.5, colour = "red") +
  geom_col(position = "dodge", show.legend = FALSE, width = 0.8, 
           colour = "black", size = 0.1) +
  scale_fill_grey() +
  facet_grid("B[lim]~risk" ~ fhist, scales = "free", space = "free_x", switch = "y",
             labeller = "label_parsed") +
  labs(y = "", x = "fitness function") +
  theme_bw(base_size = 8, base_family = "sans") +
  theme(panel.spacing.x = unit(0, units = "cm"),
        strip.text.x = element_blank(),
        strip.placement.y = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        plot.margin = unit(x = c(0, 3, 0, 3), units = "pt"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_y_continuous(trans = trans_from(0), limits = c(0, 1))
p_pol_stats_ICV <- pol_plot_stats %>% 
  ggplot(aes(x = label, y = ICV, fill = label)) +
  geom_hline(yintercept = 0, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = "dodge", show.legend = FALSE, width = 0.8, 
           colour = "black", size = 0.1) +
  scale_fill_grey() +
  facet_grid("ICV" ~ fhist, scales = "free", space = "free_x", switch = "y",
             labeller = "label_parsed") +
  labs(y = "", x = "\nrfb-rule parameters included in optimisation") +
  theme_bw(base_size = 8, base_family = "sans") +
  theme(panel.spacing.x = unit(0, units = "cm"),
        strip.text.x = element_blank(),
        strip.placement.y = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5,
                                   lineheight = 0.7),
        plot.margin = unit(x = c(0, 3, 3, 3), units = "pt"),
        axis.title.y = element_blank()) +
  scale_y_continuous(trans = trans_from(0), limits = c(0, 1))
p_pol_stats_comb <- 
  plot_grid(plot_grid(p_pol_stats_SSB, p_pol_stats_C,
                      p_pol_stats_risk, p_pol_stats_ICV,
                      ncol = 1, align = "v",
                      rel_heights = c(1.25, 1, 1, 2)), 
            plot_grid(p_pol_fitness + theme(legend.position = "none"), 
                      plot_grid(NULL, get_legend(p_pol_fitness), 
                                ncol = 2, rel_widths = c(1, 0.45)),
                      ncol = 1, 
                      rel_heights = c(1, 0.4)),
            ncol = 2, labels = c("(a)", "(b)"), label_size = 10)
p_pol_stats_comb
ggsave(filename = "output/plots/PA/pol_components_stats.png", 
       plot = p_pol_stats_comb,
       width = 17, height = 11, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/PA/pol_components_stats.pdf", 
       plot = p_pol_stats_comb,
       width = 17, height = 11, units = "cm", dpi = 600)

### ------------------------------------------------------------------------ ###
### plot - all stocks stats PA - default vs. optimised - not used ####
### ------------------------------------------------------------------------ ###

### load data
all_GA <- readRDS("output/all_stocks_GA_optimised_stats.rds")
stocks_sorted <- stocks %>%
  select(stock, k) %>%
  arrange(k)

### format for plotting
stats_plot <- all_GA %>% 
  filter(scenario == "PA" & optimised != "all") %>%
  mutate(ICV = ifelse(optimised == "zero-fishing", 0, ICV)) %>%
  mutate(stock = factor(stock, levels = stocks_sorted$stock)) %>%
  pivot_longer(c(SSB_rel, Fbar_rel, Catch_rel, risk_Blim, ICV, 
                 fitness)) %>%
  mutate(stat = name) %>%
  mutate(stat = ifelse(stat == "SSB_rel", "SSB/B[MSY]", stat),
         stat = ifelse(stat == "Fbar_rel", "F/F[MSY]", stat),
         stat = ifelse(stat == "Catch_rel", "Catch/MSY", stat),
         stat = ifelse(stat == "risk_Blim", "B[lim]~risk", stat),
         stat = ifelse(stat == "ICV", "ICV", stat),
         stat = ifelse(stat == "fitness", "fitness~value", 
                       stat)) %>%
  mutate(stat = factor(stat, levels = unique(stat)[c(1, 2, 3, 4, 5, 6)]),
         rule = factor(optimised, 
                       levels = c("zero-fishing", "default", "mult", "all_cap"),
                       labels = c("zero fishing", "not optimised", 
                                  "GA: multiplier", "GA: all")))
saveRDS(stats_plot, file = "output/plots/PA/data_stocks_stats.rds")
stats_plot <- readRDS("output/plots/PA/data_stocks_stats.rds")

### remove zero fishing
stats_plot <- stats_plot %>%
  filter(rule != "zero fishing")
### individual plots
p_stats_SSB <- stats_plot %>%
  filter(stat %in% c("SSB/B[MSY]")) %>%
  ggplot(aes(x = stock, y = value, fill = rule, colour = rule)) +
  geom_hline(yintercept = 1, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = position_dodge2(), width = 0.8, 
           show.legend = FALSE, size = 0.1) +
  scale_fill_manual(values = c("not optimised" = "#1B9E77", 
                               "GA: multiplier" =  "#D95F02", 
                               "GA: all" = "#7570B3")) +
  scale_colour_manual(values = c("not optimised" = "#1B9E77", 
                                 "GA: multiplier" =  "#D95F02", 
                                 "GA: all" = "#7570B3")) +
  facet_grid(stat ~ fhist, 
             switch = "y", scales = "free", labeller = "label_parsed") +
  labs(x = "stock", y = "") +
  theme_bw(base_size = 8, base_family = "sans") +
  theme(strip.placement.y = "outside",
        strip.background.y = element_blank(), 
        panel.spacing.x = unit(0, units = "cm"),
        strip.text.y = element_text(size = 8),
        axis.title.y = element_blank(),
        plot.margin = unit(x = c(1, 3, 0, 3), units = "pt"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) +
  scale_y_continuous(trans = trans_from(), limits = c(0, NA))
p_stats_F <- stats_plot %>%
  filter(stat %in% c("F/F[MSY]")) %>%
  ggplot(aes(x = stock, y = value, fill = rule, colour = rule)) +
  geom_hline(yintercept = 1, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = position_dodge2(), width = 0.8, 
           show.legend = FALSE, size = 0.1) +
  scale_fill_manual(values = c("not optimised" = "#1B9E77", 
                               "GA: multiplier" =  "#D95F02", 
                               "GA: all" = "#7570B3")) +
  scale_colour_manual(values = c("not optimised" = "#1B9E77", 
                                 "GA: multiplier" =  "#D95F02", 
                                 "GA: all" = "#7570B3")) +
  facet_grid(stat ~ fhist, 
             switch = "y", scales = "free", labeller = "label_parsed") +
  labs(x = "stock", y = "") +
  theme_bw(base_size = 8, base_family = "sans") +
  theme(strip.placement.y = "outside",
        strip.background.y = element_blank(), 
        panel.spacing.x = unit(0, units = "cm"),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 8),
        axis.title.y = element_blank(),
        plot.margin = unit(x = c(1, 3, 0, 3), units = "pt"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) +
  scale_y_continuous(trans = trans_from(), limits = c(0, NA))
p_stats_C <- stats_plot %>%
  filter(stat %in% c("Catch/MSY")) %>%
  ggplot(aes(x = stock, y = value, fill = rule, colour = rule)) +
  geom_hline(yintercept = 1, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = position_dodge2(), width = 0.8, 
           show.legend = FALSE, size = 0.1) +
  scale_fill_manual(values = c("not optimised" = "#1B9E77", 
                               "GA: multiplier" =  "#D95F02", 
                               "GA: all" = "#7570B3")) +
  scale_colour_manual(values = c("not optimised" = "#1B9E77", 
                                 "GA: multiplier" =  "#D95F02", 
                                 "GA: all" = "#7570B3")) +
  facet_grid(stat ~ fhist, 
             switch = "y", scales = "free", labeller = "label_parsed") +
  labs(x = "stock", y = "") +
  theme_bw(base_size = 8, base_family = "sans") +
  theme(strip.placement.y = "outside",
        strip.background.y = element_blank(), 
        panel.spacing.x = unit(0, units = "cm"),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 8),
        axis.title.y = element_blank(),
        plot.margin = unit(x = c(1, 3, 0, 3), units = "pt"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) +
  scale_y_continuous(trans = trans_from(), limits = c(0, NA))
p_stats_risk <- stats_plot %>%
  filter(stat %in% c("B[lim]~risk")) %>%
  ggplot(aes(x = stock, y = value, fill = rule, colour = rule)) +
  geom_hline(yintercept = 0, linetype = "solid", size = 0.5, colour = "grey") +
  geom_hline(yintercept = 0.05, linetype = "solid", size = 0.5, colour = "red") +
  geom_col(position = position_dodge2(), width = 0.8, 
           show.legend = FALSE, size = 0.1) +
  scale_fill_manual(values = c("not optimised" = "#1B9E77", 
                               "GA: multiplier" =  "#D95F02", 
                               "GA: all" = "#7570B3")) +
  scale_colour_manual(values = c("not optimised" = "#1B9E77", 
                                 "GA: multiplier" =  "#D95F02", 
                                 "GA: all" = "#7570B3")) +
  facet_grid(stat ~ fhist, 
             switch = "y", scales = "free", labeller = "label_parsed") +
  labs(x = "stock", y = "") +
  theme_bw(base_size = 8, base_family = "sans") +
  theme(strip.placement.y = "outside",
        strip.background.y = element_blank(), 
        panel.spacing.x = unit(0, units = "cm"),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 8),
        axis.title.y = element_blank(),
        plot.margin = unit(x = c(1, 3, 0, 3), units = "pt"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) +
  scale_y_continuous(trans = trans_from(0), limits = c(0, 1))
p_stats_ICV <- stats_plot %>%
  filter(stat %in% c("ICV")) %>%
  ggplot(aes(x = stock, y = value, fill = rule, colour = rule)) +
  geom_hline(yintercept = 0, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = position_dodge2(), width = 0.8, 
           show.legend = FALSE, size = 0.1) +
  scale_fill_manual(values = c("not optimised" = "#1B9E77", 
                               "GA: multiplier" =  "#D95F02", 
                               "GA: all" = "#7570B3")) +
  scale_colour_manual(values = c("not optimised" = "#1B9E77", 
                                 "GA: multiplier" =  "#D95F02", 
                                 "GA: all" = "#7570B3")) +
  facet_grid(stat ~ fhist, 
             switch = "y", scales = "free", labeller = "label_parsed") +
  labs(x = "stock", y = "") +
  theme_bw(base_size = 8, base_family = "sans") +
  theme(strip.placement.y = "outside",
        strip.background.y = element_blank(), 
        panel.spacing.x = unit(0, units = "cm"),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 8),
        axis.title.y = element_blank(),
        plot.margin = unit(x = c(1, 3, 0, 3), units = "pt"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) +
  scale_y_continuous(trans = trans_from(0), limits = c(0, 1))
p_stats_fitness <- stats_plot %>%
  filter(stat %in% c("fitness~value")) %>%
  ggplot(aes(x = stock, y = value, fill = rule, colour = rule)) +
  geom_hline(yintercept = 0, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = position_dodge2(), width = 0.8, 
           size = 0.1, show.legend = TRUE) +
  scale_fill_manual("", values = c("not optimised" = "#1B9E77", 
                               "GA: multiplier" =  "#D95F02", 
                               "GA: all" = "#7570B3")) +
  scale_colour_manual("", values = c("not optimised" = "#1B9E77", 
                                 "GA: multiplier" =  "#D95F02", 
                                 "GA: all" = "#7570B3")) +
  facet_grid(stat ~ fhist, 
             switch = "y", scales = "free", labeller = "label_parsed") +
  labs(x = "stock", y = "") +
  theme_bw(base_size = 8, base_family = "sans") +
  theme(strip.placement.y = "outside",
        strip.background.y = element_blank(), 
        panel.spacing.x = unit(0, units = "cm"),
        plot.margin = unit(x = c(0, 3, 3, 3), units = "pt"),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 8),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.title.y = element_blank(),
        legend.key.width = unit(0.2, units = "lines"),
        legend.key.height = unit(1, units = "lines")) +
  scale_y_continuous(trans = trans_from(0), limits = c(NA, NA))
p_stats_combined <- plot_grid(
  plot_grid(p_stats_SSB, 
            #p_stats_F, 
            p_stats_C,
            p_stats_risk, p_stats_ICV,
            p_stats_fitness + theme(legend.position = "none"),
            ncol = 1, align = "v",
            rel_heights = c(1.25, 1, 1, 1, 1.5)),
  get_legend(p_stats_fitness), rel_widths = c(1, 0.15))
p_stats_combined
ggsave(filename = "output/plots/PA/all_stocks_stats.png", 
       plot = p_stats_combined,
       width = 17, height = 10, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/PA/all_stocks_stats.pdf", 
       plot = p_stats_combined,
       width = 17, height = 10, units = "cm", dpi = 600)

### ------------------------------------------------------------------------ ###
### plot - all stocks stats PA vs. 2 over 3 - not used ####
### ------------------------------------------------------------------------ ###

### load data
all_GA <- readRDS("output/all_stocks_GA_optimised_stats.rds")

stats_2over3 <- foreach(stock = stocks$stock, .combine = "rbind") %:%
  foreach(fhist = c("one-way", "random"), .combine = "rbind") %do% {
    stats_i <- readRDS(paste0("output/500_50/2over3/", fhist, "/", stock, 
                              "_stats.rds"))
    stats_i <- as.data.frame(lapply(as.data.frame(t(stats_i)), unlist))
    stats_i$fhist <- fhist
    stats_i$stock = stock
    stats_i$catch_rule = "2 over 3"
    stats_i$group = "2 over 3"
    ### generate fitness
    stats_i_MSY <- stats_i_PA <- stats_i
    stats_i_MSY$fitness <- -sum(abs(stats_i_MSY$SSB_rel - 1),
                                abs(stats_i_MSY$Catch_rel - 1),
                                stats_i_MSY$ICV,
                                stats_i_MSY$risk_Blim)
    stats_i_PA$fitness <- -sum(abs(stats_i_PA$SSB_rel - 1),
                               abs(stats_i_PA$Catch_rel - 1),
                               stats_i_PA$ICV,
                               penalty(x = stats_i_PA$risk_Blim, negative = FALSE, 
                                       max = 5, inflection = 0.06, 
                                       steepness = 0.5e+3))
    stats_i_MSY$scenario <- "MSY"
    stats_i_PA$scenario <- "PA"
    stats_i_MSY$optimised <- stats_i_PA$optimised <- "default"
    return(rbind(stats_i_MSY, stats_i_PA))
}
saveRDS(stats_2over3, "output/all_stocks_2over_stats.rds")
stats_2over3 <- readRDS("output/all_stocks_2over_stats.rds")

### combine 
stats_plot <- bind_rows(
  all_GA %>% 
    filter(optimised == "all" & scenario == "MSY") %>%
    select(fhist, stock, SSB_rel, Fbar_rel, Catch_rel, risk_Blim, ICV) %>%
    mutate(catch_rule = "rfb", group = "rfb: MSY"),
  all_GA %>% 
    filter(optimised == "all_cap" & scenario == "PA") %>%
    select(fhist, stock, SSB_rel, Fbar_rel, Catch_rel, risk_Blim, ICV) %>%
    mutate(catch_rule = "rfb", group = "rfb: MSY-PA"),
  stats_2over3 %>%
    filter(scenario == "PA") %>%
    select(fhist, stock, SSB_rel, Fbar_rel, Catch_rel, risk_Blim, ICV) %>%
    mutate(catch_rule = "2 over 3", group = "2 over 3")
)

stocks_sorted <- stocks %>%
  select(stock, k) %>%
  arrange(k)

### format for plotting
stats_plot <- stats_plot %>% 
  mutate(stock = factor(stock, levels = stocks_sorted$stock)) %>%
  pivot_longer(c(SSB_rel, Fbar_rel, Catch_rel, risk_Blim, ICV)) %>%
  mutate(stat = name) %>%
  mutate(stat = ifelse(stat == "SSB_rel", "SSB/B[MSY]", stat),
         stat = ifelse(stat == "Fbar_rel", "F/F[MSY]", stat),
         stat = ifelse(stat == "Catch_rel", "Catch/MSY", stat),
         stat = ifelse(stat == "risk_Blim", "B[lim]~risk", stat),
         stat = ifelse(stat == "ICV", "ICV", stat)) %>%
  mutate(stat = factor(stat, levels = unique(stat)[c(1, 2, 3, 4, 5)]),
         group = factor(group, 
                        levels = c("2 over 3", "rfb: MSY", "rfb: MSY-PA"),
                        labels = c("2 over 3", "rfb: MSY", "rfb: MSY-PA")))
saveRDS(stats_plot, file = "output/plots/PA/data_stocks_2over3_stats.rds")
stats_plot <- readRDS("output/plots/PA/data_stocks_2over3_stats.rds")

### individual plots
p_stats_SSB <- stats_plot %>%
  filter(stat %in% c("SSB/B[MSY]")) %>%
  ggplot(aes(x = stock, y = value, fill = group, colour = group)) +
  geom_hline(yintercept = 1, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = position_dodge2(), width = 0.8, 
           show.legend = FALSE, size = 0.1) +
  scale_fill_manual(values = c("2 over 3" = "#E6AB02", 
                               "rfb: MSY" = "#66A61E", 
                               "rfb: MSY-PA" =  "#7570B3")) +
  scale_colour_manual(values = c("2 over 3" = "#E6AB02", 
                                 "rfb: MSY" = "#66A61E", 
                                 "rfb: MSY-PA" =  "#7570B3")) +
  facet_grid(stat ~ fhist, 
             switch = "y", scales = "free", labeller = "label_parsed") +
  labs(x = "stock", y = "") +
  theme_bw(base_size = 8, base_family = "sans") +
  theme(strip.placement.y = "outside",
        strip.background.y = element_blank(), 
        panel.spacing.x = unit(0, units = "cm"),
        strip.text.y = element_text(size = 8),
        axis.title.y = element_blank(),
        plot.margin = unit(x = c(1, 3, 0, 3), units = "pt"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) +
  scale_y_continuous(trans = trans_from(), limits = c(0, NA))
p_stats_F <- stats_plot %>%
  filter(stat %in% c("F/F[MSY]")) %>%
  ggplot(aes(x = stock, y = value, fill = group, colour = group)) +
  geom_hline(yintercept = 1, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = position_dodge2(), width = 0.8, 
           show.legend = FALSE, size = 0.1) +
  scale_fill_manual(values = c("2 over 3" = "#E6AB02", 
                               "rfb: MSY" = "#66A61E", 
                               "rfb: MSY-PA" =  "#7570B3")) +
  scale_colour_manual(values = c("2 over 3" = "#E6AB02", 
                                 "rfb: MSY" = "#66A61E", 
                                 "rfb: MSY-PA" =  "#7570B3")) +
  facet_grid(stat ~ fhist, 
             switch = "y", scales = "free", labeller = "label_parsed") +
  labs(x = "stock", y = "") +
  theme_bw(base_size = 8, base_family = "sans") +
  theme(strip.placement.y = "outside",
        strip.background.y = element_blank(), 
        panel.spacing.x = unit(0, units = "cm"),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 8),
        axis.title.y = element_blank(),
        plot.margin = unit(x = c(1, 3, 0, 3), units = "pt"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) +
  scale_y_continuous(trans = trans_from(), limits = c(0, NA))
p_stats_C <- stats_plot %>%
  filter(stat %in% c("Catch/MSY")) %>%
  ggplot(aes(x = stock, y = value, fill = group, colour = group)) +
  geom_hline(yintercept = 1, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = position_dodge2(), width = 0.8, 
           show.legend = FALSE, size = 0.1) +
  scale_fill_manual(values = c("2 over 3" = "#E6AB02", 
                               "rfb: MSY" = "#66A61E", 
                               "rfb: MSY-PA" =  "#7570B3")) +
  scale_colour_manual(values = c("2 over 3" = "#E6AB02", 
                                 "rfb: MSY" = "#66A61E", 
                                 "rfb: MSY-PA" =  "#7570B3")) +
  facet_grid(stat ~ fhist, 
             switch = "y", scales = "free", labeller = "label_parsed") +
  labs(x = "stock", y = "") +
  theme_bw(base_size = 8, base_family = "sans") +
  theme(strip.placement.y = "outside",
        strip.background.y = element_blank(), 
        panel.spacing.x = unit(0, units = "cm"),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 8),
        axis.title.y = element_blank(),
        plot.margin = unit(x = c(1, 3, 0, 3), units = "pt"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) +
  scale_y_continuous(trans = trans_from(), limits = c(0, NA))
p_stats_risk <- stats_plot %>%
  filter(stat %in% c("B[lim]~risk")) %>%
  ggplot(aes(x = stock, y = value, fill = group, colour = group)) +
  geom_hline(yintercept = 0, linetype = "solid", size = 0.5, colour = "grey") +
  geom_hline(yintercept = 0.05, linetype = "solid", size = 0.5, colour = "red") +
  geom_col(position = position_dodge2(), width = 0.8, 
           show.legend = FALSE, size = 0.1) +
  scale_fill_manual(values = c("2 over 3" = "#E6AB02", 
                               "rfb: MSY" = "#66A61E", 
                               "rfb: MSY-PA" =  "#7570B3")) +
  scale_colour_manual(values = c("2 over 3" = "#E6AB02", 
                                 "rfb: MSY" = "#66A61E", 
                                 "rfb: MSY-PA" =  "#7570B3")) +
  facet_grid(stat ~ fhist, 
             switch = "y", scales = "free", labeller = "label_parsed") +
  labs(x = "stock", y = "") +
  theme_bw(base_size = 8, base_family = "sans") +
  theme(strip.placement.y = "outside",
        strip.background.y = element_blank(), 
        panel.spacing.x = unit(0, units = "cm"),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 8),
        axis.title.y = element_blank(),
        plot.margin = unit(x = c(1, 3, 0, 3), units = "pt"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) +
  scale_y_continuous(trans = trans_from(0), limits = c(0, 1))
p_stats_ICV <- stats_plot %>%
  filter(stat %in% c("ICV")) %>%
  ggplot(aes(x = stock, y = value, fill = group, colour = group)) +
  geom_hline(yintercept = 0, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = position_dodge2(), width = 0.8, 
           show.legend = TRUE, size = 0.1) +
  scale_fill_manual("catch rule",
                    values = c("2 over 3" = "#E6AB02", 
                               "rfb: MSY" = "#66A61E", 
                               "rfb: MSY-PA" =  "#7570B3")) +
  scale_colour_manual("catch rule",
                      values = c("2 over 3" = "#E6AB02", 
                                 "rfb: MSY" = "#66A61E", 
                                 "rfb: MSY-PA" =  "#7570B3")) +
  # scale_fill_discrete("catch rule") + scale_colour_discrete("catch rule") +
  facet_grid(stat ~ fhist, 
             switch = "y", scales = "free", labeller = "label_parsed") +
  labs(x = "stock", y = "") +
  theme_bw(base_size = 8, base_family = "sans") +
  theme(strip.placement.y = "outside",
        strip.background.y = element_blank(), 
        panel.spacing.x = unit(0, units = "cm"),
        plot.margin = unit(x = c(0, 3, 3, 3), units = "pt"),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 8),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.title.y = element_blank(),
        legend.key.width = unit(0.2, units = "lines"),
        legend.key.height = unit(1, units = "lines")) +
  scale_y_continuous(trans = trans_from(0), limits = c(NA, 1))
p_stats_combined <- plot_grid(
  plot_grid(p_stats_SSB, 
            #p_stats_F, 
            p_stats_C,
            p_stats_risk, 
            p_stats_ICV + theme(legend.position = "none"),
            ncol = 1, align = "v",
            rel_heights = c(1.25, 1, 1, 1.5)),
  get_legend(p_stats_ICV), rel_widths = c(1, 0.15))
p_stats_combined
ggsave(filename = "output/plots/PA/all_stocks_2over3_stats.png", 
       plot = p_stats_combined,
       width = 17, height = 10, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/PA/all_stocks_2over3_stats.pdf", 
       plot = p_stats_combined,
       width = 17, height = 10, units = "cm", dpi = 600)

### ------------------------------------------------------------------------ ###
### plot comparison of rules & optimisations ####
### ------------------------------------------------------------------------ ###

### stats from 2 over 3 rule
stats_2over3 <- readRDS("output/all_stocks_2over_stats.rds")

### stats from rfb-rule
stats_rfb <- readRDS("output/all_stocks_GA_optimised_stats.rds")

### combine 
stats_plot <- bind_rows(
  ### zero fishing
  stats_rfb %>% 
    filter(capped == FALSE, optimised == "zero-fishing" & multiplier == 0 & 
             scenario == "PA") %>%
    mutate(catch_rule = "zero-fishing", group = "zero-fishing"),
  ### 2 over 3 rule
  stats_2over3 %>%
    filter(scenario == "PA") %>%
    mutate(catch_rule = "2 over 3", group = "2 over 3"),
  ### rfb-rule, default
  stats_rfb %>% 
    filter(capped == FALSE, optimised == "default" & scenario == "PA") %>%
    mutate(catch_rule = "rfb", target = "none", group = "rfb: default"),
  ### rfb-rule, optimisation with multiplier, target MSY
  stats_rfb %>% 
    filter(capped == FALSE, optimised == "mult" & scenario == "MSY") %>%
    mutate(catch_rule = "rfb", target = "MSY", group = "rfb: MSY - mult"),
  ### rfb-rule, optimisation with all parameters, target MSY
  stats_rfb %>% 
    filter(capped == FALSE, optimised == "all" & scenario == "MSY") %>%
    mutate(catch_rule = "rfb", target = "MSY", group = "rfb: MSY - all"),
  ### rfb-rule, optimisation with multiplier, target PA
  stats_rfb %>% 
    filter(capped == FALSE, optimised == "mult" & scenario == "PA") %>%
    mutate(catch_rule = "rfb", target = "PA", group = "rfb: PA - mult"),
  ### rfb-rule, optimisation with all parameters, target PA
  stats_rfb %>% 
    filter(capped == FALSE, optimised == "all_cap" & scenario == "PA") %>%
    mutate(catch_rule = "rfb", target = "PA", group = "rfb: PA - all"),
  ### rfb-rule, always capped, optimisation with multiplier, target PA
  stats_rfb %>% 
    filter(capped == TRUE, optimised == "mult" & scenario == "PA") %>%
    mutate(catch_rule = "rfb", target = "PA", 
           group = "rfb (capped):\nPA - mult"),
  ### rfb-rule, capped below Itrigger, optimisation with multiplier, target PA
  stats_rfb %>% 
    filter(capped == TRUE, optimised == "mult" & scenario == "PA_capped") %>%
    mutate(catch_rule = "rfb", target = "PA", 
           group = "rfb (cond. capped):\nPA - mult"),
  ### rfb-rule, capped below Itrigger, optimisation with all, target PA
  stats_rfb %>% 
    filter(capped == TRUE, optimised == "all" & scenario == "PA_capped") %>%
    mutate(catch_rule = "rfb", target = "PA", 
           group = "rfb (cond. capped):\nPA - all")
)
stats_plot <- stats_plot %>%
  select(fhist, stock, catch_rule, target, group, 
         SSB_rel, Catch_rel, ICV, risk_Blim) %>%
  mutate(penalty = penalty(x = risk_Blim, 
                           negative = FALSE, max = 5, 
                           inflection = 0.06, 
                           steepness = 0.5e+3)) %>%
  mutate(comp_Catch = Catch_rel - 1,
         comp_SSB = SSB_rel - 1,
         comp_ICV = ICV,
         comp_risk_penalty =  penalty) %>%
  pivot_longer(c(comp_SSB, comp_Catch, comp_ICV, comp_risk_penalty),
               names_prefix = "comp_") %>%
  mutate(name = factor(name,
                       levels = rev(c("SSB", "Catch", "ICV",
                                      "risk_penalty")),
                       labels = rev(c("SSB", "Catch", "ICV",
                                      "risk-PA\n(penalty)")))) %>% 
  full_join(stocks %>%
    select(stock, k) %>%
    mutate(stock = factor(stock, levels = stock),
           stock_k = paste0(stock, "~(italic(k)==", sprintf(k, fmt =  "%.2f"), 
                            "*year^-1)")) %>%
    mutate(stock_k = factor(stock_k, levels = stock_k))) %>%
  mutate(group = factor(group, 
    levels = c("zero-fishing", "2 over 3", 
               "rfb: default",
               "rfb: MSY - mult", "rfb: MSY - all",
               "rfb: PA - mult", "rfb: PA - all",
               "rfb (capped):\nPA - mult",
               "rfb (cond. capped):\nPA - mult",
               "rfb (cond. capped):\nPA - all"),
    labels = c("(a) zero-fishing", "(b) 2 over 3", 
               "(c) rfb: default*",
               "(d) rfb: MSY - mult", "(e) rfb: MSY - all",
               "(f) rfb: MSY-PA - mult*", 
               "(g) rfb: MSY-PA - all*",
               "(h) rfb (capped):\n     MSY-PA - mult",
               "(i) rfb (cond. capped):\n     MSY-PA - mult",
               "(j) rfb (cond. capped):\n     MSY-PA - all"))) %>%
  mutate(group_numeric = factor(group, labels = c(1, 2.5, 4, 5.5, 6.5, 8, 9,
                                                  10.5, 11.5, 12.5))) %>%
  mutate(group_numeric = as.numeric(as.character(group_numeric)))
stats_plot_dev <- stats_plot %>%
  filter(name %in% c("SSB", "Catch")) %>%
  mutate(value_sign = ifelse(name == "SSB",
                             -abs(SSB_rel - 1)/2,
                             -abs(SSB_rel - 1) - abs(Catch_rel - 1)/2)) %>%
  mutate(sign = ifelse(name == "SSB",
                       ifelse(SSB_rel > 1, "+", "-"),
                       ifelse(Catch_rel > 1, "+", "-"))) %>%
  filter(abs(value) >= 0.02)

saveRDS(stats_plot, "output/plots/PA/data_all_comparison.rds")
stats_plot <- readRDS("output/plots/PA/data_all_comparison.rds")
saveRDS(stats_plot_dev, "output/plots/PA/data_all_comparison_dev.rds")
stats_plot_dev <- readRDS("output/plots/PA/data_all_comparison_dev.rds")


plot_comparison <- function(stock, data, data_dev, ylim) {
  data %>% 
    filter(stock %in% !!stock) %>%
    mutate(value = -abs(value)) %>%
    ggplot(aes(x = group_numeric, y = value, fill = name)) +
    geom_hline(yintercept = 0, linetype = "solid", size = 0.5, colour = "grey") +
    geom_col(position = "stack", width = 1, 
             colour = "black", size = 0.1) +
    ### for reversing order in legend (not in plot)
    scale_fill_discrete("fitness\nelements",
                        breaks = rev(levels(stats_plot$name))) +
    geom_text(data = data_dev %>%
                filter(stock %in% !!stock),
              aes(x = group_numeric, y = value_sign, label = sign),
              vjust = 0.5, colour = "grey20") +
    facet_grid(fhist ~ stock_k, scales = "free", space = "free_x", 
               labeller = label_parsed) +
    labs(y = "fitness", x = "") +
    scale_x_continuous(breaks = c(1, 2.5, 4, 5.5, 6.5, 8, 9, 10.5, 11.5, 12.5), 
                       labels = levels(data$group), 
                       minor_breaks = NULL) +
    ylim(ylim) +
    theme_bw(base_size = 8, base_family = "sans") +
    theme(panel.spacing.x = unit(0, units = "cm"),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5,
                                     lineheight = 0.7),
          #plot.margin = unit(x = c(1, 3, 3, 3), units = "pt"),
          legend.key.width = unit(1, "lines"),
          legend.key.height = unit(1, "lines"))
}

plot_comparison_six <- function(stocks, data, data_dev, ylim) {
  #browser()
  if (isTRUE(length(stocks) == 5)) {
    p1 <- NULL
    stocks <- c(NA, stocks)
  } else {
    p1 <- plot_comparison(stock = stocks[1], data = data, 
                          data_dev = data_dev, ylim = ylim) + 
      theme(strip.text.y = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            legend.position = "none",) +
      labs(y = "")
  }
  p2 <- plot_comparison(stock = stocks[2], data = data, 
                        data_dev = data_dev, ylim = ylim) + 
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank())
  p3 <- plot_comparison(stock = stocks[3], data = data, 
                        data_dev = data_dev, ylim = ylim) + 
    theme(strip.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none") +
    labs(y = "fitness")
  p4 <- plot_comparison(stock = stocks[4], data = data, 
                        data_dev = data_dev, ylim = ylim) + 
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none",
          axis.title.y = element_blank())
  p5 <- plot_comparison(stock = stocks[5], data = data, 
                        data_dev = data_dev, ylim = ylim) + 
    theme(strip.text.y = element_blank(),
          legend.position = "none") +
    labs(y = "")
  p6 <- plot_comparison(stock = stocks[6], data = data, 
                        data_dev = data_dev, ylim = ylim) + 
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),legend.position = "none",
          axis.title.y = element_blank())
  ### combine all plots
  p <- plot_grid(p1, p2 + theme(legend.position = "none"), get_legend(p2),
            p3, p4, NULL,
            p5, p6, NULL,
            ncol = 3, rel_heights = c(1, 1, 1.25), rel_widths = c(1, 1, 0.3))
  return(p)
}

### plot for manuscript
plot_comparison_six(stocks = c("ang3", "pol", "ple", "tur", "jnd", "her"),
                    data = stats_plot, data_dev = stats_plot_dev, 
                    ylim = c(-8.5, 0))
ggsave(filename = "output/plots/PA/all_comparison.png",
       width = 17, height = 18, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/PA/all_comparison.pdf",
       width = 17, height = 18, units = "cm", dpi = 600)

### plot all stocks for supplementary material
for (i in seq_along(split(stocks$stock, ceiling(seq_along(stocks$stock)/6)))) {
  #browser()
  plot_comparison_six(stocks = split(stocks$stock, 
                                     ceiling(seq_along(stocks$stock)/6))[[i]],
                      data = stats_plot, data_dev = stats_plot_dev, 
                      ylim = c(-8.5, 0))
  ggsave(filename = paste0("output/plots/PA/all_comparison", i , ".png"),
         width = 17, height = 18, units = "cm", dpi = 600, type = "cairo")
  ggsave(filename = paste0("output/plots/PA/all_comparison", i , ".pdf"),
         width = 17, height = 18, units = "cm", dpi = 600)
}
  

### ------------------------------------------------------------------------ ###
### table for Supplementary Material ####
### ------------------------------------------------------------------------ ###

all_GA <- readRDS("output/all_stocks_GA_optimised_stats.rds")
pol_PA <- readRDS("output/pol_PA_components_stats.rds")

### default fitness
fitness_default <- all_GA %>% 
  filter(optimised == "default" & capped == FALSE & scenario == "PA") %>%
  group_by(stock, fhist) %>%
  mutate(fitness_default = -sum(abs(SSB_rel - 1),
                        abs(Catch_rel - 1),
                        ICV,
                        penalty(x = risk_Blim, negative = FALSE, 
                                max = 5, inflection = 0.06, 
                                steepness = 0.5e+3))) %>%
  select(fhist, stock, fitness_default) %>%
  ungroup() %>%
  arrange(fhist)

res <- bind_rows(
  ### default
  all_GA %>% 
    filter(optimised == "default" & capped == FALSE & scenario == "PA" &
             stock == "pol" & fhist == "one-way") %>%
    mutate(stock = NA, fhist = NA) %>%
    mutate(scenario_name = "default", scenario_no = 1),
  ### Pollack parameter exploration
  pol_PA %>%
    filter(optimised != "default") %>%
    mutate(optimised = factor(optimised, 
                                 levels = c("mult", "cap", "mult_cap",
                                            "all", "all_cap"))) %>%
    arrange(as.numeric(optimised)) %>% 
    mutate(scenario_name = "parameter explorations", scenario_no = 2),
  ### multiplier
  all_GA %>% 
    filter(optimised == "mult" & capped == FALSE & scenario == "PA") %>%
    mutate(scenario_name = "rfb: PA - mult", scenario_no = 3),
  ### all parameters
  all_GA %>% 
    filter(optimised == "all_cap" & capped == FALSE & scenario == "PA") %>%
    mutate(scenario_name = "rfb: PA - all", scenario_no = 4),
  ### cap & multiplier
  all_GA %>% 
    filter(optimised == "mult" & capped == TRUE & scenario == "PA") %>%
    mutate(scenario_name = "rfb (capped): PA - mult", scenario_no = 5),
  ### conditional cap & multiplier
  all_GA %>% 
    filter(optimised == "mult" & capped == TRUE & scenario == "PA_capped") %>%
    mutate(scenario_name = "rfb (cond. capped): PA - mult", scenario_no = 6),
  ### all parameters
  all_GA %>% 
    filter(optimised == "all_cap" & capped == TRUE & scenario == "PA_capped") %>%
    mutate(scenario_name = "rfb (cond. capped): PA - all", scenario_no = 7)
) %>%
  full_join(fitness_default) %>%
  mutate(fitness_improvement = round((1 - fitness/fitness_default)*100)) %>%
  select(scenario_name, scenario_no, fhist, stock, iter, 
         lag_idx, range_idx_1, range_idx_2, exp_r, exp_f, 
         exp_b, interval, multiplier, upper_constraint, lower_constraint, 
         fitness_improvement) %>% 
  arrange(scenario_no, fhist) %>%
  select(-scenario_no)

write.csv(res, file = "output/PA_summary_table_parameters.csv", 
          row.names = FALSE)




