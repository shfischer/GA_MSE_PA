### ------------------------------------------------------------------------ ###
### analyse MSY GA runs with PA fitness function ####
### ------------------------------------------------------------------------ ###

req_pckgs <- c("doParallel", "doRNG", "mse", "GA", "ggplot2", "cowplot", 
               "scales", "Cairo", "tidyr", "dplyr", "FLCore", "FLash", "FLBRP")
for (i in req_pckgs) library(package = i, character.only = TRUE)

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
### collate results - all stocks ####
### ------------------------------------------------------------------------ ###
n_yrs <- 50
n_iter <- 500

### get optimised parameterisation
all_PA <- foreach(stock = stocks$stock, .combine = rbind) %:%
  foreach(optimised = c("zero-fishing", "default", "mult", "all_cap"),
          .combine = rbind) %:%
  foreach(scenario = "PA", .combine = rbind) %:%
  foreach(stat_yrs = "more", .combine = rbind) %:%
  foreach(fhist = c("one-way", "random"), .combine = rbind) %do% {#browser()
    #browser()
    ### find files
    file_name <- switch(optimised,
      "zero-fishing" = "multiplier--obj_ICES_MSYPA",
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
    if (isTRUE(optimised %in% c("default", "zero-fishing"))) {
      trials <- trials[1, ]
      res_lst <- res_lst[1]
      res_par <- res_par[1]
      res_par[[1]][] <- c(1, 2, 3, 1, 1, 1, 1, 2, 1, Inf, 0)
      if (identical(optimised, "zero-fishing")) res_par[[1]][9] <- 0
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
saveRDS(all_PA, file = "output/all_stocks_PA_stats.rds")
write.csv(all_PA, file = "output/all_stocks_PA_stats.csv", row.names = FALSE)

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
### correct ICV when multiplier = 0
all_mult <- all_mult %>% 
  mutate(ICV = ifelse(multiplier == 0, 0, ICV)) %>%
  select(1:22, stock, fhist)

saveRDS(all_mult, file = "output/all_stocks_PA_multiplier_stats.rds")
write.csv(all_mult, file = "output/all_stocks_PA_multiplier_stats.csv", 
          row.names = FALSE)

### ------------------------------------------------------------------------ ###
### plot - all stocks multipliers ####
### ------------------------------------------------------------------------ ###

all_mult <- readRDS("output/all_stocks_PA_multiplier_stats.rds")

### format for plotting
all_mult <- all_mult %>%
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
        file = "output/plots/PA/data_all_stocks_multiplier_stats.rds")
all_mult_refs <- readRDS("output/plots/PA/data_all_stocks_multiplier_stats.rds")


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
             aes(shape = point), size = 0.8) +
  scale_shape_manual("risk limit", values = 2:4, 
                     guide = guide_legend(order = 2)) + 
  scale_colour_discrete("fishing\nhistory", guide = guide_legend(order = 1)) +
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
               aes(shape = point), size = 0.8) +
    scale_shape_manual("risk limit", values = 2:4, 
                       guide = guide_legend(order = 2)) + 
    scale_colour_discrete("fishing\nhistory", guide = guide_legend(order = 1)) +
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
### plot - pollack rfb-rule component explorations ####
### ------------------------------------------------------------------------ ###

pol <- readRDS("output/pol_PA_components_stats.rds")

pol$label <- as.character(pol$optimised)
pol$label[pol$label == "default"] <- "not\noptimised"
pol$label[pol$label == "mult"] <- "GA: multi-\nplier"
pol$label[pol$label == "cap"] <- "GA: uncer-\ntainty cap"
pol$label[pol$label == "mult_cap"] <- "GA: multi-\nplier and cap"
pol$label[pol$label == "all"] <- "GA: all with-\nout cap"
pol$label[pol$label == "all_cap"] <- "GA: all"
pol$label <- as.factor(pol$label)
pol$label <- factor(pol$label, 
                              levels = levels(pol$label)[c(6, 3, 5, 4, 2, 1)])


pol_plot <- pol %>%
  pivot_longer(c(SSB_rel, Fbar_rel, Catch_rel, risk_Blim, ICV, fitness), 
               names_to = "key", values_to = "value") %>%
  mutate(stat = factor(key, levels = c("SSB_rel", "Fbar_rel", "Catch_rel", "risk_Blim", 
                                       "ICV", "fitness"), 
                       labels = c("SSB/B[MSY]", "F/F[MSY]", "Catch/MSY", 
                                  "B[lim]~risk", "ICV", "fitness~value")))
stats_targets <- data.frame(stat = c("SSB/B[MSY]", "F/F[MSY]", "Catch/MSY", 
                                     "B[lim]~risk", "ICV", "fitness~value"),
                            target = c(1, 1, 1, 0, 0, NA))

saveRDS(pol_plot, file = "output/plots/PA/data_pol_components_stats.rds")
pol_plot <- readRDS("output/plots/PA/data_pol_components_stats.rds")
saveRDS(stats_targets, file = "output/plots/PA/data_pol_components_targets.rds")
stats_targets <- readRDS("output/plots/PA/data_pol_components_targets.rds")

### individual plots
y_max <- 3.25
p_pol_stats_SSB <- pol_plot %>% 
  filter(stat %in% c("SSB/B[MSY]")) %>%
  ggplot(aes(x = label, y = value, fill = label,
             colour = label)) +
  geom_hline(yintercept = 1, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = position_dodge2(preserve = "single"), width = 0.8, 
           show.legend = FALSE, colour = "black", size = 0.1) +
  facet_grid(stat ~ fhist, scales = "free", space = "free_x", switch = "y",
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
  scale_y_continuous(trans = trans_from(), limits = c(0, y_max))
p_pol_stats_F <- pol_plot %>% 
  filter(stat %in% c("F/F[MSY]")) %>%
  ggplot(aes(x = label, y = value, fill = label,
             colour = label)) +
  geom_hline(yintercept = 1, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = "dodge", show.legend = FALSE, width = 0.8, 
           colour = "black", size = 0.1) +
  facet_grid(stat ~ fhist, scales = "free", space = "free_x", switch = "y",
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
  scale_y_continuous(trans = trans_from(), limits = c(0, y_max))
p_pol_stats_C <- pol_plot %>% 
  filter(stat %in% c("Catch/MSY")) %>%
  ggplot(aes(x = label, y = value, fill = label,
             colour = label)) +
  geom_hline(yintercept = 1, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = "dodge", show.legend = FALSE, width = 0.8, 
           colour = "black", size = 0.1) +
  facet_grid(stat ~ fhist, scales = "free", space = "free_x", switch = "y",
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
  scale_y_continuous(trans = trans_from(), limits = c(0, y_max))
p_pol_stats_risk <- pol_plot %>% 
  filter(stat %in% c("B[lim]~risk")) %>%
  ggplot(aes(x = label, y = value, fill = label,
             colour = label)) +
  geom_hline(yintercept = 0, linetype = "solid", size = 0.5, colour = "grey") +
  geom_hline(yintercept = 0.05, linetype = "solid", size = 0.5, colour = "red") +
  geom_col(position = "dodge", show.legend = FALSE, width = 0.8, 
           colour = "black", size = 0.1) +
  facet_grid(stat ~ fhist, scales = "free", space = "free_x", switch = "y",
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
p_pol_stats_ICV <- pol_plot %>% 
  filter(stat %in% c("ICV")) %>%
  ggplot(aes(x = label, y = value, fill = label,
             colour = label)) +
  geom_hline(yintercept = 0, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = "dodge", show.legend = FALSE, width = 0.8, 
           colour = "black", size = 0.1) +
  facet_grid(stat ~ fhist, scales = "free", space = "free_x", switch = "y",
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
p_pol_stats_fitness <- pol_plot %>% 
  filter(stat %in% c("fitness~value")) %>%
  ggplot(aes(x = label, y = value, fill = label,
             colour = label)) +
  geom_hline(yintercept = 0, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = "dodge", show.legend = FALSE, width = 0.8, 
           colour = "black", size = 0.1) +
  facet_grid(stat ~ fhist, scales = "free", space = "free_x", switch = "y",
             labeller = "label_parsed") +
  labs(y = "", x = "\nrfb-rule components") +
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
  scale_y_continuous(trans = trans_from(0), limits = c(-NA, NA))#,
                     #breaks = c(0, -0.5, -1), 
                     #minor_breaks = c(-0.25, -0.75, -1.25))
p_pol_stats_comb <- plot_grid(p_pol_stats_SSB, 
                              #p_pol_stats_F, 
                              p_pol_stats_C,
                              p_pol_stats_risk, 
                              p_pol_stats_ICV,
                              p_pol_stats_fitness,
                              ncol = 1, align = "v",
                              rel_heights = c(1.25, 1, 1, 1, 2.15))
ggsave(filename = "output/plots/PA/pol_components_stats.png", 
       plot = p_pol_stats_comb,
       width = 8.5, height = 11, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/PA/pol_components_stats.pdf", 
       plot = p_pol_stats_comb,
       width = 8.5, height = 11, units = "cm", dpi = 600)

### ------------------------------------------------------------------------ ###
### plot - all stocks stats PA - default vs. optimised ####
### ------------------------------------------------------------------------ ###

### load data
all_PA <- readRDS("output/all_stocks_PA_stats.rds")
stocks_sorted <- stocks %>%
  select(stock, k) %>%
  arrange(k)

### format for plotting
stats_plot <- all_PA %>% 
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

### individual plots
p_stats_SSB <- stats_plot %>%
  filter(stat %in% c("SSB/B[MSY]")) %>%
  ggplot(aes(x = stock, y = value, fill = rule, colour = rule)) +
  geom_hline(yintercept = 1, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = position_dodge2(), width = 0.8, 
           show.legend = FALSE, size = 0.1) +
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
  scale_fill_discrete("") + scale_colour_discrete("") +
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


ggsave(filename = "output/plots/PA/all_stocks_stats.png", 
       plot = p_stats_combined,
       width = 17, height = 10, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/PA/all_stocks_stats.pdf", 
       plot = p_stats_combined,
       width = 17, height = 10, units = "cm", dpi = 600)


### ------------------------------------------------------------------------ ###
### plot - all stocks stats PA vs. 2 over 3 ####
### ------------------------------------------------------------------------ ###

### load data
all_PA <- readRDS("output/all_stocks_PA_stats.rds")
all_MSY <- readRDS("output/all_stocks_MSY_stats.rds")

stats_2over3 <- foreach(stock = stocks$stock, .combine = "rbind") %:%
  foreach(fhist = c("one-way", "random"), .combine = "rbind") %do% {
  
    stats_i <- readRDS(paste0("output/500_50/2over3/", fhist, "/", stock, 
                              "_stats.rds"))
    stats_i <- as.data.frame(lapply(as.data.frame(t(stats_i)), unlist))
    stats_i$fhist <- fhist
    stats_i$stock = stock
    stats_i$catch_rule = "2 over 3"
    stats_i$group = "2 over 3"
    return(stats_i)
}

### combine 
stats_plot <- bind_rows(
  all_MSY %>% 
    filter(optimised == TRUE) %>%
    select(fhist, stock, SSB_rel, Fbar_rel, Catch_rel, risk_Blim, ICV) %>%
    mutate(catch_rule = "rfb", group = "rfb: MSY"),
  all_PA %>% 
    filter(optimised == "all_cap") %>%
    select(fhist, stock, SSB_rel, Fbar_rel, Catch_rel, risk_Blim, ICV) %>%
    mutate(catch_rule = "rfb", group = "rfb: MSY-PA"),
  stats_2over3 %>%
    select(fhist, stock, SSB_rel, Fbar_rel, Catch_rel, risk_Blim, ICV) %>%
    mutate(catch_rule = "2 over 3", group = "2 over 3")
)
saveRDS(stats_plot, "output/all_stocks_rfb_opt_stats.rds")
stats_plot <- readRDS("output/all_stocks_rfb_opt_stats.rds")


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
  scale_fill_discrete("catch rule") + scale_colour_discrete("catch rule") +
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
            p_stats_risk, 
            p_stats_ICV + theme(legend.position = "none"),
            ncol = 1, align = "v",
            rel_heights = c(1.25, 1, 1, 1.5)),
  get_legend(p_stats_ICV), rel_widths = c(1, 0.15))


ggsave(filename = "output/plots/PA/all_stocks_2over3_stats.png", 
       plot = p_stats_combined,
       width = 17, height = 10, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/PA/all_stocks_2over3_stats.pdf", 
       plot = p_stats_combined,
       width = 17, height = 10, units = "cm", dpi = 600)


### ------------------------------------------------------------------------ ###
### table with rfb-rule parameters ####
### ------------------------------------------------------------------------ ###

all_PA %>%
  filter(optimised %in% c("mult", "all_cap")) %>%
  select("optimised", "fhist", "stock", "iter", "lag_idx", "range_idx_1", "range_idx_2",
         "exp_r", "exp_f", "exp_b", "interval", "multiplier", 
         "upper_constraint", "lower_constraint") %>%
  arrange(desc(optimised), fhist) %>%
  View()

all_PA %>% 
  filter(optimised == "default") %>%
  filter(risk_Blim <= 0.055)


