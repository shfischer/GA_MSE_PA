### ------------------------------------------------------------------------ ###
### analyse "multi-species" GA runs ####
### ------------------------------------------------------------------------ ###

library(doParallel)
library(doRNG)
library(GA)
library(ggplot2)
library(scales)
library(cowplot)
library(Cairo)
library(tidyr)
library(dplyr)
library(FLCore)
library(FLash)
library(FLBRP)
library(mseDL)

source("funs_GA.R")
source("funs.R")

trans_from <- function(from = 1) {
  trans <- function(x) x - from
  inv <- function(x) x + from
  trans_new("from", trans, inv, 
            domain = c(from, Inf))
}


stocks <- read.csv("input/stocks.csv", stringsAsFactors = FALSE)

### ------------------------------------------------------------------------ ###
### uncertainty cap ####
### ------------------------------------------------------------------------ ###
### GA runs for pollack

ga_solution <- function(object) {
  res <- tail(object@bestSol, 1)[[1]][1, ]
  names(res) <- object@names
  res[c(1:4, 8)] <- round(res[c(1:4, 8)])
  res[c(5:7)] <- round(res[c(5:7)], 1)
  res[c(9:11)] <- round(res[c(9:11)], 2)
  return(res)
}

res <- foreach(fhist = c("one-way", "random"), .combine = bind_rows) %:%
  foreach(obj = c("PA", "MSY", "MSYPA"), .combine = bind_rows) %:%
  foreach(params = c("default", "multiplier", "cap", "cap_multiplier", "full", 
                     "full_cap"), .combine = bind_rows) %:%
  foreach(stat_yrs = c("all", "last10"), .combine = bind_rows) %do% {#browser()
    #browser()
    #print(paste(fhist, obj, params, stat_yrs))
    ### load data
    par_file <- switch(params,
      "default" = "multiplier-upper_constraint-lower_constraint", 
      "multiplier" = "multiplier", 
      "cap" = "upper_constraint-lower_constraint", 
      "cap_multiplier" = "multiplier-upper_constraint-lower_constraint", 
      "full" = paste0("lag_idx-range_idx_1-range_idx_2-exp_r-exp_f-exp_b-",
                      "interval-multiplier"), 
      "full_cap" = paste0("lag_idx-range_idx_1-range_idx_2-exp_r-exp_f-exp_b-",
                      "interval-multiplier-upper_constraint-lower_constraint")
    )
    obj_file <- switch(obj,
      "PA" = "obj_ICES_PA2",
      "MSY" = "obj_SSB_C_risk_ICV",
      "MSYPA" = "obj_ICES_MSYPA"
    )
    path <- paste0("output/500_50/uncertainty_cap/", fhist, "/pol/",
                   par_file, "--", obj_file)
    path_runs <- paste0(path, "_runs", 
                        ifelse(stat_yrs == "all", "", paste0("_", stat_yrs)),
                        ".rds")
    path_res <- paste0(path, "_res", 
                        ifelse(stat_yrs == "all", "", paste0("_", stat_yrs)),
                        ".rds")
    ### use GA paper results for "full" GA
    # if (isTRUE(obj == "MSY" & params == "full")) {
    #   path <- paste0("output/500_50/ms/trial/", fhist, "/pol/",
    #                par_file, "--", obj_file)
    # }
    if (!file.exists(path_runs)) return(NULL)
    print("found something")
    ga_res <- readRDS(path_res)
    ga_runs <- readRDS(path_runs)
    
    ### optimised parameters
    if (isFALSE(params == "default")) {
      pars <- ga_solution(ga_res)
    } else {
      pars <- c(1, 2, 3, 1, 1, 1, 1, 2, 1, Inf, 0)
    }
    pars[which(is.nan(pars))] <- Inf
    tmp <- as.data.frame(t(pars))
    names(tmp) <- c("lag_idx", "range_idx_1", "range_idx_2", "range_catch",
                    "exp_r", "exp_f", "exp_b", "interval", "multiplier",
                    "upper_constraint", "lower_constraint")
    #if (is.nan(tmp$upper_constraint)) tmp$upper_constraint <- Inf
    tmp$obj <- obj
    tmp$fhist <- fhist
    tmp$stat_yrs_obj <- stat_yrs
    tmp$ga_obj <- params
    
    ### stats
    par_scn <- pars
    if (isTRUE(length(which(is.na(par_scn))) > 0)) {
      par_scn <- par_scn[-which(is.na(par_scn))]
    }
    par_scn <- paste0(par_scn, collapse = "_")
    stats_tmp <- ga_runs[[par_scn]]
    stats_tmp <- as.data.frame(lapply(as.data.frame(t(stats_tmp$stats)), unlist))
    
    ### combine pars and stats
    stats_tmp <- cbind(tmp, stats_tmp)
    
    ### if different stat_yrs period used, extract also default stats
    if (isFALSE(stat_yrs == "all")) {
      stats_tmp <- rbind(stats_tmp, stats_tmp)
      stats_tmp$stat_yrs <- c("all", stat_yrs)
      stats_names <- c("risk_Blim", "risk_Bmsy", "risk_halfBmsy", 
                       "risk_collapse", "SSB", "Fbar", "Catch", "SSB_rel", 
                       "Fbar_rel", "Catch_rel", "ICV")
      stats_tmp[2, stats_names] <- stats_tmp[2, paste0(stats_names, "_", stat_yrs)]
      stats_tmp[, paste0(stats_names, "_", stat_yrs)] <- NULL
    } else {
      stats_tmp$stat_yrs <- "all"
      ### remove redundant stats
      stats_tmp[, grep(x = names(stats_tmp), pattern = "_last10")] <- NULL
    }
    
    ### recreate fitness
    yr_suffix <- ""
    if (isTRUE(obj == "PA")) {
      stats_tmp$fitness <- sapply(seq(nrow(stats_tmp)), function(x) {
        sum(stats_tmp[x, "Catch_rel"]) -
          sum(penalty(x = stats_tmp[x, "risk_Blim"], negative = FALSE, max = 5,
                      inflection = 0.06, steepness = 0.5e+3))
      })
    } else if (isTRUE(obj == "MSY")) {
      stats_tmp$fitness <- sapply(seq(nrow(stats_tmp)), function(x) {
        -sum(abs(stats_tmp[x, "SSB_rel"] - 1),
                          abs(stats_tmp[x, "Catch_rel"] - 1),
                          stats_tmp[x, "ICV"], 
                          stats_tmp[x, "risk_Blim"])
      })
    } else if (isTRUE(obj == "MSYPA")) {
      stats_tmp$fitness <- sapply(seq(nrow(stats_tmp)), function(x) {
        -sum(abs(stats_tmp[x, "SSB_rel"] - 1),
                          abs(stats_tmp[x, "Catch_rel"] - 1),
                          stats_tmp[x, "ICV"], 
                          penalty(x = stats_tmp[x, "risk_Blim"], 
                                  negative = FALSE, max = 5, 
                                  inflection = 0.06, steepness = 0.5e+3))
      })
    }
    return(stats_tmp)
}

res %>%
  filter(obj == "PA")
saveRDS(res, file = "output/500_50/uncertainty_cap/results.rds")
res <- readRDS("output/500_50/uncertainty_cap/results.rds")


### format for plotting
stats_pol <- res %>%
  select(obj, fhist, stat_yrs_obj, stat_yrs, ga_obj, risk_Blim, SSB_rel, Fbar_rel, Catch_rel,
         ICV, fitness) %>%
  pivot_longer(c(SSB_rel, Fbar_rel, Catch_rel, risk_Blim, ICV, fitness), 
               names_to = "key", values_to = "value") %>%
  mutate(stat = factor(key, levels = c("SSB_rel", "Fbar_rel", "Catch_rel", "risk_Blim", 
                                       "ICV", "fitness"), 
                       labels = c("SSB/B[MSY]", "F/F[MSY]", "Catch/MSY", 
                                  "B[lim]~risk", "ICV", "fitness~value"))) %>%
  mutate(scenario = factor(ga_obj,
                           levels = c("default", "multiplier", "cap", 
                                      "cap_multiplier", "full", "full_cap"),
                           labels = c("default", "GA multiplier", "GA cap", 
                                      "GA cap+\nmultiplier", 
                                      "GA all w/o cap", "GA all")))
stats_targets <- data.frame(stat = c("SSB/B[MSY]", "F/F[MSY]", "Catch/MSY", 
                                     "B[lim]~risk", "ICV", "fitness~value"),
                            target = c(1, 1, 1, 0, 0, NA))


plot_stats <- function(SSB_min = 0, SSB_max = NA,
                       F_min = 0, F_max = NA,
                       C_min = 0, C_max = NA,
                       risk_min = 0, risk_max = NA,
                       ICV_min = 0, ICV_max = NA,
                       fitness_min = NA, fitness_max = NA,
                       obj = "MSY",
                       stat_yrs_obj = "all",
                       stat_yrs = "all",
                       data,
                       risk_line = FALSE
) {
  
  data <- data[data$obj == obj & data$stat_yrs_obj %in% stat_yrs_obj &
                 data$stat_yrs %in% stat_yrs, ]
  
  p_theme <- theme_bw(base_size = 8, base_family = "sans") +
    theme(panel.spacing.x = unit(0, units = "cm"),
          strip.placement.y = "outside",
          strip.background.y = element_blank(),
          strip.text.y = element_text(size = 8),
          plot.margin = unit(x = c(1, 3, 0, 3), units = "pt"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank())
  
  p_pol_stats_SSB <- data %>% 
    filter(stat %in% c("SSB/B[MSY]")) %>%
    ggplot(aes(x = scenario, y = value, fill = scenario,
               colour = scenario)) +
    geom_hline(yintercept = 1, linetype = "solid", size = 0.5, colour = "grey") +
    geom_col(position = position_dodge2(preserve = "single"), width = 0.8, 
             show.legend = FALSE, colour = "black", size = 0.1) +
    facet_grid(stat ~ fhist, scales = "free", space = "free_x", switch = "y",
               labeller = "label_parsed") +
    labs(y = "", x = "fitness function") +
    p_theme +
    scale_y_continuous(trans = trans_from(), limits = c(SSB_min, SSB_max))
  
  p_pol_stats_SSB <- data %>% 
    filter(stat %in% c("SSB/B[MSY]")) %>%
    ggplot(aes(x = scenario, y = value, fill = scenario,
               colour = scenario)) +
    geom_hline(yintercept = 1, linetype = "solid", size = 0.5, colour = "grey") +
    geom_col(position = position_dodge2(preserve = "single"), width = 0.8, 
             show.legend = FALSE, colour = "black", size = 0.1) +
    facet_grid(stat ~ fhist, scales = "free", space = "free_x", switch = "y",
               labeller = "label_parsed") +
    labs(y = "", x = "fitness function") +
    p_theme +
    scale_y_continuous(trans = trans_from(), limits = c(SSB_min, SSB_max))
  
  p_pol_stats_F <- data %>% 
    filter(stat %in% c("F/F[MSY]")) %>%
    ggplot(aes(x = scenario, y = value, fill = scenario,
               colour = scenario)) +
    geom_hline(yintercept = 1, linetype = "solid", size = 0.5, colour = "grey") +
    geom_col(position = "dodge", show.legend = FALSE, width = 0.8, 
             colour = "black", size = 0.1) +
    facet_grid(stat ~ fhist, scales = "free", space = "free_x", switch = "y",
               labeller = "label_parsed") +
    labs(y = "", x = "fitness function") +
    p_theme +
    theme(strip.text.x = element_blank(),
          plot.margin = unit(x = c(0, 3, 0, 3), units = "pt")) + 
    scale_y_continuous(trans = trans_from(), limits = c(F_min, F_max))
  p_pol_stats_C <- data %>% 
    filter(stat %in% c("Catch/MSY")) %>%
    ggplot(aes(x = scenario, y = value, fill = scenario,
               colour = scenario)) +
    geom_hline(yintercept = 1, linetype = "solid", size = 0.5, colour = "grey") +
    geom_col(position = "dodge", show.legend = FALSE, width = 0.8, 
             colour = "black", size = 0.1) +
    facet_grid(stat ~ fhist, scales = "free", space = "free_x", switch = "y",
               labeller = "label_parsed") +
    labs(y = "", x = "fitness function") +
    p_theme +
    theme(strip.text.x = element_blank(),
          plot.margin = unit(x = c(0, 3, 0, 3), units = "pt")) + 
    scale_y_continuous(trans = trans_from(), limits = c(C_min, C_max))
  p_pol_stats_risk <- data %>% 
    filter(stat %in% c("B[lim]~risk")) %>%
    ggplot(aes(x = scenario, y = value, fill = scenario,
               colour = scenario)) +
    geom_hline(yintercept = ifelse(isTRUE(risk_line), 0.05, 0), 
               linetype = "solid", size = 0.5, 
               colour = ifelse(isTRUE(risk_line), "red", "grey")) +
    geom_col(position = "dodge", show.legend = FALSE, width = 0.8, 
             colour = "black", size = 0.1) +
    facet_grid(stat ~ fhist, scales = "free", space = "free_x", switch = "y",
               labeller = "label_parsed") +
    labs(y = "", x = "fitness function") +
    p_theme +
    theme(strip.text.x = element_blank(),
          plot.margin = unit(x = c(0, 3, 0, 3), units = "pt")) + 
    scale_y_continuous(trans = trans_from(0), limits = c(risk_min, risk_max))
  p_pol_stats_ICV <- data %>% 
    filter(stat %in% c("ICV")) %>%
    ggplot(aes(x = scenario, y = value, fill = scenario,
               colour = scenario)) +
    geom_hline(yintercept = 0, linetype = "solid", size = 0.5, colour = "grey") +
    geom_col(position = "dodge", show.legend = FALSE, width = 0.8, 
             colour = "black", size = 0.1) +
    facet_grid(stat ~ fhist, scales = "free", space = "free_x", switch = "y",
               labeller = "label_parsed") +
    labs(y = "", x = "fitness function") +
    p_theme +
    theme(strip.text.x = element_blank(),
          plot.margin = unit(x = c(0, 3, 0, 3), units = "pt")) + 
    scale_y_continuous(trans = trans_from(0), limits = c(ICV_min, ICV_max))
  p_pol_stats_fitness <- data %>% 
    filter(stat %in% c("fitness~value")) %>%
    ggplot(aes(x = scenario, y = value, fill = scenario,
               colour = scenario)) +
    geom_hline(yintercept = 0, linetype = "solid", size = 0.5, colour = "grey") +
    geom_col(position = "dodge", show.legend = FALSE, width = 0.8, 
             colour = "black", size = 0.1) +
    facet_grid(stat ~ fhist, scales = "free", space = "free_x", switch = "y",
               labeller = "label_parsed") +
    labs(y = "", x = "") +
    theme_bw(base_size = 8, base_family = "sans") +
    theme(panel.spacing.x = unit(0, units = "cm"),
          strip.text.x = element_blank(),
          strip.placement.y = "outside",
          strip.background.y = element_blank(),
          strip.text.y = element_text(size = 8),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          plot.margin = unit(x = c(0, 3, 3, 3), units = "pt")) +
    scale_y_continuous(trans = trans_from(0), 
                       limits = c(fitness_min, fitness_max)#,
                       #breaks = c(0, -0.5, -1, -1.5), 
                       #minor_breaks = c(-0.25, -0.75, -1.25)
                       )
  p_pol_stats_comb <- plot_grid(p_pol_stats_SSB, p_pol_stats_F, p_pol_stats_C,
                                p_pol_stats_risk, p_pol_stats_ICV,
                                p_pol_stats_fitness,
                                ncol = 1, align = "v",
                                rel_heights = c(1.25, 1, 1, 1, 1, 2.1))
  return(p_pol_stats_comb)
}


### plots for MSY fitness function
plot_stats(obj = "MSY", stat_yrs_obj = "all", data = stats_pol,
           SSB_max = 1.5, F_max = 1.5, C_max = 1.5, risk_max = 1.5, 
           ICV_max = 1.5, fitness_min = -1.5)

ggsave(filename = "output/plots/PA/pol_GA_params_MSY.png",
       width = 17, height = 13, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/PA/pol_GA_params_MSY.pdf",
      width = 17, height = 13, units = "cm", dpi = 600)

### PA fitness function
plot_stats(obj = "PA", stat_yrs_obj = "all", data = stats_pol,
           SSB_max = NA, F_max = 1, C_max = 1, risk_max = 1, 
           ICV_max = 1, fitness_min = NA, fitness_max = 0.5, 
           risk_line = TRUE)

ggsave(filename = "output/plots/PA/pol_GA_params_PA.png",
       width = 17, height = 13, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/PA/pol_GA_params_PA.pdf",
      width = 17, height = 13, units = "cm", dpi = 600)



### MSY & PA fitness function
plot_stats(obj = "MSYPA", stat_yrs_obj = "all", data = stats_pol,
           SSB_max = NA, F_max = 1, C_max = 1, risk_max = 1, 
           ICV_max = 1, fitness_min = NA, fitness_max = NA, 
           risk_line = TRUE)

ggsave(filename = "output/plots/PA/pol_GA_params_MSYPA.png",
       width = 17, height = 13, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/PA/pol_GA_params_MSYPA.pdf",
      width = 17, height = 13, units = "cm", dpi = 600)

### MSY & PA fitness function & last 10 years for stats
plot_stats_2(obj = "MSYPA", 
             stat_yrs = c("all", "last10"), 
             stat_yrs_obj = c("all", "last10"),
             data = stats_pol,
             SSB_max = NA, F_max = 1.2, C_max = 1.2, risk_max = 1, 
             ICV_max = 1, fitness_min = NA, fitness_max = NA, 
             risk_line = TRUE)

ggsave(filename = "output/plots/PA/pol_GA_params_MSYPA_last10.png",
       width = 17, height = 13, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/PA/pol_GA_params_MSYPA_last10.pdf",
      width = 17, height = 13, units = "cm", dpi = 600)


### plot function, now for plotting several options
plot_stats_2 <- function(SSB_min = 0, SSB_max = NA,
                       F_min = 0, F_max = NA,
                       C_min = 0, C_max = NA,
                       risk_min = 0, risk_max = NA,
                       ICV_min = 0, ICV_max = NA,
                       fitness_min = NA, fitness_max = NA,
                       obj = "MSY",
                       stat_yrs_obj = "all",
                       stat_yrs = "all",
                       data,
                       risk_line = FALSE
) {
  
  data <- data[data$obj == obj & data$stat_yrs_obj %in% stat_yrs_obj &
                 data$stat_yrs %in% stat_yrs, ]
  
  data <- data %>%
    mutate(scenario2 = paste0("GA yrs: ", stat_yrs_obj, "\n",
                             "stat yrs: ", stat_yrs))
  
  p_theme <- theme_bw(base_size = 8, base_family = "sans") +
    theme(panel.spacing.x = unit(0, units = "cm"),
          strip.placement.y = "outside",
          strip.background.y = element_blank(),
          strip.text.y = element_text(size = 8),
          plot.margin = unit(x = c(1, 3, 0, 3), units = "pt"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank())
  
  p_pol_stats_SSB <- data %>% 
    filter(stat %in% c("SSB/B[MSY]")) %>%
    ggplot(aes(x = scenario, y = value, fill = scenario2,
               colour = scenario2)) +
    geom_hline(yintercept = 1, linetype = "solid", size = 0.5, colour = "grey") +
    geom_col(position = position_dodge2(preserve = "single"), width = 0.8, 
             show.legend = TRUE, colour = "black", size = 0.1) +
    scale_fill_discrete("") +
    facet_grid(stat ~ fhist, scales = "free", space = "free_x", switch = "y",
               labeller = "label_parsed") +
    labs(y = "", x = "fitness function") +
    p_theme +
    theme(legend.key.height = unit(1, "lines"),
          legend.key.width = unit(0.5, "lines")) +
    scale_y_continuous(trans = trans_from(), limits = c(SSB_min, SSB_max))
  
  p_pol_stats_F <- data %>% 
    filter(stat %in% c("F/F[MSY]")) %>%
    ggplot(aes(x = scenario, y = value, fill = scenario2,
               colour = scenario2)) +
    geom_hline(yintercept = 1, linetype = "solid", size = 0.5, colour = "grey") +
    geom_col(position = "dodge", show.legend = FALSE, width = 0.8, 
             colour = "black", size = 0.1) +
    facet_grid(stat ~ fhist, scales = "free", space = "free_x", switch = "y",
               labeller = "label_parsed") +
    labs(y = "", x = "fitness function") +
    p_theme +
    theme(strip.text.x = element_blank(),
          plot.margin = unit(x = c(0, 3, 0, 3), units = "pt")) + 
    scale_y_continuous(trans = trans_from(), limits = c(F_min, F_max))
  p_pol_stats_C <- data %>% 
    filter(stat %in% c("Catch/MSY")) %>%
    ggplot(aes(x = scenario, y = value, fill = scenario2,
               colour = scenario2)) +
    geom_hline(yintercept = 1, linetype = "solid", size = 0.5, colour = "grey") +
    geom_col(position = "dodge", show.legend = FALSE, width = 0.8, 
             colour = "black", size = 0.1) +
    facet_grid(stat ~ fhist, scales = "free", space = "free_x", switch = "y",
               labeller = "label_parsed") +
    labs(y = "", x = "fitness function") +
    p_theme +
    theme(strip.text.x = element_blank(),
          plot.margin = unit(x = c(0, 3, 0, 3), units = "pt")) + 
    scale_y_continuous(trans = trans_from(), limits = c(C_min, C_max))
  p_pol_stats_risk <- data %>% 
    filter(stat %in% c("B[lim]~risk")) %>%
    ggplot(aes(x = scenario, y = value, fill = scenario2,
               colour = scenario2)) +
    geom_hline(yintercept = ifelse(isTRUE(risk_line), 0.05, 0), 
               linetype = "solid", size = 0.5, 
               colour = ifelse(isTRUE(risk_line), "red", "grey")) +
    geom_col(position = "dodge", show.legend = FALSE, width = 0.8, 
             colour = "black", size = 0.1) +
    facet_grid(stat ~ fhist, scales = "free", space = "free_x", switch = "y",
               labeller = "label_parsed") +
    labs(y = "", x = "fitness function") +
    p_theme +
    theme(strip.text.x = element_blank(),
          plot.margin = unit(x = c(0, 3, 0, 3), units = "pt")) + 
    scale_y_continuous(trans = trans_from(0), limits = c(risk_min, risk_max))
  p_pol_stats_ICV <- data %>% 
    filter(stat %in% c("ICV")) %>%
    ggplot(aes(x = scenario, y = value, fill = scenario2,
               colour = scenario2)) +
    geom_hline(yintercept = 0, linetype = "solid", size = 0.5, colour = "grey") +
    geom_col(position = "dodge", show.legend = FALSE, width = 0.8, 
             colour = "black", size = 0.1) +
    facet_grid(stat ~ fhist, scales = "free", space = "free_x", switch = "y",
               labeller = "label_parsed") +
    labs(y = "", x = "fitness function") +
    p_theme +
    theme(strip.text.x = element_blank(),
          plot.margin = unit(x = c(0, 3, 0, 3), units = "pt")) + 
    scale_y_continuous(trans = trans_from(0), limits = c(ICV_min, ICV_max))
  p_pol_stats_fitness <- data %>% 
    filter(stat %in% c("fitness~value")) %>%
    ggplot(aes(x = scenario, y = value, fill = scenario2,
               colour = scenario2)) +
    geom_hline(yintercept = 0, linetype = "solid", size = 0.5, colour = "grey") +
    geom_col(position = "dodge", show.legend = FALSE, width = 0.8, 
             colour = "black", size = 0.1) +
    facet_grid(stat ~ fhist, scales = "free", space = "free_x", switch = "y",
               labeller = "label_parsed") +
    labs(y = "", x = "") +
    theme_bw(base_size = 8, base_family = "sans") +
    theme(panel.spacing.x = unit(0, units = "cm"),
          strip.text.x = element_blank(),
          strip.placement.y = "outside",
          strip.background.y = element_blank(),
          strip.text.y = element_text(size = 8),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          plot.margin = unit(x = c(0, 3, 3, 3), units = "pt")) +
    scale_y_continuous(trans = trans_from(0), 
                       limits = c(fitness_min, fitness_max)#,
                       #breaks = c(0, -0.5, -1, -1.5), 
                       #minor_breaks = c(-0.25, -0.75, -1.25)
                       )
  p_pol_stats_comb <- plot_grid(
    plot_grid(p_pol_stats_SSB + theme(legend.position = "none"), 
              p_pol_stats_F, p_pol_stats_C, p_pol_stats_risk, p_pol_stats_ICV,
              p_pol_stats_fitness,
              ncol = 1, align = "v", rel_heights = c(1.25, 1, 1, 1, 1, 2.1)),
    get_legend(p_pol_stats_SSB), rel_widths = c(1, 0.2), ncol = 2
  )
  return(p_pol_stats_comb)
}


### ------------------------------------------------------------------------ ###
### fitness penalty visualisation ####
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
  coord_cartesian(xlim = c(0, 0.3)) +
  labs(x = expression(B[lim]~risk),
       y = "fitness penalty")
p
ggsave(filename = "output/plots/PA/Blim_penalty_curve.png",
       width = 8.5, height = 6, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/PA/Blim_penalty_curve.pdf",
      width = 8.5, height = 6, units = "cm", dpi = 600)



### ------------------------------------------------------------------------ ###
### multiplier all stocks ####
### ------------------------------------------------------------------------ ###

stocks_subset <- stocks$stock[1:29]

mult_all <- foreach(fhist = c("one-way", "random"), .combine = bind_rows) %:%
  foreach(stock = stocks_subset, .combine = bind_rows) %do% {
    #browser()
    runs <- readRDS(paste0("output/500_50/uncertainty_cap/", fhist, "/", stock, "/",
                   "multiplier--obj_ICES_MSYPA_runs_last10.rds"))
    
    runs <- lapply(runs, function(x) {
      s <- as.data.frame(rbind(x$stats[1:11], x$stats[12:22]))
      names(s) <- rownames(x$stats)[1:11]
      tmp <- cbind(t(x$pars), s)
      tmp <- as.data.frame(lapply(as.data.frame(tmp), unlist))
      tmp$stat_yrs <- c("all", "last10")
      return(tmp)
    })
    runs <- do.call(rbind, runs)
    row.names(runs) <- NULL
    runs$fhist <- fhist
    runs$stock <- stock
    
    return(runs)
}
mult_all <- mult_all %>% 
    left_join(stocks[, c("stock", "k")]) %>%
    mutate(stock_k = paste0(stock, "~(italic(k)==", k, ")")) %>%
    mutate(stock_k = factor(stock_k, levels = unique(stock_k)))



mult_all %>% 
  filter(stat_yrs == "last10" & fhist == "random") %>%
  ggplot(aes(x = multiplier, y = risk_Blim)) +
  geom_line() +
  geom_hline(yintercept = 0.05, colour = "red") +
  facet_wrap(~ stock_k, labeller = "label_parsed", scales = "free") +
  coord_cartesian(xlim = c(0, 1))



### ------------------------------------------------------------------------ ###
### MSYPA fitness function - all stocks ####
### ------------------------------------------------------------------------ ###

fhist <- "one-way"
stocks_subset <- stocks$stock[1:29]
res <- foreach(stock = stocks_subset, .combine = bind_rows) %:%
  foreach(obj = c("MSYPA"), .combine = bind_rows) %:%
  foreach(params = c("default", "multiplier", "full_cap"), 
          .combine = bind_rows) %:%
  foreach(stat_yrs = c("all"), .combine = bind_rows) %do% {#browser()
    #browser()
    ### load data
    par_file <- switch(params,
      "default" = "multiplier", 
      "multiplier" = "multiplier", 
      "cap" = "upper_constraint-lower_constraint", 
      "cap_multiplier" = "multiplier-upper_constraint-lower_constraint", 
      "full" = paste0("lag_idx-range_idx_1-range_idx_2-exp_r-exp_f-exp_b-",
                      "interval-multiplier"), 
      "full_cap" = paste0("lag_idx-range_idx_1-range_idx_2-exp_r-exp_f-exp_b-",
                      "interval-multiplier-upper_constraint-lower_constraint")
    )
    obj_file <- switch(obj,
      "PA" = "obj_ICES_PA2",
      "MSY" = "obj_SSB_C_risk_ICV",
      "MSYPA" = "obj_ICES_MSYPA"
    )
    path <- paste0("output/500_50/uncertainty_cap/", fhist, "/", stock, "/",
                   par_file, "--", obj_file)
    path_runs <- paste0(path, "_runs", 
                        ifelse(stat_yrs == "all", "", paste0("_", stat_yrs)),
                        ".rds")
    path_res <- paste0(path, "_res", 
                        ifelse(stat_yrs == "all", "", paste0("_", stat_yrs)),
                        ".rds")
    ### use GA paper results for "full" GA
    # if (isTRUE(obj == "MSY" & params == "full")) {
    #   path <- paste0("output/500_50/ms/trial/", fhist, "/pol/",
    #                par_file, "--", obj_file)
    # }
    if (!file.exists(path_runs)) return(NULL)
    print("found something")
    ga_res <- readRDS(path_res)
    ga_runs <- readRDS(path_runs)
    
    ### optimised parameters
    if (isFALSE(params == "default")) {
      pars <- ga_solution(ga_res)
    } else {
      pars <- c(1, 2, 3, 1, 1, 1, 1, 2, 1, Inf, 0)
    }
    pars[which(is.nan(pars))] <- Inf
    tmp <- as.data.frame(t(pars))
    names(tmp) <- c("lag_idx", "range_idx_1", "range_idx_2", "range_catch",
                    "exp_r", "exp_f", "exp_b", "interval", "multiplier",
                    "upper_constraint", "lower_constraint")
    #if (is.nan(tmp$upper_constraint)) tmp$upper_constraint <- Inf
    tmp$obj <- obj
    tmp$fhist <- fhist
    tmp$stat_yrs_obj <- stat_yrs
    tmp$ga_obj <- params
    tmp$stock <- stock
    
    ### stats
    par_scn <- pars
    if (isTRUE(length(which(is.na(par_scn))) > 0)) {
      par_scn <- par_scn[-which(is.na(par_scn))]
    }
    par_scn <- paste0(par_scn, collapse = "_")
    stats_tmp <- ga_runs[[par_scn]]
    stats_tmp <- as.data.frame(lapply(as.data.frame(t(stats_tmp$stats)), unlist))
    
    ### combine pars and stats
    stats_tmp <- cbind(tmp, stats_tmp)
    
    ### if different stat_yrs period used, extract also default stats
    if (isFALSE(stat_yrs == "all")) {
      stats_tmp <- rbind(stats_tmp, stats_tmp)
      stats_tmp$stat_yrs <- c("all", stat_yrs)
      stats_names <- c("risk_Blim", "risk_Bmsy", "risk_halfBmsy", 
                       "risk_collapse", "SSB", "Fbar", "Catch", "SSB_rel", 
                       "Fbar_rel", "Catch_rel", "ICV")
      stats_tmp[2, stats_names] <- stats_tmp[2, paste0(stats_names, "_", stat_yrs)]
      stats_tmp[, paste0(stats_names, "_", stat_yrs)] <- NULL
    } else {
      stats_tmp$stat_yrs <- "all"
      ### remove redundant stats
      stats_tmp[, grep(x = names(stats_tmp), pattern = "_last10")] <- NULL
    }
    
    ### recreate fitness
    yr_suffix <- ""
    if (isTRUE(obj == "PA")) {
      stats_tmp$fitness <- sapply(seq(nrow(stats_tmp)), function(x) {
        sum(stats_tmp[x, "Catch_rel"]) -
          sum(penalty(x = stats_tmp[x, "risk_Blim"], negative = FALSE, max = 5,
                      inflection = 0.06, steepness = 0.5e+3))
      })
    } else if (isTRUE(obj == "MSY")) {
      stats_tmp$fitness <- sapply(seq(nrow(stats_tmp)), function(x) {
        -sum(abs(stats_tmp[x, "SSB_rel"] - 1),
                          abs(stats_tmp[x, "Catch_rel"] - 1),
                          stats_tmp[x, "ICV"], 
                          stats_tmp[x, "risk_Blim"])
      })
    } else if (isTRUE(obj == "MSYPA")) {
      stats_tmp$fitness <- sapply(seq(nrow(stats_tmp)), function(x) {
        -sum(abs(stats_tmp[x, "SSB_rel"] - 1),
                          abs(stats_tmp[x, "Catch_rel"] - 1),
                          stats_tmp[x, "ICV"], 
                          penalty(x = stats_tmp[x, "risk_Blim"], 
                                  negative = FALSE, max = 5, 
                                  inflection = 0.06, steepness = 0.5e+3))
      })
    }
    return(stats_tmp)
}
res <- res %>%
  mutate(ga_obj = factor(ga_obj, levels = c("default", "multiplier", "full_cap"),
                         labels = c("default", "GA multiplier", 
                                    "GA all"))) %>% 
  left_join(stocks[, c("stock", "k")]) %>%
  mutate(stock_k = paste0(stock, "~(italic(k)==", k, ")")) %>%
  mutate(stock_k = factor(stock_k, levels = unique(stock_k))) %>%
  mutate(stock = factor(stock, levels = stocks$stock))
  


saveRDS(res, file = "output/500_50/uncertainty_cap/all_stocks_mult_full.rds")
res <- readRDS("output/500_50/uncertainty_cap/all_stocks_mult_full.rds")

### plot
res_plot <- res %>%
  select(obj, fhist, stat_yrs_obj, stat_yrs, ga_obj, risk_Blim, SSB_rel, 
         Fbar_rel, Catch_rel, ICV, fitness, stock) %>%
  pivot_longer(c(SSB_rel, Fbar_rel, Catch_rel, risk_Blim, ICV, fitness), 
               names_to = "key", values_to = "value") %>%
  mutate(stat = factor(key, levels = c("SSB_rel", "Fbar_rel", "Catch_rel",
                                       "risk_Blim", "ICV", "fitness"), 
                       labels = c("SSB/B[MSY]", "F/F[MSY]", "Catch/MSY", 
                                  "B[lim]~risk", "ICV", "fitness~value")))
stats_targets <- data.frame(stat = c("SSB/B[MSY]", "F/F[MSY]", "Catch/MSY", 
                                     "B[lim]~risk", "ICV", "fitness~value"),
                            target = c(1, 1, 1, 0, 0, NA))

p_theme <- theme_bw(base_size = 8, base_family = "sans") +
  theme(panel.spacing.x = unit(0, units = "cm"),
        strip.placement.y = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        plot.margin = unit(x = c(1, 3, 0, 3), units = "pt"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())

p_stats_SSB <- res_plot %>% 
  filter(stat %in% c("SSB/B[MSY]")) %>%
  ggplot(aes(x = stock, y = value, fill = ga_obj,
             colour = ga_obj)) +
  geom_hline(yintercept = 1, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = position_dodge2(preserve = "single"), width = 0.8, 
           show.legend = TRUE, colour = "black", size = 0.1) +
  scale_fill_discrete("") +
  facet_grid(stat ~ fhist, scales = "free", space = "free_x", switch = "y",
             labeller = "label_parsed") +
  labs(y = "", x = "fitness function") +
  p_theme +
  theme(legend.key.height = unit(1, "lines"),
        legend.key.width = unit(0.5, "lines")) +
  scale_y_continuous(trans = trans_from(), limits = c(0, NA))

p_stats_F <- res_plot %>% 
  filter(stat %in% c("F/F[MSY]")) %>%
  ggplot(aes(x = stock, y = value, fill = ga_obj,
             colour = ga_obj)) +
  geom_hline(yintercept = 1, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = "dodge", show.legend = FALSE, width = 0.8, 
           colour = "black", size = 0.1) +
  facet_grid(stat ~ fhist, scales = "free", space = "free_x", switch = "y",
             labeller = "label_parsed") +
  labs(y = "", x = "fitness function") +
  p_theme +
  theme(strip.text.x = element_blank(),
        plot.margin = unit(x = c(0, 3, 0, 3), units = "pt")) + 
  scale_y_continuous(trans = trans_from(), limits = c(0, NA))
p_stats_C <- res_plot %>% 
  filter(stat %in% c("Catch/MSY")) %>%
  ggplot(aes(x = stock, y = value, fill = ga_obj,
             colour = ga_obj)) +
  geom_hline(yintercept = 1, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = "dodge", show.legend = FALSE, width = 0.8, 
           colour = "black", size = 0.1) +
  facet_grid(stat ~ fhist, scales = "free", space = "free_x", switch = "y",
             labeller = "label_parsed") +
  labs(y = "", x = "fitness function") +
  p_theme +
  theme(strip.text.x = element_blank(),
        plot.margin = unit(x = c(0, 3, 0, 3), units = "pt")) + 
  scale_y_continuous(trans = trans_from(), limits = c(0, NA))
p_stats_risk <- res_plot %>% 
  filter(stat %in% c("B[lim]~risk")) %>%
  ggplot(aes(x = stock, y = value, fill = ga_obj,
             colour = ga_obj)) +
  geom_hline(yintercept = ifelse(isTRUE(TRUE), 0.05, 0), 
             linetype = "solid", size = 0.5, 
             colour = ifelse(isTRUE(TRUE), "red", "grey")) +
  geom_col(position = "dodge", show.legend = FALSE, width = 0.8, 
           colour = "black", size = 0.1) +
  facet_grid(stat ~ fhist, scales = "free", space = "free_x", switch = "y",
             labeller = "label_parsed") +
  labs(y = "", x = "fitness function") +
  p_theme +
  theme(strip.text.x = element_blank(),
        plot.margin = unit(x = c(0, 3, 0, 3), units = "pt")) + 
  scale_y_continuous(trans = trans_from(0), limits = c(0, 1))
p_stats_ICV <- res_plot %>% 
  filter(stat %in% c("ICV")) %>%
  ggplot(aes(x = stock, y = value, fill = ga_obj,
             colour = ga_obj)) +
  geom_hline(yintercept = 0, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = "dodge", show.legend = FALSE, width = 0.8, 
           colour = "black", size = 0.1) +
  facet_grid(stat ~ fhist, scales = "free", space = "free_x", switch = "y",
             labeller = "label_parsed") +
  labs(y = "", x = "fitness function") +
  p_theme +
  theme(strip.text.x = element_blank(),
        plot.margin = unit(x = c(0, 3, 0, 3), units = "pt")) + 
  scale_y_continuous(trans = trans_from(0), limits = c(0, 1))
p_stats_fitness <- res_plot %>% 
  filter(stat %in% c("fitness~value")) %>%
  ggplot(aes(x = stock, y = value, fill = ga_obj,
             colour = ga_obj)) +
  geom_hline(yintercept = 0, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = "dodge", show.legend = FALSE, width = 0.8, 
           colour = "black", size = 0.1) +
  facet_grid(stat ~ fhist, scales = "free", space = "free_x", switch = "y",
             labeller = "label_parsed") +
  labs(y = "", x = "") +
  theme_bw(base_size = 8, base_family = "sans") +
  theme(panel.spacing.x = unit(0, units = "cm"),
        strip.text.x = element_blank(),
        strip.placement.y = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        plot.margin = unit(x = c(0, 3, 3, 3), units = "pt")) +
  scale_y_continuous(trans = trans_from(0), 
                     limits = c(NA, NA)#,
                     #breaks = c(0, -0.5, -1, -1.5), 
                     #minor_breaks = c(-0.25, -0.75, -1.25)
                     )
p_stats_comb <- plot_grid(
  plot_grid(p_stats_SSB + theme(legend.position = "none"), 
            p_stats_F, p_stats_C, p_stats_risk, p_stats_ICV,
            p_stats_fitness,
            ncol = 1, align = "v", rel_heights = c(1.25, 1, 1, 1, 1, 1.5)),
  get_legend(p_stats_SSB), rel_widths = c(1, 0.2), ncol = 2
)

ggsave(filename = "output/plots/PA/all_GA_params_MSYPA.png",
       width = 17, height = 13, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/PA/all_GA_params_MSYPA.pdf",
      width = 17, height = 13, units = "cm", dpi = 600)
