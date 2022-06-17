### ------------------------------------------------------------------------ ###
### analysis of constant harvest rate rule ####
### ------------------------------------------------------------------------ ###

library(mse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(doParallel)
library(scales)
library(cowplot)
library(patchwork)
library(RColorBrewer)
library(akima)
library(GA)
source("funs.R")
source("funs_GA.R")

# cl <- makeCluster(3)
# registerDoParallel(cl)
# clusterEvalQ(cl = cl, expr = {library(mse)})

### stock list
stocks <- read.csv("input/stocks.csv", stringsAsFactors = FALSE)
stocks_subset <- stocks$stock[21:29]
names(stocks_subset) <- stocks_subset
### brps
brps <- readRDS("input/brps.rds")

### ------------------------------------------------------------------------ ###
### collate data ####
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### collate data: pollack HR components ####
pol_GA <- foreach(stock = "pol", .combine = bind_rows) %:%
  #foreach(optimised = c(TRUE, FALSE), .combine = rbind) %:%
  foreach(objective = c("MSY", "MSY-PA"), .combine = bind_rows) %:%
  foreach(fhist = c("one-way", "random"), .combine = bind_rows) %:% 
  foreach(parameters = c("none",
                         "multiplier", 
                         "idx_lag", "idx_range", "b_multiplier", 
                         "interval",
                         "upper_cap", "lower_cap", "caps", 
                         "all_no_caps", "all", 
                         "mult_cond_cap", "all_cond_cap"), 
          .combine = bind_rows) %:%
  foreach(scenario = c("GA"), 
          .combine = bind_rows) %:%
  foreach(catch_rule = "hr", .combine = bind_rows) %:%
  foreach(stat_yrs = "all", .combine = bind_rows) %do% {#browser()
    
    ### find files
    file_prefix <- switch(parameters,
      "all" = paste0("idxB_lag-idxB_range_3-comp_b_multiplier-interval-",
                     "multiplier-upper_constraint-lower_constraint"),
      "all_no_caps" = paste0("idxB_lag-idxB_range_3-comp_b_multiplier-interval",
                             "-multiplier"),
      "caps" = "upper_constraint-lower_constraint",
      "upper_cap" = "upper_constraint",
      "lower_cap" = "lower_constraint",
      "multiplier" = "multiplier",
      "interval" = "interval",
      "idx_range" = "idxB_range_3",
      "idx_lag" = "idxB_lag",
      "b_multiplier" = "comp_b_multiplier",
      "none" = "multiplier",
      "mult_cond_cap" = "multiplier-upper_constraint1.2-lower_constraint0.7", 
      "all_cond_cap" = paste0("idxB_lag-idxB_range_3-comp_b_multiplier-",
                              "interval-multiplier-upper_constraint1.2-",
                              "lower_constraint0.7")
    )
    file_suffix <- switch(objective,
                          "MSY" = paste0("obj_SSB_C_risk_ICV"),
                          "MSY-PA" = "obj_ICES_MSYPA")
    file_res <- paste0(file_prefix, "--", file_suffix, "_res.rds")
    file_runs <- paste0(file_prefix, "--", file_suffix, "_runs.rds")
    scenario <- ifelse(parameters %in% c("mult_cond_cap", "all_cond_cap"),
                       "GA_cond_cap", scenario)
    
    if (!file.exists(paste0("output/", catch_rule, "/500_50/", 
                            scenario, "/", fhist, "/", stock, "/", file_res))) { 
      return(NULL)
    }
    
    df_res <- data.frame(parameters = parameters, objective = objective, 
                         scenario = scenario, catch_rule = catch_rule, 
                         stat_yrs = stat_yrs, fhist = fhist)
    ### load GA results
    res <- readRDS(paste0("output/", catch_rule, "/500_50/", 
                          scenario, "/", fhist, "/", stock, "/", file_res))
    runs <- readRDS(paste0("output/", catch_rule, "/500_50/", 
                           scenario, "/", fhist, "/", stock, "/", file_runs))
    solution <- res@solution
    ### if several equal best solutions, pick first
    if (isTRUE(nrow(solution) > 1)) {
      solution <- solution[1,, drop = FALSE]
      
    }
    solution[c(1, 2, 5)] <- round(solution[c(1, 2, 5)])
    solution[c(3, 4)] <- round(solution[c(3, 4)], 1)
    solution[c(6, 7, 8)] <- round(solution[c(6, 7, 8)], 2)
    
    ### default (non-optimised?)
    if (identical(parameters, "none")) {
      solution[] <- c(1, 1, 1, 1.4, 1, 1, Inf, 0)
    }
    df_res <- cbind(df_res, solution)
    solution_c <- paste0(solution, collapse = "_")
    ### replace NaN for upper_constraint with Inf
    solution_c <- gsub(x = solution_c, pattern = "NaN", replacement = "Inf")
    stats <- runs[[solution_c]]$stats
    stats <- as.data.frame(t(data.frame(unlist(stats), 
                                        row.names = rownames(stats))))
    stats <- stats %>% select(risk_Blim:ICV)
    row.names(stats) <- NULL
    df_res <- cbind(df_res, stats)
    ### prepare risk penalty
    penalty_tmp <- penalty(x = stats$risk_Blim, 
                           negative = FALSE, max = 5, 
                           inflection = 0.05 + 0.01, 
                           steepness = 0.5e+3)
    if (!identical(parameters, "none")) {
      df_res$fitness <- res@fitnessValue
      df_res$iter <- res@iter
    } else {
      if (isTRUE(objective == "MSY")) {
        df_res$fitness <- -sum(abs(stats$SSB_rel - 1), abs(stats$Catch_rel - 1), 
                               stats$risk_Blim, stats$ICV)
      } else if (isTRUE(objective == "MSY-PA")) {
        
        df_res$fitness <- -sum(abs(stats$SSB_rel - 1), 
                               abs(stats$Catch_rel - 1), 
                               stats$ICV,
                               penalty_tmp)
      }
      df_res$iter <- NA
    }
    df_res$risk_penalty <- ifelse(isTRUE(objective == "MSY-PA"),
                                  penalty_tmp, NA)
    return(df_res)
}
pol_GA <- pol_GA %>%
  mutate(parameters = factor(as.character(parameters), 
    levels = c("none", 
               "multiplier", 
               "idx_lag", "idx_range", "b_multiplier", 
               "interval", 
               "upper_cap", "lower_cap", "caps",
               "all_no_caps", "all", "mult_cond_cap", "all_cond_cap"),
    labels = c("not optimised", 
               "multiplier", 
               "time lag", "index range", "index trigger\nbuffer",
               "interval", 
               "upper cap", "lower cap", "both caps", 
               "all without caps", "all",
               "multiplier (cond. cap)",
               "all (cond. cap)")),
  ) %>%
  mutate(fhist = factor(as.character(fhist), 
                        levels = c("one-way", "random")))

saveRDS(pol_GA, file = "output/hr_pol_comps_stats.rds")
write.csv(pol_GA, file = "output/hr_pol_comps_stats.csv", 
          row.names = FALSE)

### ------------------------------------------------------------------------ ###
### collate data: all stocks ####
all_GA <- foreach(stock = stocks$stock, .combine = bind_rows) %:%
  foreach(parameters = c("none", "mult", "all"), .combine = bind_rows) %:%
  foreach(objective = c("MSY", "MSY-PA"), .combine = bind_rows) %:%
  foreach(scenario = c("GA", "GA_cond_cap"), .combine = bind_rows) %:%
  foreach(catch_rule = "hr", .combine = bind_rows) %:%
  foreach(stat_yrs = "all", .combine = bind_rows) %:%
  foreach(fhist = c("one-way", "random"), .combine = bind_rows) %do% {#browser()
    
    ### find files
    file_prefix <- switch(parameters,
      "none" = "multiplier",
      "mult" = "multiplier",
      "all" = paste0("idxB_lag-idxB_range_3-comp_b_multiplier-interval-",
                     "multiplier-upper_constraint-lower_constraint")
    )
    if (identical(scenario, "GA_cond_cap")) {
      if (identical(parameters, "all")) {
        file_prefix <- gsub(x = file_prefix, 
                            pattern = "upper_constraint-lower_constraint",
                        replacement = "upper_constraint1.2-lower_constraint0.7")
      } else if (isTRUE(parameters %in% c("none", "mult"))) {
        file_prefix <- paste0(file_prefix, "-", 
                              "upper_constraint1.2-lower_constraint0.7")
      }
    }
    file_suffix <- switch(objective,
                          "MSY" = paste0("obj_SSB_C_risk_ICV"),
                          "MSY-PA" = "obj_ICES_MSYPA")
    file_res <- paste0(file_prefix, "--", file_suffix, "_res.rds")
    file_runs <- paste0(file_prefix, "--", file_suffix, "_runs.rds")
    
    if (!file.exists(paste0("output/", catch_rule, "/500_50/", 
                            scenario, "/", fhist, "/", stock, "/", file_res))) { 
      return(NULL)
    }
    
    df_res <- data.frame(parameters = parameters, objective = objective, 
                         scenario = scenario, catch_rule = catch_rule, 
                         stat_yrs = stat_yrs, fhist = fhist, stock = stock)
    ### load GA results
    res <- readRDS(paste0("output/", catch_rule, "/500_50/", 
                          scenario, "/", fhist, "/", stock, "/", file_res))
    runs <- readRDS(paste0("output/", catch_rule, "/500_50/", 
                           scenario, "/", fhist, "/", stock, "/", file_runs))
    solution <- res@solution
    ### if several equal best solutions, pick first
    if (isTRUE(nrow(solution) > 1)) {
      solution <- solution[1,, drop = FALSE]
      
    }
    solution[c(1, 2, 5)] <- round(solution[c(1, 2, 5)])
    solution[c(3, 4)] <- round(solution[c(3, 4)], 1)
    solution[c(6, 7, 8)] <- round(solution[c(6, 7, 8)], 2)
    
    ### default (non-optimised?)
    if (identical(parameters, "none")) {
      solution[6] <- 1 ### multiplier
    }
    df_res <- cbind(df_res, solution)
    solution_c <- paste0(solution, collapse = "_")
    ### replace NaN for upper_constraint with Inf
    solution_c <- gsub(x = solution_c, pattern = "NaN", replacement = "Inf")
    stats <- runs[[solution_c]]$stats
    stats <- as.data.frame(t(data.frame(unlist(stats), 
                                        row.names = rownames(stats))))
    stats <- stats %>% select(risk_Blim:ICV)
    row.names(stats) <- NULL
    df_res <- cbind(df_res, stats)
    ### prepare risk penalty
    penalty_tmp <- penalty(x = stats$risk_Blim, 
                           negative = FALSE, max = 5, 
                           inflection = 0.05 + 0.01, 
                           steepness = 0.5e+3)
    if (!identical(parameters, "none")) {
      df_res$fitness <- res@fitnessValue
      df_res$iter <- res@iter
    } else {
      if (isTRUE(objective == "MSY")) {
        df_res$fitness <- -sum(abs(stats$SSB_rel - 1), abs(stats$Catch_rel - 1), 
                               stats$risk_Blim, stats$ICV)
      } else if (isTRUE(objective == "MSY-PA")) {
        df_res$fitness <- -sum(abs(stats$SSB_rel - 1), 
                               abs(stats$Catch_rel - 1), 
                               stats$ICV,
                               penalty_tmp)
      }
      df_res$iter <- NA
    }
    df_res$risk_penalty <- ifelse(isTRUE(objective == "MSY-PA"),
                                  penalty_tmp, NA)
    return(df_res)
}
saveRDS(all_GA, file = "output/hr_all_GA_stats.rds")
write.csv(all_GA, file = "output/hr_all_GA_stats.csv", row.names = FALSE)

### ------------------------------------------------------------------------ ###
### collate data: all stocks - multiplier with cond. cap - all runs ####
all_mult <- foreach(stock = stocks$stock, .combine = bind_rows) %:%
  foreach(parameters = c("mult"), .combine = bind_rows) %:%
  foreach(objective = c("MSY-PA"), .combine = bind_rows) %:%
  foreach(scenario = c("GA_cond_cap", "GA_cond_cap_mult"), .combine = bind_rows) %:%
  foreach(catch_rule = "hr", .combine = bind_rows) %:%
  foreach(stat_yrs = "all", .combine = bind_rows) %:%
  foreach(fhist = c("one-way", "random"), .combine = bind_rows) %:%
  foreach(years = c(50, 100), .combine = bind_rows) %do% {#browser()
    
    ### find files
    file_prefix <- "multiplier-upper_constraint1.2-lower_constraint0.7"
    file_suffix <- "obj_ICES_MSYPA"
    file_res <- paste0(file_prefix, "--", file_suffix, "_res.rds")
    file_runs <- paste0(file_prefix, "--", file_suffix, "_runs.rds")
    
    if (!file.exists(paste0("output/", catch_rule, "/500_", years, "/", 
                            scenario, "/", fhist, "/", stock, "/", file_res))) { 
      return(NULL)
    }
    
    df_res <- data.frame(parameters = parameters, objective = objective, 
                         scenario = scenario, catch_rule = catch_rule, 
                         stat_yrs = stat_yrs, fhist = fhist, stock = stock,
                         years = years)
    ### load GA results
    res <- readRDS(paste0("output/", catch_rule, "/500_", years, "/", 
                          scenario, "/", fhist, "/", stock, "/", file_res))
    runs <- readRDS(paste0("output/", catch_rule, "/500_", years, "/", 
                           scenario, "/", fhist, "/", stock, "/", file_runs))
    solution <- c(res@solution[, "multiplier"])
    
    runs_df <- lapply(runs, function(x) {
      tmp <- cbind(as.data.frame(t(x$pars)), 
                   as.data.frame(t(unlist(x$stats[1:11, ]))))
      penalty_tmp <- penalty(x = tmp$risk_Blim, 
                             negative = FALSE, max = 5, 
                             inflection = 0.05 + 0.01, 
                             steepness = 0.5e+3)
      if (isTRUE(objective == "MSY")) {
        tmp$fitness <- -sum(abs(tmp$SSB_rel - 1), abs(tmp$Catch_rel - 1), 
                            tmp$risk_Blim, tmp$ICV)
      } else if (isTRUE(objective == "MSY-PA")) {
        tmp$fitness <- -sum(abs(tmp$SSB_rel - 1),
                            abs(tmp$Catch_rel - 1), 
                            tmp$ICV,
                            penalty_tmp)
      }
      return(tmp)
    })
    runs_df <- do.call(rbind, runs_df)
    row.names(runs_df) <- NULL
    runs_df$optimum <- ifelse(round(runs_df$multiplier, 2) == round(solution, 2), 
                              TRUE, FALSE)
    
    df_res <- cbind(df_res, runs_df)
    return(df_res)
}
saveRDS(all_mult, file = "output/hr_all_mult_runs_stats.rds")
write.csv(all_mult, file = "output/hr_all_mult_runs_stats.csv", row.names = FALSE)

### ------------------------------------------------------------------------ ###
### collate data: all stocks - multiplier without cap - all runs ####
all_mult <- foreach(stock = stocks$stock, .combine = bind_rows) %:%
  foreach(parameters = c("mult"), .combine = bind_rows) %:%
  foreach(objective = c("MSY-PA"), .combine = bind_rows) %:%
  foreach(scenario = c("GA"), .combine = bind_rows) %:%
  foreach(catch_rule = "hr", .combine = bind_rows) %:%
  foreach(stat_yrs = "all", .combine = bind_rows) %:%
  foreach(fhist = c("one-way", "random", "roller-coaster"), 
          .combine = bind_rows) %:%
  foreach(years = c(50), .combine = bind_rows) %do% {#browser()
    
    ### find files
    file_prefix <- "multiplier"
    file_suffix <- "obj_ICES_MSYPA"
    file_res <- paste0(file_prefix, "--", file_suffix, "_res.rds")
    file_runs <- paste0(file_prefix, "--", file_suffix, "_runs.rds")
    
    if (!file.exists(paste0("output/", catch_rule, "/500_", years, "/", 
                            scenario, "/", fhist, "/", stock, "/", file_res))) { 
      return(NULL)
    }
    
    df_res <- data.frame(parameters = parameters, objective = objective, 
                         scenario = scenario, catch_rule = catch_rule, 
                         stat_yrs = stat_yrs, fhist = fhist, stock = stock,
                         years = years)
    ### load GA results
    res <- readRDS(paste0("output/", catch_rule, "/500_", years, "/", 
                          scenario, "/", fhist, "/", stock, "/", file_res))
    runs <- readRDS(paste0("output/", catch_rule, "/500_", years, "/", 
                           scenario, "/", fhist, "/", stock, "/", file_runs))
    solution <- c(res@solution[, "multiplier"])
    
    runs_df <- lapply(runs, function(x) {
      tmp <- cbind(as.data.frame(t(x$pars)), 
                   as.data.frame(t(unlist(x$stats[1:11, ]))))
      penalty_tmp <- penalty(x = tmp$risk_Blim, 
                             negative = FALSE, max = 5, 
                             inflection = 0.05 + 0.01, 
                             steepness = 0.5e+3)
      if (isTRUE(objective == "MSY")) {
        tmp$fitness <- -sum(abs(tmp$SSB_rel - 1), abs(tmp$Catch_rel - 1), 
                            tmp$risk_Blim, tmp$ICV)
      } else if (isTRUE(objective == "MSY-PA")) {
        tmp$fitness <- -sum(abs(tmp$SSB_rel - 1),
                            abs(tmp$Catch_rel - 1), 
                            tmp$ICV,
                            penalty_tmp)
      }
      return(tmp)
    })
    runs_df <- do.call(rbind, runs_df)
    row.names(runs_df) <- NULL
    runs_df$optimum <- ifelse(round(runs_df$multiplier, 2) == round(solution, 2), 
                              TRUE, FALSE)
    
    df_res <- cbind(df_res, runs_df)
    return(df_res)
}
saveRDS(all_mult, file = "output/hr_all_mult_runs_stats_nocaps.rds")
write.csv(all_mult, file = "output/hr_all_mult_runs_stats_nocaps.csv", 
          row.names = FALSE)

all_mult %>% filter(stock == "pol") %>%
  ggplot(aes(x = multiplier, y = risk_Blim, colour = as.factor(fhist))) +
  geom_line()
all_mult %>% filter(stock == "pol") %>%
  ggplot(aes(x = multiplier, y = Catch_rel, colour = as.factor(fhist))) +
  geom_line()


### ------------------------------------------------------------------------ ###
### figures for paper ####
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### visualise calculation of target HR - pollack ####
### ------------------------------------------------------------------------ ###

stock <- "pol"
input <- readRDS(paste0("input/hr/500_50/OM_2_mp_input/random/", stock, ".rds"))
lhist <- stocks[stocks$stock == stock, ]
brp <- brps$pol
refpts(brp)

### find iteration with increasing F
pos <- which(fbar(input$om@stock)[, ac(50)] < 0.05 * c(refpts(brp)["crash","harvest"]) &
               fbar(input$om@stock)[, ac(100)] > 0.95 * c(refpts(brp)["crash","harvest"]))
### 478
plot(input$om@stock, iter = 478)

### recreate mean catch length without noise
Lc <- calc_lc(stk = input$om@stock[, ac(50:100)], 
              a = lhist$a, b = lhist$b)
pars_l <- FLPar(a = lhist$a,  b = lhist$b, Lc = Lc)
idxL = lmean(stk = input$om@stock[, ac(50:100)], params = pars_l)
LFeM <- (lhist$linf + 2*1.5*c(Lc)) / (1 + 2*1.5)

plot(idxL, iter = 478) + geom_hline(yintercept = LFeM)

df_idxL <- as.data.frame(idxL) %>%
  filter(iter == pos)  %>%
  mutate(year = year - 100)
df_catch <- as.data.frame(catch(input$om@stock[, ac(50:100)])) %>%
  filter(iter == pos)  %>%
  mutate(year = year - 100)
df_idxB <- as.data.frame(input$oem@observations$idx$idxB[, ac(50:100)]) %>%
  filter(iter == pos)  %>%
  mutate(year = year - 100)
df_cr <- as.data.frame(catch(input$om@stock[, ac(50:100)]) / 
                         input$oem@observations$idx$idxB[, ac(50:100)]) %>%
  filter(iter == pos)  %>%
  mutate(year = year - 100)
### combine dfs and save
df_list <- list(idxL = df_idxL, catch = df_catch, idxB = df_idxB, cr = df_cr)
saveRDS(df_list, file = "output/hr_data_pol_visualisation.rds")
df_list <- readRDS("output/hr_data_pol_visualisation.rds")
df_idxL <- df_list$idxL
df_catch <- df_list$catch
df_idxB <- df_list$idxB
df_cr <- df_list$cr


### find years where L > LFeM
df_years <- df_idxL %>% filter(data >= LFeM) %>% select(year) %>% unlist()

### mean length
p_idxL1 <- df_idxL %>% 
  ggplot(aes(x = year, y = data)) +
  geom_line() +
  facet_wrap(~ "step 1", strip.position = "top") +
  labs(x = "year", y = "mean catch length [cm]") +
  theme_bw(base_size = 8) +
  theme(strip.text = element_text(size = 8))
p_idxL1
### add reference length
p_idxL2 <- df_idxL %>% 
  ggplot(aes(x = year, y = data)) +
  geom_line() +
  geom_hline(yintercept = LFeM, linetype = "dashed") + 
  facet_wrap(~ "step 2", strip.position = "top") +
  labs(x = "year", y = "mean catch length [cm]") +
  theme_bw(base_size = 8) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_text(size = 8))
p_idxL2
### add points
p_idxL3 <- df_idxL %>% 
  ggplot(aes(x = year, y = data)) +
  geom_line() +
  geom_hline(yintercept = LFeM, linetype = "dashed") + 
  geom_point(data = df_idxL %>% filter(year %in% df_years),
             colour = "red", size = 0.4) +
  facet_wrap(~ "step 3", strip.position = "top") +
  labs(x = "year", y = "mean catch length [cm]") +
  theme_bw(base_size = 8) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_text(size = 8))
p_idxL3

### catch
p_catch <- df_catch %>% 
  ggplot(aes(x = year, y = data)) +
  geom_line() +
  labs(x = "year", y = "catch") +
  theme_bw(base_size = 8) +
  theme(axis.title.x = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent", colour = NA))
p_catch
### biomass index
p_idxB <- df_idxB %>% 
  ggplot(aes(x = year, y = data)) +
  geom_line() +
  labs(x = "year", y = "index") +
  theme_bw(base_size = 8) +
  theme(axis.title.x = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent", colour = NA))
p_idxB

### catch rate
p_cr1 <- df_cr %>% 
  ggplot(aes(x = year, y = data)) +
  geom_line() +
  facet_wrap(~ "step 4", strip.position = "top") +
  labs(x = "year", y = "catch/index") +
  ylim(c(0, 0.4)) +
  theme_bw(base_size = 8) +
  theme(strip.text = element_text(size = 8)) +
  annotation_custom(grob = ggplotGrob(p_catch),
                    xmin = -50, xmax = -35,
                    ymin = 0.3, ymax = 0.4) +
  annotation_custom(grob = ggplotGrob(p_idxB),
                    xmin = -30, xmax = -15,
                    ymin = 0.3, ymax = 0.4) +
  annotate("segment", x = -34, xend = -31, y = 0.3, yend = 0.4,
           colour = "black", size = 1)
p_cr1
### catch rate & points
p_cr2 <- df_cr %>% 
  ggplot(aes(x = year, y = data)) +
  geom_line() +
  geom_point(data = df_cr %>% filter(year %in% df_years),
             colour = "red", size = 0.4) +
  facet_wrap(~ "step 5", strip.position = "top") +
  labs(x = "year", y = "catch/index") +
  ylim(c(0, 0.4)) +
  theme_bw(base_size = 8) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_text(size = 8))
p_cr2
### catch rate & points & reference
p_cr3 <- df_cr %>% 
  ggplot(aes(x = year, y = data)) +
  geom_line() +
  geom_point(data = df_cr %>% filter(year %in% df_years),
             colour = "red", size = 0.4) +
  geom_hline(yintercept = df_cr %>% 
               filter(year %in% df_years) %>% 
               summarise(mean(data)) %>% unlist(),
             colour = "red") +
  facet_wrap(~ "step 6", strip.position = "top") +
  labs(x = "year", y = "catch/index") +
  ylim(c(0, 0.4)) +
  theme_bw(base_size = 8) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_text(size = 8))
p_cr3


### combine all plots
plot_grid(p_idxL1, p_idxL2, p_idxL3,
          ncol = 3, rel_widths = c(1.12, 1, 1))
ggsave(filename = "output/plots/length_procedure_combined_1.png",
       width = 17, height = 6, units = "cm", dpi = 600,
       type = "cairo")
ggsave(filename = "output/plots/length_procedure_combined_1.pdf",
       width = 17, height = 6, units = "cm")
ggsave(filename = "output/plots/length_procedure_combined_1_.pdf",
       width = 17, height = 6, units = "cm")
plot_grid(p_cr1, p_cr2, p_cr3,
          ncol = 3, rel_widths = c(1.12, 1, 1))
ggsave(filename = "output/plots/length_procedure_combined_2.png",
       width = 17, height = 6, units = "cm", dpi = 600,
       type = "cairo")
ggsave(filename = "output/plots/length_procedure_combined_2.pdf",
       width = 17, height = 6, units = "cm")

### ------------------------------------------------------------------------ ###
### HR principle visualisation ####
### ------------------------------------------------------------------------ ###

### no colours
data.frame(x = c(0, 1, 2),
           y = c(0, 1, 1)) %>%
  ggplot(aes(x = x, y = y)) +
  geom_line() +
  scale_x_continuous("Index I", expand = c(0, 0),
                     breaks = c(0, 1), labels = c(0, expression(I[trigger]))) +
  scale_y_continuous("Harvest rate", limits = c(0, 1.2), expand = c(0, 0),
                     breaks = c(0, 1), 
                     labels = c(0, expression(H[target]))) +
  annotate(geom = "segment", x = 1, xend = 1, y = 0, yend = 1,
           linetype = "dotted") +
  annotate(geom = "segment", x = 0, xend = 1, y = 1, yend = 1,
           linetype = "dotted") +
  theme_classic()
ggsave(filename = "output/plots/paper/HR_principle_bw.png",
       width = 8.5, height = 5, units = "cm", dpi = 600,
       type = "cairo")
ggsave(filename = "output/plots/paper/HR_principle_bw.pdf",
       width = 8.5, height = 5, units = "cm")

### ------------------------------------------------------------------------ ###
### HR principle: plot catch and index selectivity ####
### ------------------------------------------------------------------------ ###

### extract fishery selectivity and maturity for all stocks
df_sel <- foreach(stock = stocks$stock, .combine = bind_rows) %do% {
  #browser()
  brp_i <- brps[[stock]]
  df_i <- as.data.frame(FLQuants(maturity = brp_i@mat,
                                 selectivity = brp_i@landings.sel/
                                   max(brp_i@landings.sel))) %>%
    select(age, data, qname) %>%
    mutate(stock = stock)
  return(df_i)
}
df_sel <- df_sel %>%
  full_join(stocks %>% 
              select(stock, k) %>%
              mutate(ID = 1:29)) %>%
  mutate(stock_id = paste0(stock, "~(italic(k)==", sprintf(k, fmt =  "%.2f"), 
                           ")")) 
df_sel <- df_sel %>%
  mutate(stock_id = factor(stock_id, levels = unique(df_sel$stock_id))) %>%
  mutate(qname = factor(qname, levels = c("maturity", "selectivity"))) %>%
  rename(value = data)

df_sel %>%
  ggplot(aes(x = age, y = value, colour = qname)) +
  geom_line(size = 0.3) +
  facet_wrap(~ stock_id, scales = "free_x", labeller = "label_parsed", 
             ncol = 5) +
  labs(x = "Age [years]", y = "Maturity | Selectivity") + 
  scale_x_continuous(breaks = pretty_breaks(n = 3), limits = c(0, NA)) +
  scale_colour_discrete("") +
  theme_bw(base_size = 8) +
  theme(legend.position = c(0.9, 0.075),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.5, "lines"))
ggsave(filename = "output/plots/paper/HR_OMs_sel_mat.png",
       width = 16, height = 16, units = "cm", dpi = 600,
       type = "cairo")
ggsave(filename = "output/plots/paper/HR_OMs_sel_mat.pdf",
       width = 16, height = 16, units = "cm")


### ------------------------------------------------------------------------ ###
### sensitivity to simulation assumptions ####
### ------------------------------------------------------------------------ ###
### use pollack as example
### default: 50 yrs, 500 iterations

### some generic parameters
brp <- readRDS("input/brps.rds")[["pol"]]
Blim <- brp@Blim
Bmsy <- c(refpts(brp)["msy", "ssb"])
MSY <- c(refpts(brp)["msy", "yield"])

### stats over time
stats_sens_time <- foreach(fhist = c("random", "one-way", "roller-coaster"),
                           .combine = rbind) %do% {
                             #browser()
  file <- ifelse(fhist %in% c("random", "one-way"),
                 "mp_length_1_TRUE_1_1_1_Inf_0",
                 "mp_length_1_TRUE_1_1_1_Inf_0_0.2_0.2_0_0_0.6_0_0.75")
  res <- readRDS(paste0("output/hr/500_100/sensitivity/", fhist, 
                        "/pol/", file, ".rds"))
  
  ### collapse correction
  res_corrected <- collapse_correction(stk = res@stock, yrs = 101:200)
  ### template
  tmp <- data.frame(year = 1:100)
  ### Blim risk
  tmp$risk_average <- sapply(1:100, function(x) {
    mean(c(res_corrected$ssb[, ac(seq(from = 101, length.out = x))] < Blim), 
         na.rm = TRUE)
  })
  tmp$risk_annual <- c(apply(res_corrected$ssb < Blim, 2, mean, na.rm = TRUE))
  ### SSB
  tmp$SSB_annual <- sapply(1:100, function(x) {
    median(c(res_corrected$ssb[, x]/Bmsy), na.rm = TRUE)
  })
  tmp$SSB_average <-  sapply(1:100, function(x) {
    median(c(res_corrected$ssb[, ac(seq(from = 101, length.out = x))]/Bmsy), 
           na.rm = TRUE)
  })
  ### Catch
  tmp$Catch_annual <- sapply(1:100, function(x) {
    median(c(res_corrected$catch[, x]/MSY), na.rm = TRUE)
  })
  tmp$Catch_average <-  sapply(1:100, function(x) {
    median(c(res_corrected$catch[, ac(seq(from = 101, length.out = x))]/MSY), 
           na.rm = TRUE)
  })
  tmp <- tmp %>%
    pivot_longer(2:7, names_to = c(".value", "period"), names_sep = "_")
  ### full data.frame
  df_i <- data.frame(
    stock = "pol", multiplier = 1, comp_b = TRUE, idxB_lag = 1, 
    idxB_range_3 = 1, interval = 1, upper_constraint = Inf, 
    lower_constraint = 0, sigmaL = 0.2, sigmaB = 0.2, sigmaR = 0.6,
    sigmaR_rho = 0,
    risk_Blim = tmp$risk,
    SSB_rel = tmp$SSB,
    Catch_rel = tmp$Catch,
    stat_metric = tmp$period,
    fhist = fhist,
    n_yrs = tmp$year,
    n_iter = 10000,
    steepness = 0.75,
    sensitivity = "period") %>%
    arrange(stat_metric, n_yrs)
  return(df_i)
}

### stock status
stats_sens_status <- foreach(fhist = c("random", "one-way"),
                             .combine = rbind) %do% {
                               #browser()
  res <- readRDS(paste0("output/hr/10000_50/sensitivity/", fhist, 
                        "/pol/mp_length_1_TRUE_1_1_1_Inf_0.rds"))
  ### collapse correction
  res_corrected <- collapse_correction(stk = res@stock, yrs = 101:150)
  ### starting condition
  SSBs0 <- ssb(res@stock)[, ac(100)]
  SSBs0 <- SSBs0/Bmsy
  SSBs0 <- c(SSBs0)
  SSB_breaks <- seq(from = 0, to = max(SSBs0), by = 0.1)
  SSB_groups <- cut(SSBs0, breaks = SSB_breaks)
  SSB_levels <- unique(as.character(SSB_groups))
  ### number of replicates per group
  group_n <- sapply(SSB_levels, function(x) {
    length(which(SSB_groups %in% x))
  })
  # group_n[sort(names(group_n))]
  ### Blim risk per group
  ### SSB is on absolute scale 
  risk_group <- sapply(SSB_levels, function(x) {
    tmp <- res_corrected$ssb[,,,,, which(SSB_groups %in% x)]
    mean(tmp < (Blim))
  })
  ### SSB (long-term median) per group
  SSB_group <- sapply(SSB_levels, function(x) {
    tmp <- res_corrected$ssb[,,,,, which(SSB_groups %in% x)]/Bmsy
    median(tmp)
  })
  ### Catch (long-term median) per group
  Catch_group <- sapply(SSB_levels, function(x) {
    tmp <- res_corrected$catch[,,,,, which(SSB_groups %in% x)]/MSY
    median(tmp)
  })
  ### get starting conditions
  SSB_levels <- sapply(SSB_levels, function(x) {
    x <- gsub(x = x, pattern = "\\(|\\]", replacement = "")
    x <- unlist(strsplit(x, split = ","))
    mean(as.numeric(x))
  })
  pos_remove <- which(is.na(SSB_levels))
  df_i <- data.frame(
    stock = "pol", multiplier = 1, comp_b = TRUE, idxB_lag = 1, 
    idxB_range_3 = 1, interval = 1, upper_constraint = Inf, 
    lower_constraint = 0, sigmaL = 0.2, sigmaB = 0.2, sigmaR = 0.6,
    sigmaR_rho = 0,
    risk_Blim = unlist(risk_group)[-pos_remove],
    SSB_rel = unlist(SSB_group)[-pos_remove],
    Catch_rel = unlist(Catch_group)[-pos_remove],
    SSB0_rel = unlist(SSB_levels)[-pos_remove],
    n_iter_part = unlist(group_n)[-pos_remove],
    fhist = fhist,
    n_yrs = 50,
    n_iter = 10000,
    steepness = 0.75,
    sensitivity = "stock_status")
  row.names(df_i) <- NULL
  df_i <- df_i[order(df_i$SSB0_rel), ]
  return(df_i)
}

### recruitment & observations
stats_sens_more <- foreach(sensitivity = c("rec_var", "rec_rho",
                                           "rec_steepness",
                                           #"obs_idx_sd", "obs_lngth_sd",
                                           "obs_idx_lngth_sd",
                                           "obs_idx_lngth_rho"),
                           .combine = bind_rows) %:%
  foreach(fhist = c("one-way", "roller-coaster", "random"), 
          .combine = bind_rows) %do% {
            #browser()
  if (fhist %in% c("one-way", "random")) {
    file_name <- switch(sensitivity,
      "rec_var" = "1_1_1_1_1_Inf_0_0.2_0.2_0-1_0.rds",
      "rec_rho" = "1_1_1_1_1_Inf_0_0.2_0.2_0.6_0-1_0.75.rds",
      "rec_steepness" = "1_1_1_1_1_Inf_0_0.2_0.2_0.6_0_0-1.rds",
      "obs_idx_sd" = "1_1_1_1_1_Inf_0_0.2_0-1_0.6_0.rds",
      "obs_lngth_sd" = "1_1_1_1_1_Inf_0_0-1_0.2_0.6_0.rds",
      "obs_idx_lngth_sd" = "1_1_1_1_1_Inf_0_0-1_0-1_0.6_0.rds",
      "obs_idx_lngth_rho" = "1_1_1_1_1_Inf_0_0.2_0.2_0-1_0-1_0.6_0_0.75.rds"
    )
  } else {
    file_name <- switch(sensitivity,
      "rec_var" = "1_1_1_1_1_Inf_0_0.2_0.2_0_0_0-1_0_0.75.rds",
      "rec_rho" = "1_1_1_1_1_Inf_0_0.2_0.2_0_0_0.6_0-1_0.75.rds",
      "rec_steepness" = "1_1_1_1_1_Inf_0_0.2_0.2_0_0_0.6_0_0-1.rds",
      "obs_idx_sd" = "1_1_1_1_1_Inf_0_0.2_0-1_0_0_0.6_0_0.75.rds",
      "obs_lngth_sd" = "1_1_1_1_1_Inf_0_0-1_0.2_0_0_0.6_0_0.75.rds",
      "obs_idx_lngth_sd" = "1_1_1_1_1_Inf_0_0-1_0-1_0_0_0.6_0_0.75.rds",
      "obs_idx_lngth_rho" = "1_1_1_1_1_Inf_0_0.2_0.2_0-1_0-1_0.6_0_0.75.rds"
    )
  }
  
  file_name <- paste0("collated_stats_length_", file_name)
  tmp <- readRDS(paste0("output/hr/500_50/sensitivity/", fhist, "/pol/",
                        file_name))
  tmp <- as.data.frame(lapply(tmp, unlist))
  tmp$sensitivity <- sensitivity
  tmp$fhist <- fhist
  tmp$n_yrs <- 50
  tmp$n_iter <- 500
  return(tmp)
}

### combine all sensitivity runs
stats_sens <- bind_rows(stats_sens_more, stats_sens_status, stats_sens_time)
saveRDS(stats_sens, "output/hr_pol_sensitivity.rds")
write.csv(stats_sens, "output/hr_pol_sensitivity.csv", row.names = FALSE)
stats_sens <- readRDS("output/hr_pol_sensitivity.rds")

### plot sensitivity analysis
stats_sens_plot <- stats_sens %>%
  pivot_longer(c(risk_Blim, SSB_rel, Catch_rel)) %>%
  mutate(name = factor(name, levels = c("SSB_rel", "Catch_rel", "risk_Blim"),
                       labels = c("SSB/B[MSY]", "Catch/MSY", "B[lim]~risk")),
         fhist = factor(fhist, levels = c("one-way", "roller-coaster",
                                          "random")))
df_blank <- data.frame(name = rep(c("SSB/B[MSY]", "Catch/MSY", "B[lim]~risk"),
                                  each = 2),
                       x = c(0, 1, 0, 1, 0, 1),
                       value = c(0, 1.8, 0, 1.4, 0, 1),
                       fhist = NA)
res_def_colours <- c("one-way" = brewer.pal(n = 4, name = "Set1")[1], 
                     "roller-coaster" = brewer.pal(n = 4, name = "Set1")[4], 
                     "random" = brewer.pal(n = 4, name = "Set1")[2])
res_def_linetype <- c("one-way" = "solid", 
                      "roller-coaster" = "1212", 
                      "random" = "3232")
p_sens_rec <- stats_sens_plot %>%
  filter(sensitivity == "rec_var") %>%
  ggplot(aes(x = sigmaR, y = value, fill = fhist, colour = fhist, 
             linetype = fhist)) +
  geom_vline(xintercept = 0.6, size = 0.4, colour = "grey") +
  stat_smooth(n = 50, span = 0.2, se = FALSE, geom = "line", size = 0.4,
              show.legend = FALSE) + 
  geom_point(size = 0.15, stroke = 0, shape = 21, show.legend = FALSE) +
  geom_blank(data = df_blank, aes(x = x, y = value)) +
  facet_grid(name ~ "'Recruitment\n\ \ variability'", scales = "free", 
             labeller = "label_parsed",
             switch = "y") +
  scale_linetype_manual("fishing history", values = res_def_linetype) +
  scale_colour_manual("fishing history", values = res_def_colours) +
  scale_fill_manual("fishing history", values = res_def_colours) +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  labs(x = expression(sigma[R])) +
  theme_bw(base_size = 8) +
  theme(strip.placement = "outside",
        strip.text.y = element_text(size = 8),
        strip.text.x = element_text(margin = margin(8, 0, 0, 0)),
        strip.background.y = element_blank(),
        axis.title.y = element_blank(),
        strip.switch.pad.grid = unit(0, "pt"),
        plot.margin = unit(c(2, 2, 4, 4), "pt"))
p_sens_rec_rho <- stats_sens_plot %>%
  filter(sensitivity == "rec_rho" &
           sigmaR_rho < 0.99) %>%
  ggplot(aes(x = sigmaR_rho, y = value, fill = fhist, colour = fhist, 
             linetype = fhist)) +
  geom_vline(xintercept = 0.0, size = 0.4, colour = "grey") +
  stat_smooth(n = 50, span = 0.2, se = FALSE, geom = "line", size = 0.4,
              show.legend = FALSE) + 
  geom_point(size = 0.15, stroke = 0, shape = 21, show.legend = FALSE) +
  geom_blank(data = df_blank, aes(x = x, y = value)) +
  facet_grid(name ~ "'\ \ \ Recruitment\nauto-correlation'", scales = "free", 
             labeller = "label_parsed",
             switch = "y") +
  scale_linetype_manual("fishing history", values = res_def_linetype) +
  scale_colour_manual("fishing history", values = res_def_colours) +
  scale_fill_manual("fishing history", values = res_def_colours) +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  labs(x = expression(italic(rho)[R])) +
  theme_bw(base_size = 8) +
  theme(strip.placement = "outside",
        strip.text.y = element_blank(),
        strip.text.x = element_text(margin = margin(8, 0, 1.2, 0)),
        strip.background.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        strip.switch.pad.grid = unit(0, "pt"),
        plot.margin = unit(c(2, 2, 4, 0), "pt"))
p_sens_rec_steepness <- stats_sens_plot %>%
  filter(sensitivity == "rec_steepness" &
           steepness > 0.2) %>%
  ggplot(aes(x = steepness, y = value, fill = fhist, colour = fhist, 
             linetype = fhist)) +
  geom_vline(xintercept = 0.75, size = 0.4, colour = "grey") +
  stat_smooth(n = 50, span = 0.2, se = FALSE, geom = "line", size = 0.4,
              show.legend = FALSE) + 
  geom_point(size = 0.15, stroke = 0, shape = 21, show.legend = FALSE) +
  geom_blank(data = df_blank, aes(x = x, y = value)) +
  facet_grid(name ~ "'Recruitment\n\ steepness'", scales = "free", 
             labeller = "label_parsed",
             switch = "y") +
  scale_linetype_manual("fishing history", values = res_def_linetype) +
  scale_colour_manual("fishing history", values = res_def_colours) +
  scale_fill_manual("fishing history", values = res_def_colours) +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  labs(x = expression(italic(h))) +
  theme_bw(base_size = 8) +
  theme(strip.placement = "outside",
        strip.text.y = element_blank(),
        strip.text.x = element_text(margin = margin(8, 0, 0, 0)),
        strip.background.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        strip.switch.pad.grid = unit(0, "pt"),
        plot.margin = unit(c(2, 2, 4, 0), "pt"))
p_sens_obs <- stats_sens_plot %>%
  filter(sensitivity == "obs_idx_lngth_sd") %>%
  ggplot(aes(x = sigmaB, y = value, fill = fhist, colour = fhist, 
             linetype = fhist)) +
  geom_vline(xintercept = 0.2, size = 0.4, colour = "grey") +
  stat_smooth(n = 50, span = 0.2, se = FALSE, geom = "line", size = 0.4,
              show.legend = FALSE) + 
  geom_point(size = 0.15, stroke = 0, shape = 21, show.legend = FALSE) +
  geom_blank(data = df_blank, aes(x = x, y = value)) +
  facet_grid(name ~ "'Observation\n\ uncertainty'", scales = "free", 
             labeller = "label_parsed",
             switch = "y") +
  scale_linetype_manual("fishing history", values = res_def_linetype) +
  scale_colour_manual("fishing history", values = res_def_colours) +
  scale_fill_manual("fishing history", values = res_def_colours) +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  labs(x = expression(italic(sigma)[obs])) +
  theme_bw(base_size = 8) +
  theme(strip.placement = "outside",
        strip.text.y = element_blank(),
        strip.text.x = element_text(margin = margin(8, 0, 0, 0)),
        strip.background.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        strip.switch.pad.grid = unit(0, "pt"),
        plot.margin = unit(c(2, 2, 4, 0), "pt"))
p_sens_obs_rho <- stats_sens_plot %>%
  filter(sensitivity == "obs_idx_lngth_rho" & sigmaB_rho < 1) %>%
  ggplot(aes(x = sigmaB_rho, y = value, fill = fhist, colour = fhist, 
             linetype = fhist)) +
  geom_vline(xintercept = 0.0, size = 0.4, colour = "grey") +
  stat_smooth(n = 50, span = 0.2, se = FALSE, geom = "line", size = 0.4,
              show.legend = FALSE) + 
  geom_point(size = 0.15, stroke = 0, shape = 21, show.legend = FALSE) +
  geom_blank(data = df_blank, aes(x = x, y = value)) +
  facet_grid(name ~ "'\ \ \ Observation\nauto-correlation'", scales = "free", 
             labeller = "label_parsed",
             switch = "y") +
  scale_linetype_manual("fishing history", values = res_def_linetype) +
  scale_colour_manual("fishing history", values = res_def_colours) +
  scale_fill_manual("fishing history", values = res_def_colours) +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  labs(x = expression(italic(rho)[obs])) +
  theme_bw(base_size = 8) +
  theme(strip.placement = "outside",
        strip.text.y = element_blank(),
        strip.text.x = element_text(margin = margin(8, 0, 1, 0)),
        strip.background.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        strip.switch.pad.grid = unit(0, "pt"),
        plot.margin = unit(c(2, 2, 4, 0), "pt"))
p_sens_status <- stats_sens_plot %>%
  filter(sensitivity == "stock_status" &
           SSB0_rel <= 2) %>%
  mutate(value = ifelse(fhist == "one-way", NA, value)) %>%
  ggplot(aes(x = SSB0_rel, y = value, fill = fhist, colour = fhist, 
             linetype = fhist)) +
  stat_smooth(n = 50, span = 0.4, se = FALSE, geom = "line", size = 0.4,
              show.legend = FALSE) + 
  geom_point(size = 0.15, stroke = 0, shape = 21, show.legend = FALSE) +
  geom_blank(data = df_blank, aes(x = x, y = value)) +
  facet_grid(name ~ "'\ \ \ \ \ \ Initial\nstock status'", scales = "free",
             labeller = "label_parsed",
             switch = "y") +
  scale_linetype_manual("fishing history", values = res_def_linetype) +
  scale_colour_manual("fishing history", values = res_def_colours) +
  scale_fill_manual("fishing history", values = res_def_colours) +
  scale_x_continuous(limits = c(-0.05, 2.05)) +
  labs(x = expression(SSB[y == 0]/B[MSY])) +
  theme_bw(base_size = 8) +
  theme(strip.placement = "outside",
        strip.text.y = element_blank(),
        ### manual margins because no letter goes below base
        #strip.text.x = element_text(margin = unit(c(3.8, 0, 3.8, 0), "pt")),
        strip.text.x = element_text(margin = margin(8, 0, 1.2, 0)),
        strip.background.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        strip.switch.pad.grid = unit(0, "pt"),
        plot.margin = unit(c(2, 2, 4, 0), "pt"))
p_sens_period <- stats_sens_plot %>%
  filter(sensitivity == "period" &
           stat_metric == "average") %>%
  ggplot(aes(x = n_yrs, y = value, fill = fhist, colour = fhist, 
             linetype = fhist)) +
  geom_vline(xintercept = 50, size = 0.4, colour = "grey") +
  stat_smooth(n = 50, span = 0.1, se = FALSE, geom = "line", size = 0.4) + 
  geom_point(size = 0.15, stroke = 0, shape = 21) +
  geom_blank(data = df_blank, aes(x = x, y = value)) +
  facet_grid(name ~ "'Implementation\n\ \ \ \ \ \ \ period'", scales = "free", 
             labeller = "label_parsed",
             switch = "y") +
  scale_linetype_manual("fishing history", values = res_def_linetype) +
  scale_colour_manual("fishing history", values = res_def_colours) +
  scale_fill_manual("fishing history", values = res_def_colours) +
  labs(x = "years") +
  theme_bw(base_size = 8) +
  theme(strip.placement = "outside",
        strip.text.y = element_blank(),
        strip.text.x = element_text(margin = margin(8, 0, 0, 0)),
        strip.background.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        strip.switch.pad.grid = unit(0, "pt"),
        plot.margin = unit(c(2, 4, 4, 0), "pt"),
        legend.position = c(0.55, 0.26),
        legend.background = element_blank(),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.6, "lines"),
        legend.title = element_blank(),
        legend.key = element_blank())

p <- p_sens_rec + p_sens_rec_rho + p_sens_rec_steepness +
  p_sens_obs + p_sens_obs_rho +
  p_sens_status + p_sens_period +
  plot_layout(nrow = 1)
p

ggsave(filename = "output/plots/paper/pol_sensitivity_stats.png", 
       type = "cairo", plot = p,
       width = 17, height = 8, units = "cm", dpi = 600)
ggsave(filename = "output/plots/paper/pol_sensitivity_stats.pdf", plot = p,
       width = 17, height = 8, units = "cm")

### ------------------------------------------------------------------------ ###
### sensitivity - index selectivity for pollack ####
### ------------------------------------------------------------------------ ###

### collate results
df_sel <- foreach(data = c("projection", "summary")) %:%
  foreach(stock = c("pol"), .combine = rbind) %:%
  foreach(fhist = c("random", "one-way", "roller-coaster"), .combine = rbind) %:%
  foreach(idx_sel = c("tsb", "ssb", "catch", "dome_shaped"), 
          .combine = rbind) %do% {
  #browser()
  res <- readRDS(paste0("output/hr/500_50/sensitivity_idx_sel/", fhist, "/",
                        stock, "/mp_length_1_TRUE_1_1_1_Inf_0_0.2_0.2_0_0_",
                        "0.6_0_0.75", 
                        ifelse(identical(idx_sel, "tsb"), "", 
                               paste0("_", idx_sel)),
                        ".rds"))
  ### collapse correction
  res_corrected <- collapse_correction(stk = res@stock, yrs = 100:150)
  ### get catch and ssb
  ssb_rel <- res_corrected$ssb/c(refpts(brps[[stock]])["msy", "ssb"])
  catch_rel <- res_corrected$catch/c(refpts(brps[[stock]])["msy", "yield"])
  ssb_rel_median <- iterMedians(ssb_rel)
  catch_rel_median <- iterMedians(catch_rel)
  
  res_desc <- data.frame(stock = stock, fhist = fhist, idx_sel = idx_sel, 
                         stringsAsFactors = FALSE)
  res_projection <- data.frame(year = as.numeric(dimnames(ssb_rel_median)$year),
                               catch = c(catch_rel_median),
                               ssb = c(ssb_rel_median))
  res_projection <- bind_cols(res_desc, res_projection)
  res_summary <- data.frame(catch = c(catch_rel), ssb = c(ssb_rel))
  res_summary <- bind_cols(res_desc, res_summary)
  
  if (identical(data, "projection")) {
    return(res_projection)
  } else if (identical(data, "summary")) {
    return(res_summary)
  }
}
saveRDS(df_sel, file = "output/hr_sensitivity_sel.rds")
df_sel <- readRDS("output/hr_sensitivity_sel.rds")

df_sel <- lapply(df_sel, function(x) {
  x <- x %>%
    mutate(idx_sel = factor(idx_sel, 
                            levels = c("ssb", "dome_shaped", "catch", "tsb"),
                            labels = c("maturity", "dome shaped", "commercial", 
                                       "total biomass")))
  return(x)
})
names(df_sel) <- c("projection", "summary")

#brewer.pal(n = 3, name = "Dark2")
cols_sel <- c("maturity" = brewer.pal(n = 3, name = "Dark2")[1], 
              "dome shaped" = brewer.pal(n = 3, name = "Dark2")[2], 
              "commercial" = brewer.pal(n = 3, name = "Dark2")[3], 
              "total biomass" = "azure4")

p_ssb_proj_ow <- df_sel$projection %>%
  filter(fhist == "one-way") %>%
  ggplot(aes(x = year - 100, y = ssb, colour = idx_sel)) +
  geom_line(size = 0.3) +
  facet_wrap(~ "Projection") +
  scale_colour_manual("", values = cols_sel) + 
  scale_x_continuous(breaks = c(0, 5, 10)) + 
  coord_cartesian(ylim = c(0, 2.5), xlim = c(0, 10.5),
                  expand = FALSE) +
  labs(x = "year",
       y = expression(SSB/B[MSY])) + 
  theme_bw(base_size = 8) +
  theme(legend.position = c(0.22, 0.80),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.6, "lines"),
        legend.background = element_blank(),
        legend.key = element_blank(),
        plot.margin = unit(c(4, 0.1, 4, 4), "pt"),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
p_ssb_smry_ow <- df_sel$summary %>%
  filter(fhist == "one-way") %>%
  mutate(x = 1) %>%
  ggplot(aes(y = ssb, x = x, fill = idx_sel, group = idx_sel)) +
  geom_boxplot(outlier.size = 0.3, size = 0.3,
               width = 1, position = position_dodge(1),
               show.legend = FALSE) +
  facet_wrap(~ "Summary") +
  scale_fill_manual(values = cols_sel) + 
  coord_cartesian(ylim = c(0, 2.5), xlim = c(0.4, 1.6),  
                  expand = FALSE) +
  theme_bw(base_size = 8) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.margin = unit(c(4, 4, 4, 0), "pt"))
p_catch_proj_ow <- df_sel$projection %>%
  filter(fhist == "one-way") %>%
  ggplot(aes(x = year - 100, y = catch, colour = idx_sel)) +
  geom_line(size = 0.3, show.legend = FALSE) +
  scale_x_continuous(breaks = c(0, 5, 10)) + 
  scale_colour_manual(values = cols_sel) + 
  coord_cartesian(ylim = c(0, 2.05), xlim = c(0, 10.5),
                  expand = FALSE) +
  labs(x = "year",
       y = expression(Catch/MSY)) + 
  theme_bw(base_size = 8) +
  theme(plot.margin = unit(c(0, 0.1, 4, 4), "pt"))
p_catch_smry_ow <- df_sel$summary %>%
  filter(fhist == "one-way") %>%
  mutate(x = 1) %>%
  ggplot(aes(y = catch, x = x, fill = idx_sel, group = idx_sel)) +
  geom_boxplot(outlier.size = 0.3, size = 0.3,
               width = 1, position = position_dodge(1),
               show.legend = FALSE) +
  scale_fill_manual(values = cols_sel) + 
  coord_cartesian(ylim = c(0, 2.05), xlim = c(0.4, 1.6),  
                  expand = FALSE) +
  theme_bw(base_size = 8) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.margin = unit(c(0, 4, 4, 0), "pt"))
    
p_ow <- p_ssb_proj_ow + p_ssb_smry_ow + 
  p_catch_proj_ow + p_catch_smry_ow + 
  plot_layout(ncol = 2, widths = c(1, 0.5))
p_ow
ggsave(filename = "output/plots/paper/pol_sensitivity_idx_sel_ow.png", 
       type = "cairo", plot = p_ow,
       width = 8.5, height = 6, units = "cm", dpi = 600)
ggsave(filename = "output/plots/paper/pol_sensitivity_idx_sel_ow.pdf", 
       plot = p_ow,
       width = 8.5, height = 6, units = "cm")

### same for roller-coaster
p_ssb_proj_rc <- df_sel$projection %>%
  filter(fhist == "roller-coaster") %>%
  ggplot(aes(x = year - 100, y = ssb, colour = idx_sel)) +
  geom_line(size = 0.3) +
  facet_wrap(~ "Projection") +
  scale_colour_manual("", values = cols_sel) + 
  scale_x_continuous(breaks = c(0, 5, 10)) + 
  coord_cartesian(ylim = c(0, 2.5), xlim = c(0, 10.5),
                  expand = FALSE) +
  labs(x = "year",
       y = expression(SSB/B[MSY])) + 
  theme_bw(base_size = 8) +
  theme(legend.position = c(0.22, 0.80),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.6, "lines"),
        legend.background = element_blank(),
        legend.key = element_blank(),
        plot.margin = unit(c(4, 0.1, 4, 4), "pt"),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
p_ssb_smry_rc <- df_sel$summary %>%
  filter(fhist == "roller-coaster") %>%
  mutate(x = 1) %>%
  ggplot(aes(y = ssb, x = x, fill = idx_sel, group = idx_sel)) +
  geom_boxplot(outlier.size = 0.3, size = 0.3,
               width = 1, position = position_dodge(1),
               show.legend = FALSE) +
  facet_wrap(~ "Summary") +
  scale_fill_manual(values = cols_sel) + 
  coord_cartesian(ylim = c(0, 2.5), xlim = c(0.4, 1.6),  
                  expand = FALSE) +
  theme_bw(base_size = 8) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.margin = unit(c(4, 4, 4, 0), "pt"))
p_catch_proj_rc <- df_sel$projection %>%
  filter(fhist == "roller-coaster") %>%
  ggplot(aes(x = year - 100, y = catch, colour = idx_sel)) +
  geom_line(size = 0.3, show.legend = FALSE) +
  scale_x_continuous(breaks = c(0, 5, 10)) + 
  scale_colour_manual(values = cols_sel) + 
  coord_cartesian(ylim = c(0, 2.05), xlim = c(0, 10.5),
                  expand = FALSE) +
  labs(x = "year",
       y = expression(Catch/MSY)) + 
  theme_bw(base_size = 8) +
  theme(plot.margin = unit(c(0, 0.1, 4, 4), "pt"))
p_catch_smry_rc <- df_sel$summary %>%
  filter(fhist == "roller-coaster") %>%
  mutate(x = 1) %>%
  ggplot(aes(y = catch, x = x, fill = idx_sel, group = idx_sel)) +
  geom_boxplot(outlier.size = 0.3, size = 0.3,
               width = 1, position = position_dodge(1),
               show.legend = FALSE) +
  scale_fill_manual(values = cols_sel) + 
  coord_cartesian(ylim = c(0, 2.05), xlim = c(0.4, 1.6),  
                  expand = FALSE) +
  theme_bw(base_size = 8) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.margin = unit(c(0, 4, 4, 0), "pt"))
p_rc <- p_ssb_proj_rc + p_ssb_smry_rc + 
  p_catch_proj_rc + p_catch_smry_rc + 
  plot_layout(ncol = 2, widths = c(1, 0.5))

### and for random
p_ssb_proj_rnd <- df_sel$projection %>%
  filter(fhist == "random") %>%
  ggplot(aes(x = year - 100, y = ssb, colour = idx_sel)) +
  geom_line(size = 0.3) +
  facet_wrap(~ "Projection") +
  scale_colour_manual("", values = cols_sel) + 
  scale_x_continuous(breaks = c(0, 5, 10)) + 
  coord_cartesian(ylim = c(0, 2.5), xlim = c(0, 10.5),
                  expand = FALSE) +
  labs(x = "year",
       y = expression(SSB/B[MSY])) + 
  theme_bw(base_size = 8) +
  theme(legend.position = c(0.22, 0.80),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.6, "lines"),
        legend.background = element_blank(),
        legend.key = element_blank(),
        plot.margin = unit(c(4, 0.1, 4, 4), "pt"),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
p_ssb_smry_rnd <- df_sel$summary %>%
  filter(fhist == "random") %>%
  mutate(x = 1) %>%
  ggplot(aes(y = ssb, x = x, fill = idx_sel, group = idx_sel)) +
  geom_boxplot(outlier.size = 0.3, size = 0.3,
               width = 1, position = position_dodge(1),
               show.legend = FALSE) +
  facet_wrap(~ "Summary") +
  scale_fill_manual(values = cols_sel) + 
  coord_cartesian(ylim = c(0, 2.5), xlim = c(0.4, 1.6),  
                  expand = FALSE) +
  theme_bw(base_size = 8) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.margin = unit(c(4, 4, 4, 0), "pt"))
p_catch_proj_rnd <- df_sel$projection %>%
  filter(fhist == "random") %>%
  ggplot(aes(x = year - 100, y = catch, colour = idx_sel)) +
  geom_line(size = 0.3, show.legend = FALSE) +
  scale_x_continuous(breaks = c(0, 5, 10)) + 
  scale_colour_manual(values = cols_sel) + 
  coord_cartesian(ylim = c(0, 2.05), xlim = c(0, 10.5),
                  expand = FALSE) +
  labs(x = "year",
       y = expression(Catch/MSY)) + 
  theme_bw(base_size = 8) +
  theme(plot.margin = unit(c(0, 0.1, 4, 4), "pt"))
p_catch_smry_rnd <- df_sel$summary %>%
  filter(fhist == "random") %>%
  mutate(x = 1) %>%
  ggplot(aes(y = catch, x = x, fill = idx_sel, group = idx_sel)) +
  geom_boxplot(outlier.size = 0.3, size = 0.3,
               width = 1, position = position_dodge(1),
               show.legend = FALSE) +
  scale_fill_manual(values = cols_sel) + 
  coord_cartesian(ylim = c(0, 2.05), xlim = c(0.4, 1.6),  
                  expand = FALSE) +
  theme_bw(base_size = 8) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.margin = unit(c(0, 4, 4, 0), "pt"))
p_rnd <- p_ssb_proj_rnd + p_ssb_smry_rnd + 
  p_catch_proj_rnd + p_catch_smry_rnd + 
  plot_layout(ncol = 2, widths = c(1, 0.5))

### combine all fishing histories
p <- wrap_plots(
  p_ssb_proj_ow + labs(title = "(a) one-way"), p_ssb_smry_ow,
  p_ssb_proj_rc + labs(title = "(b) roller-coaster"), p_ssb_smry_rc,
  p_catch_proj_ow, p_catch_smry_ow,
  p_catch_proj_rc, p_catch_smry_rc,
  p_ssb_proj_rnd + labs(title = "(c) random"), p_ssb_smry_rnd, 
  plot_spacer(), plot_spacer(),
  p_catch_proj_rnd, p_catch_smry_rnd,
  plot_spacer(), plot_spacer(),
  ncol = 4, widths = c(1, 0.5, 1, 0.5)
)
ggsave(filename = "output/plots/paper/pol_sensitivity_idx_sel_fhist.png", 
       type = "cairo", plot = p,
       width = 16, height = 12, units = "cm", dpi = 600)
ggsave(filename = "output/plots/paper/pol_sensitivity_idx_sel_fhist.pdf", 
       plot = p, width = 16, height = 12, units = "cm")


### plot selectivity
idx_sel <- foreach(idx_sel = c("tsb", "ssb", "catch", "dome_shaped"), 
          .combine = rbind) %do% {
  #browser()
  res <- readRDS(paste0("output/hr/500_50/sensitivity_idx_sel/", "one-way", "/",
                        "pol", "/mp_length_1_TRUE_1_1_1_Inf_0_0.2_0.2_0_0_",
                        "0.6_0_0.75", 
                        ifelse(identical(idx_sel, "tsb"), "", 
                               paste0("_", idx_sel)),
                        ".rds"))
  sel <- iterMedians(res@oem@observations$idx$sel[, 1])
  if (identical(idx_sel, "tsb")) sel[] <- 1
  
  res <- data.frame(stock = "pol", idx_sel = idx_sel, 
                    age = as.numeric(dimnames(sel)$age),
                    selectivity = c(sel),
                    stringsAsFactors = FALSE)
  return(res)
}
saveRDS(idx_sel, file = "output/hr_sensitivity_sel_curves.rds")
idx_sel <- readRDS("output/hr_sensitivity_sel_curves.rds")



idx_sel <- idx_sel %>%
  mutate(idx_sel = factor(idx_sel, 
                          levels = c("ssb", "dome_shaped", "catch", "tsb"),
                          labels = c("maturity", "dome shaped", "commercial", 
                                     "total biomass")))
p <- idx_sel %>%
  ggplot(aes(x = age, y = selectivity, colour = idx_sel)) +
  geom_line(size = 0.3) +
  scale_colour_manual(name = "", values = cols_sel) + 
  labs(x = "Age [years]", y = "Index selectivity") +
  scale_x_continuous(limits = c(0, 16)) + 
  theme_bw(base_size = 8) + 
  theme(legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.6, "lines"),
        legend.background = element_blank(),
        legend.key = element_blank())
p
ggsave(filename = "output/plots/paper/pol_sensitivity_idx_sel_age.png", 
       type = "cairo", plot = p,
       width = 8.5, height = 5, units = "cm", dpi = 600)
ggsave(filename = "output/plots/paper/pol_sensitivity_idx_sel_age.pdf", plot = p,
       width = 8.5, height = 5, units = "cm")


### ------------------------------------------------------------------------ ###
### pure harvest rate - catch ~ depletion & HR ####
### ------------------------------------------------------------------------ ###

stocks_subset <- stocks$stock#"pol"
### extract and "correct" time series (run on HPC)
# res <- foreach(stock = stocks_subset, .errorhandling = "pass") %do% {
#   #browser()
#     tmp <- readRDS(paste0("output/hr/10000_100/concept/random/", stock, 
#                           "/mp_uniform_1_FALSE_-1_1_1_Inf_0_0.2_0_0.6_0.rds"))
#     ### stock metrics
#     SSBs <- FLCore::ssb(tmp@stock)
#     Fs <- FLCore::fbar(tmp@stock)
#     Cs <- FLCore::catch(tmp@stock)
#     yrs_check <- ac(100:200)
#     yrs <- dim(tmp@stock[, yrs_check])[2]
#     its <- dim(tmp@stock[, yrs_check])[6]
#     ### collapse correction
#     if (isTRUE(TRUE)) {
#       ### find collapses
#       cd <- sapply(seq(its), function(x) {
#         min_yr <- min(which(SSBs[, yrs_check,,,, x] < 1))
#         if (is.finite(min_yr)) {
#           all_yrs <- min_yr:yrs
#         } else {
#           all_yrs <- NA
#         }
#         all_yrs + (x - 1)*yrs
#       })
#       cd <- unlist(cd)
#       cd <- cd[which(!is.na(cd))]
#       ### remove values
#       SSBs[, yrs_check]@.Data[cd] <- 0
#       Cs[, yrs_check]@.Data[cd] <- 0
#       Fs[, yrs_check]@.Data[cd] <- 0
#     }
#     brp <- brps[[stock]]
#     Bmsy <- c(refpts(brp)["msy", "ssb"])
#     Fmsy <- c(refpts(brp)["msy", "harvest"])
#     Cmsy <- c(refpts(brp)["msy", "yield"])
#     ### summarise
#     qnts <- FLQuants(SSB = SSBs/Bmsy, F = Fs/Fmsy, Catch = Cs/Cmsy)
#     rm(stk); gc()
#     return(qnts)
# }
# names(res) <- stocks_subset
# saveRDS(res, paste0("output/hr/10000_100/concept/random/smry_corrected.rds"))
# 
# res <- readRDS("output/hr/10000_100/concept/random/smry_corrected.rds")
# stocks_subset <- stocks$stock
# ### harvest rates
# set.seed(33)
# hr_rates <- runif(n = 10000, min = 0, max = 1)
# stats_yrs <- c(10, 50, 100)
# 
# ### uniform: extract stats
# stats <- foreach(stock = stocks_subset, .combine = bind_rows) %:%
#   foreach(years = stats_yrs, .combine = bind_rows) %do% {#browser()
#     
#     ### extract stats (median per iteration)
#     # SSB_i <- c(apply(res[[stock]]$SSB[, ac(seq(from = 101, to = 100 + years))], 
#     #                  6, median))
#     Catch_i <- c(apply(res[[stock]]$Catch[, ac(seq(from = 101, to = 100 + years))], 
#                      6, median))
#     HR_i <- hr_rates
#     Status_i <- c(res[[stock]]$SSB[, ac(100)])
#     
#     ### define steps for summary & interpolation
#     Status_i_steps <- seq(from = 0.0, to = max(Status_i), by = 0.05)
#     HR_i_steps <- seq(from = 0.025, to = 0.975, by = 0.05)
#     
#     ### interpolate 
#     ### Cubic interpolation & extrapolate missing values around edges
#     int_i <- interp(x = Status_i, y = HR_i, z = Catch_i,
#                     xo = Status_i_steps, yo = HR_i_steps,
#                     linear = TRUE, extrap = TRUE, duplicate = "median")
#     
#     ### format for plotting
#     Catch_int <- data.frame(int_i$z)
#     Catch_int[Catch_int < 0] <- 0
#     colnames(Catch_int) <- HR_i_steps
#     Catch_int$Status <- Status_i_steps
#     Catch_int <- Catch_int %>%
#       pivot_longer(-Status, values_to = "catch", names_to = "HR") %>%
#       mutate(HR = as.numeric(as.character(HR))) %>%
#       #filter(!is.na(catch)) %>%
#       mutate(years = years, stock = stock)
#     #str(Catch_int)
#     # Catch_int %>%
#     #   ggplot(aes(x = Status, y = HR, fill = catch)) +
#     #   geom_raster(interpolate = FALSE) +
#     #   scale_fill_gradientn(expression(catch/MSY),
#     #                        colours = c("red", "orange", "yellow", "green",
#     #                                    "darkgreen"),
#     #                        values = c(0, 0.33, 0.66, 1, max_y)/max_y
#     #   ) +
#     #   theme_bw(base_size = 7) +
#     #   labs(x = expression(initial~stock~status~(SSB/italic(B)[MSY])),
#     #        y = "harvest rate")
#     return(Catch_int)
# }
# 
# saveRDS(stats, file = "output/hr_pure_catch_vs_hr-depl-period.rds")
# write.csv(stats, file = "output/hr_pure_catch_vs_hr-depl-period.csv", 
#           row.names = FALSE)
stats <- readRDS("output/hr_pure_catch_vs_hr-depl-period.rds")

### add von Bertalanffy k
stats2 <- stats %>% 
  left_join(stocks[, c("stock", "k")]) %>%
  mutate(stock_k = paste0(stock, "~~(italic(k)==", k, ")"),
         years_label = factor(years, levels = c(10, 50, 100),
                              labels = paste0(c(10, 50, 100), "~years"))) %>%
  mutate(stock_k = factor(stock_k, levels = unique(stock_k))) %>%
  filter(!is.na(catch)) %>%
  filter(Status <= 3)
max_y_global <- max(stats2$catch, na.rm = TRUE)

### plot (in groups)
groups <- append(list(c(1, 12, 18, 22, 26, 29)), 
                 split(seq_along(stocks$stock), rep(1:5, each = 6)))
for (i in seq_along(groups)) {
  stats3 <- stats2 %>%
    filter(stock %in% stocks$stock[groups[[i]]])
  if (identical(i, 1)) {stats3 <- stats3 %>% filter(years %in% c(10, 50))}
  max_y_local <- max(stats3$catch, na.rm = TRUE)
  stats3 %>%
    ggplot(aes(x = Status, y = HR, fill = catch)) +
    facet_grid(years_label ~ stock_k, labeller = label_parsed) +
    geom_raster(interpolate = FALSE) +
    # scale_fill_gradientn(expression(catch/MSY),
    #                      colours = c("red", "orange", "yellow", "green",
    #                                  "darkgreen"),
    #                      values = c(0, 0.33, 0.66, 1, max_y_global)/max_y_local
    scale_fill_gradientn(expression(catch/MSY),
                         colours = c("brown", "white", "darkblue"),
                         values = c(0, 1, max_y_global)/max_y_local
    ) +
    # geom_contour(aes(z = catch), breaks = c(0, 1, max_y_global)/max_y_local,
    #              colour = "black", alpha = 0.5, size = 0.1) +
    theme_bw(base_size = 7) +
    labs(x = expression(initial~stock~status~(SSB[y==0]/B[MSY])),
         y = "harvest rate") +
    xlim(c(0, 3))
  ggsave(filename = paste0("output/plots/paper/HR_principle_catch_examples",
                           names(groups)[i], ".png"), 
         type = "cairo", width = 17, height = ifelse(i < 2, 6, 8), 
         units = "cm", dpi = 600)
  ggsave(filename = paste0("output/plots/paper/HR_principle_catch_examples", 
                           names(groups)[i], ".pdf"), 
         width = 17, height = ifelse(i < 2, 6, 8), 
         units = "cm")
}
### figure with examples for paper
### without grid lines
stats3 <- stats2 %>%
  filter(stock %in% stocks$stock[groups[[1]]]) %>% 
  filter(years %in% c(10, 50))
max_y_local <- max(stats3$catch, na.rm = TRUE)
p <- stats3 %>%
  ggplot(aes(x = Status, y = HR, fill = catch)) +
  facet_grid(years_label ~ stock_k, labeller = label_parsed) +
  geom_raster(interpolate = FALSE) +
  scale_fill_gradientn(expression(catch/MSY),
                       colours = c("brown", "white", "darkblue"),
                       values = c(0, 1, max_y_global)/max_y_local
  ) +
  labs(x = expression(initial~stock~status~(SSB[y==0]/B[MSY])),
       y = "harvest rate") +
  xlim(c(0, 3)) +
  theme_bw(base_size = 7) + 
  theme(panel.grid.minor = element_blank())
### save as postscript file, then convert to PDF
cairo_ps(filename = "output/plots/paper/HR_principle_catch_examples_ps.ps", 
         width = 17/2.54, height = 6/2.54, fallback_resolution = 600)
print(p)
dev.off()

### ------------------------------------------------------------------------ ###
### HR (length target) - multipliers - 50 years ####
### ------------------------------------------------------------------------ ###

### load results
res_def <- readRDS("output/hr_all_mult_runs_stats_nocaps.rds")
### add stock labels
res_def <- res_def %>%
  left_join(stocks[, c("stock", "k")]) %>%
  mutate(stock_k = paste0(stock, "~(italic(k)==", k, ")")) %>%
  mutate(stock_k = factor(stock_k, levels = unique(stock_k))) %>%
  mutate(stock = factor(stock, levels = stocks$stock))
### sort fishing histories
res_def <- res_def %>%
  mutate(fhist = factor(fhist, 
                        levels = c("one-way", "roller-coaster", "random")))

res_def_colours <- c("one-way" = brewer.pal(n = 4, name = "Set1")[1], 
                        "roller-coaster" = brewer.pal(n = 4, name = "Set1")[4], 
                        "random" = brewer.pal(n = 4, name = "Set1")[2])
res_def_linetype <- c("one-way" = "solid", 
                        "roller-coaster" = "dotted", 
                        "random" = "dashed")

### plot all stocks
res_def_p <- res_def %>%
  #filter(stat_yrs == "all") %>%
  select(multiplier, risk_Blim, SSB_rel, Fbar_rel, Catch_rel, ICV,
         stock, fhist, stock_k) %>%
  mutate(ICV = ifelse(multiplier == 0, NA, ICV)) %>%
  pivot_longer(c(SSB_rel, Fbar_rel, Catch_rel, risk_Blim, ICV), 
               names_to = "key", values_to = "value") %>%
  mutate(stat = factor(key, levels = c("SSB_rel", "Fbar_rel", "Catch_rel",
                                       "risk_Blim", "ICV", "fitness"), 
                       labels = c("SSB/B[MSY]", "F/F[MSY]", "Catch/MSY", 
                                  "B[lim]~risk", "ICV", "fitness~value")))
stats_targets <- data.frame(stat = c("SSB/B[MSY]", "F/F[MSY]", "Catch/MSY", 
                                     "B[lim]~risk", "ICV"),
                            target = c(1, 1, 1, 0, 0))
# p <- res_def_p %>% 
#   ggplot(aes(x = multiplier, y = value,
#              colour = as.factor(fhist), linetype = as.factor(fhist))) +
#   geom_line(size = 0.3) +
#   # geom_hline(data = data.frame(stat = "B[lim]~risk", y = 0.05),
#   #            aes(yintercept = y), colour = "red") +
#   facet_grid(stat ~ stock_k, labeller = "label_parsed", switch = "y",
#              scales = "free_y") +
#   scale_linetype_discrete("fishing\nhistory") +
#   scale_colour_discrete("fishing\nhistory") +
#   theme_bw(base_size = 8) +
#   theme(strip.placement.y = "outside",
#         strip.background.y = element_blank(),
#         strip.text.y = element_text(size = 8),
#         strip.text.x = element_text(size = 6)) +
#   labs(x = "multiplier", y = "") +
#   ylim(c(0, NA))# +
# # scale_x_continuous(breaks = c(0, 0.5, 1)#,
# #                    #expand = expansion(mult = c(0.1, 0.1))
# #                    )
# p
# ggsave(filename = "output/plots/length_all_def_mult.pdf",
#        width = 50, height = 10, units = "cm")
# ggsave(filename = "output/plots/length_all_def_mult.png", type = "cairo",
#        width = 50, height = 10, units = "cm", dpi = 600)

### find max catch
max_catch <- res_def %>%
  group_by(stock, fhist, k) %>%
  filter(Catch_rel == max(Catch_rel))
summary(max_catch)
table(max_catch$multiplier)

### plot subset and stats individually
p_SSB <- res_def_p %>% 
  filter(stock %in% c("ang3", "pol", "bll", "san")) %>%
  filter(key == "SSB_rel") %>%
  ggplot(aes(x = multiplier, y = value,
             colour = as.factor(fhist), linetype = as.factor(fhist))) +
  geom_vline(data = max_catch %>% 
               filter(stock %in% c("ang3", "pol", "bll", "san")),
             aes(xintercept = multiplier, linetype = fhist),
             size = 0.15, colour = "grey40", show.legend = FALSE) +
  geom_line(size = 0.3) +
  facet_grid(stat ~ stock_k, labeller = "label_parsed", switch = "y",
             scales = "free_y") +
  scale_linetype_manual("fishing history", values = res_def_linetype) +
  scale_colour_manual("fishing history", values = res_def_colours) +
  labs(x = "multiplier", y = "") +
  ylim(c(0, NA)) +
  theme_bw(base_size = 8) +
  theme(strip.placement.y = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 8),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.6, "lines"),
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.position = c(0.9, 0.75),
        axis.title.y = element_blank(),
        plot.margin = unit(x = c(1, 3, 0, 3), units = "pt"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())
p_SSB
p_Catch <- res_def_p %>% 
  filter(stock %in% c("ang3", "pol", "bll", "san")) %>%
  filter(key == "Catch_rel") %>%
  ggplot(aes(x = multiplier, y = value,
             colour = as.factor(fhist), linetype = as.factor(fhist))) +
  geom_vline(data = max_catch %>% 
               filter(stock %in% c("ang3", "pol", "bll", "san")),
             aes(xintercept = multiplier, 
                 linetype = as.factor(fhist)),
             size = 0.15, colour = "grey40", show.legend = FALSE) +
  geom_line(size = 0.3, show.legend = FALSE) +
  geom_point(data = res_def_p %>% 
               filter(stock %in% c("ang3", "pol", "bll", "san") & 
                        key == "Catch_rel") %>%
               group_by(stock, fhist) %>%
               filter(value == max(value)),
               aes(shape = as.factor(fhist)), show.legend = FALSE,
               size = 0.7) +
  scale_linetype_manual("fishing history", values = res_def_linetype) +
  scale_colour_manual("fishing history", values = res_def_colours) +
  facet_grid(stat ~ stock_k, labeller = "label_parsed", switch = "y",
             scales = "free_y") +
  labs(x = "multiplier", y = "") +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0, NA)) +
  theme_bw(base_size = 8) +
  theme(strip.placement.y = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        strip.text.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = unit(x = c(1, 3, 0, 3), units = "pt"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())
p_Catch
p_risk <- res_def_p %>% 
  filter(stock %in% c("ang3", "pol", "bll", "san")) %>%
  filter(key == "risk_Blim") %>%
  ggplot(aes(x = multiplier, y = value,
             colour = as.factor(fhist), linetype = as.factor(fhist))) +
  geom_vline(data = max_catch %>% 
               filter(stock %in% c("ang3", "pol", "bll", "san")),
             aes(xintercept = multiplier, 
                 linetype = as.factor(fhist)),
             size = 0.15, colour = "grey40", show.legend = FALSE) +
  geom_line(size = 0.3, show.legend = FALSE) +
  scale_linetype_manual("fishing history", values = res_def_linetype) +
  scale_colour_manual("fishing history", values = res_def_colours) +
  facet_grid(stat ~ stock_k, labeller = "label_parsed", switch = "y",
             scales = "free_y") +
  labs(x = "multiplier", y = "") +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1)) +
  theme_bw(base_size = 8) +
  theme(strip.placement.y = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        strip.text.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = unit(x = c(0, 3, 3, 3), units = "pt"))
p_risk
p_stats <- plot_grid(p_SSB, p_Catch, p_risk,
                     ncol = 1, align = "v",
                     rel_heights = c(1.25, 1, 1.25))
p_stats


### estimate correlation
max_catch_cor <- foreach(fhist = c("one-way", "roller-coaster", "random")) %:%
  foreach(stat = c("Catch_rel", "multiplier")) %do% {
    #browser()
    data_i <- max_catch %>%
      filter(fhist == !!fhist) %>%
      select(stock, k, !!stat) %>%
      rename(stat = !!stat)
    cor_res <- cor.test(y = data_i$stat, x = data_i$k)
    cor_res <- capture.output(cor_res)
    c(paste0("\n\nCorrelation for fhist: ", fhist, "; stat: ", stat, "\n"),
      paste(cor_res, "\n"))
}
lapply(do.call(c, max_catch_cor), cat)
# Correlation for fhist: one-way; stat: Catch_rel
# 
# Pearson's product-moment correlation 
#   
#  data:  data_i$k and data_i$stat 
#  t = -14.189, df = 27, p-value = 4.903e-14 
#  alternative hypothesis: true correlation is not equal to 0 
#  95 percent confidence interval: 
#   -0.9712578 -0.8729319 
#  sample estimates: 
#         cor  
#  -0.9390145  
#   
# 
# 
# Correlation for fhist: one-way; stat: multiplier
#   
#  	Pearson's product-moment correlation 
# 
# data:  data_i$k and data_i$stat 
# t = -10.358, df = 27, p-value = 6.636e-11 
# alternative hypothesis: true correlation is not equal to 0 
# 95 percent confidence interval: 
#   -0.9493434 -0.7842543 
# sample estimates: 
#   cor  
# -0.8938402  
# 
# 
# 
# Correlation for fhist: roller-coaster; stat: Catch_rel
# 
# Pearson's product-moment correlation 
#   
#  data:  data_i$k and data_i$stat 
#  t = -14.109, df = 27, p-value = 5.613e-14 
#  alternative hypothesis: true correlation is not equal to 0 
#  95 percent confidence interval: 
#   -0.9709556 -0.8716632 
#  sample estimates: 
#         cor  
#  -0.9383839  
#   
# 
# 
# Correlation for fhist: roller-coaster; stat: multiplier
#   
#  	Pearson's product-moment correlation 
# 
# data:  data_i$k and data_i$stat 
# t = -11.025, df = 27, p-value = 1.682e-11 
# alternative hypothesis: true correlation is not equal to 0 
# 95 percent confidence interval: 
#   -0.9545954 -0.8049094 
# sample estimates: 
#   cor  
# -0.9045649  
# 
# 
# 
# Correlation for fhist: random; stat: Catch_rel
# 
# Pearson's product-moment correlation 
#   
#  data:  data_i$k and data_i$stat 
#  t = -10.302, df = 27, p-value = 7.47e-11 
#  alternative hypothesis: true correlation is not equal to 0 
#  95 percent confidence interval: 
#   -0.9488603 -0.7823728 
#  sample estimates: 
#        cor  
#  -0.892857  
#   
# 
# 
# Correlation for fhist: random; stat: multiplier
#   
#  	Pearson's product-moment correlation 
# 
# data:  data_i$k and data_i$stat 
# t = -8.5716, df = 27, p-value = 3.473e-09 
# alternative hypothesis: true correlation is not equal to 0 
# 95 percent confidence interval: 
#   -0.9301314 -0.7116912 
# sample estimates: 
#   cor  
# -0.8551425  



### prepare for plotting
max_catch_p <- max_catch %>%
  select(multiplier, Catch_rel, risk_Blim, stock, fhist, k) %>%
  pivot_longer(c(multiplier, Catch_rel, risk_Blim), 
               names_to = "key", values_to = "value") %>%
  mutate(stat = factor(key, levels = c("multiplier", "Catch_rel", "risk_Blim"), 
                       labels = c("multiplier~(x)", "Catch/MSY", "B[lim]~risk")))
### plot
p_max_catch <- max_catch_p %>%
  filter(stat %in% c("multiplier~(x)", "Catch/MSY")) %>%
  ggplot(aes(x = k, y = value, colour = fhist, shape = fhist, 
             linetype = fhist)) +
  geom_smooth(method = lm, se = FALSE, size = 0.3) +
  geom_point(size = 0.5) +
  scale_linetype_manual("fishing history", values = res_def_linetype) +
  scale_colour_manual("fishing history", values = res_def_colours) +
  scale_shape("fishing history") +
  facet_wrap(~ stat, scales = "free_y", strip.position = "left",
             labeller = "label_parsed", ncol = 1) +
  #ylim(c(0, NA)) +  xlim(c(0, NA)) +
  coord_cartesian(ylim = c(0, NA), xlim = c(0, 1)) + 
  labs(x = expression(k~"[year"^{-1}*"]"), y = "") +
  theme_bw(base_size = 8) +
  theme(strip.placement.y = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8), 
        axis.title.y = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.6, "lines"),
        legend.position = c(0.77, 0.92))
p_max_catch

### combine plots
p <- plot_grid(p_stats, p_max_catch,
               ncol = 2, rel_widths = c(0.67, 0.33),
               labels = c("(a)", "(b)"), label_size = 9, align = "h", 
               axis = "t")
p
ggsave(filename = "output/plots/paper/HR_mult_max_catch_cor.png", 
       type = "cairo", plot = p,
       width = 17, height = 10, units = "cm", dpi = 600)
ggsave(filename = "output/plots/paper/HR_mult_max_catch_cor.pdf", plot = p,
       width = 17, height = 10, units = "cm")

### ------------------------------------------------------------------------ ###
### GA: pollack explorations with HR components ####
### ------------------------------------------------------------------------ ###
pol_GA <- readRDS("output/hr_pol_comps_stats.rds")
# pol_GA %>% View()

### function for transforming y-axis origin in plots
trans_from <- function(from = 1) {
  trans <- function(x) x - from
  inv <- function(x) x + from
  trans_new("from", trans, inv, 
            domain = c(from, Inf))
}
### plot
p_PA_SSB <- pol_GA %>% 
  filter(objective == "MSY-PA") %>%
  ggplot(aes(x = parameters, y = SSB_rel)) +
  geom_hline(yintercept = 1, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = position_dodge2(preserve = "single"), width = 0.8, 
           show.legend = FALSE, colour = "black", size = 0.1) +
  facet_grid("SSB/B[MSY]" ~ fhist, scales = "free", space = "free_x", 
             switch = "y", labeller = "label_parsed") +
  labs(y = "", x = "") +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(-0.01, units = "cm"),
        strip.placement.y = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        plot.margin = unit(x = c(1, 3, 0, 3), units = "pt"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_y_continuous(trans = trans_from(), limits = c(0, 2))
p_MSY_SSB <- pol_GA %>% 
  filter(objective == "MSY") %>%
  ggplot(aes(x = parameters, y = SSB_rel)) +
  geom_hline(yintercept = 1, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = position_dodge2(preserve = "single"), width = 0.8, 
           show.legend = FALSE, colour = "black", size = 0.1) +
  facet_grid("SSB/B[MSY]" ~ fhist, scales = "free", space = "free_x", 
             switch = "y", labeller = "label_parsed") +
  labs(y = "", x = "") +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(-0.01, units = "cm"),
        strip.placement.y = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        plot.margin = unit(x = c(1, 3, 0, 3), units = "pt"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_y_continuous(trans = trans_from(), limits = c(0, 2))
p_PA_catch <- pol_GA %>% 
  filter(objective == "MSY-PA") %>%
  ggplot(aes(x = parameters, y = Catch_rel)) +
  geom_hline(yintercept = 1, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = position_dodge2(preserve = "single"), width = 0.8, 
           show.legend = FALSE, colour = "black", size = 0.1) +
  facet_grid("Catch/MSY" ~ fhist, scales = "free", space = "free_x", 
             switch = "y", labeller = "label_parsed") +
  labs(y = "", x = "") +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(-0.01, units = "cm"),
        strip.text.x = element_blank(),
        strip.placement.y = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        plot.margin = unit(x = c(0, 3, 0, 3), units = "pt"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_y_continuous(trans = trans_from(), limits = c(0, 2))
p_MSY_catch <- pol_GA %>% 
  filter(objective == "MSY") %>%
  ggplot(aes(x = parameters, y = Catch_rel)) +
  geom_hline(yintercept = 1, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = position_dodge2(preserve = "single"), width = 0.8, 
           show.legend = FALSE, colour = "black", size = 0.1) +
  facet_grid("Catch/MSY" ~ fhist, scales = "free", space = "free_x", 
             switch = "y", labeller = "label_parsed") +
  labs(y = "", x = "") +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(-0.01, units = "cm"),
        strip.text.x = element_blank(),
        strip.placement.y = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        plot.margin = unit(x = c(0, 3, 0, 3), units = "pt"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_y_continuous(trans = trans_from(), limits = c(0, 2))
p_PA_ICV <- pol_GA %>% 
  filter(objective == "MSY-PA") %>%
  ggplot(aes(x = parameters, y = ICV)) +
  geom_hline(yintercept = 0, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = position_dodge2(preserve = "single"), width = 0.8, 
           show.legend = FALSE, colour = "black", size = 0.1) +
  facet_grid("ICV" ~ fhist, scales = "free", space = "free_x", 
             switch = "y", labeller = "label_parsed") +
  labs(y = "", x = "") +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(-0.01, units = "cm"),
        strip.text.x = element_blank(),
        strip.placement.y = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        plot.margin = unit(x = c(0, 3, 0, 3), units = "pt"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_y_continuous(limits = c(0, 0.5))
p_MSY_ICV <- pol_GA %>% 
  filter(objective == "MSY") %>%
  ggplot(aes(x = parameters, y = ICV)) +
  geom_hline(yintercept = 0, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = position_dodge2(preserve = "single"), width = 0.8, 
           show.legend = FALSE, colour = "black", size = 0.1) +
  facet_grid("ICV" ~ fhist, scales = "free", space = "free_x", 
             switch = "y", labeller = "label_parsed") +
  labs(y = "", x = "") +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(-0.01, units = "cm"),
        strip.text.x = element_blank(),
        strip.placement.y = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        plot.margin = unit(x = c(0, 3, 0, 3), units = "pt"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_y_continuous(limits = c(0, 0.5))
p_PA_risk <- pol_GA %>% 
  filter(objective == "MSY-PA") %>%
  ggplot(aes(x = parameters, y = risk_Blim)) +
  geom_hline(yintercept = 0, linetype = "solid", size = 0.5, colour = "grey") +
  geom_hline(yintercept = 0.055, linetype = "solid", size = 0.5, 
             colour = "red") +
  geom_col(position = position_dodge2(preserve = "single"), width = 0.8, 
           show.legend = FALSE, colour = "black", size = 0.1) +
  facet_grid("B[lim]~risk" ~ fhist, scales = "free", space = "free_x", 
             switch = "y", labeller = "label_parsed") +
  labs(y = "", x = "") +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(-0.01, units = "cm"),
        strip.text.x = element_blank(),
        strip.placement.y = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        plot.margin = unit(x = c(0, 3, 0, 3), units = "pt"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_y_continuous(limits = c(0, 0.5))
p_MSY_risk <- pol_GA %>% 
  filter(objective == "MSY") %>%
  ggplot(aes(x = parameters, y = risk_Blim)) +
  geom_hline(yintercept = 0, linetype = "solid", size = 0.5, colour = "grey") +
  # geom_hline(yintercept = 0.055, linetype = "solid", size = 0.5, 
  #            colour = "red") +
  geom_col(position = position_dodge2(preserve = "single"), width = 0.8, 
           show.legend = FALSE, colour = "black", size = 0.1) +
  facet_grid("B[lim]~risk" ~ fhist, scales = "free", space = "free_x", 
             switch = "y", labeller = "label_parsed") +
  labs(y = "", x = "") +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(-0.01, units = "cm"),
        strip.text.x = element_blank(),
        strip.placement.y = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        plot.margin = unit(x = c(0, 3, 0, 3), units = "pt"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_y_continuous(limits = c(0, 0.5))
p_PA_fitness <- pol_GA %>% 
  filter(objective == "MSY-PA") %>%
  mutate(fitness_colour = ifelse(risk_Blim >= 0.055, NA, fitness)) %>%
  ggplot(aes(x = parameters, y = fitness)) +
  geom_hline(yintercept = 0, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(aes(fill = fitness_colour, 
               alpha = ""), ### dummy to separate legend
           position = "dodge", show.legend = TRUE, width = 0.8, 
           colour = "black", 
           size = 0.1) +
  scale_fill_gradient("fitness",
                      high = "white", low = "#2222ff", na.value = "#ff6666",
                      limits = c(-6.5, 0), 
                      breaks = c(-6, -3, 0), 
                      labels = c("-6 (worst)", -3, "0 (best)")) +
  scale_alpha_manual(values = 1, labels = expression("risk" > 5*"%"),
                     guide = guide_legend(order = 2)) +
  guides(alpha = guide_legend("", 
                              override.aes = list(colour = "#ff6666",
                                                  fill = "#ff6666"), 
                              order = 2)) +
  facet_grid("fitness" ~ fhist, scales = "free", space = "free_x", 
             switch = "y", labeller = "label_parsed") +
  labs(y = "", x = "") +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(-0.01, units = "cm"),
        strip.text.x = element_blank(),
        strip.placement.y = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5,
                                   lineheight = 0.7),
        plot.margin = unit(x = c(0, 3, 3, 3), units = "pt"),
        axis.title = element_blank(),
        legend.key.width = unit(0.7, "lines"),
        legend.key.height = unit(0.5, "lines"),
        legend.title = element_text(size = 7),
        legend.justification = 1) +
  scale_y_continuous(limits = c(-6, 0), breaks = c(0, -3, -6))
p_MSY_fitness <- pol_GA %>% 
  filter(objective == "MSY") %>%
  mutate(fitness_colour = fitness) %>%
  ggplot(aes(x = parameters, y = fitness)) +
  geom_hline(yintercept = 0, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(aes(fill = fitness_colour), 
           position = "dodge", show.legend = FALSE, width = 0.8, 
           colour = "black", 
           size = 0.1) +
  scale_fill_gradient("fitness",
                      high = "white", low = "#2222ff", na.value = "#ff6666",
                      limits = c(-6.5, 0), 
                      breaks = c(-6, -3, 0), 
                      labels = c("-6 (worst)", -3, "0 (best)")) +
  guides(alpha = guide_legend("", 
                              override.aes = list(colour = "#ff6666",
                                                  fill = "#ff6666"), 
                              order = 2)) +
  facet_grid("fitness" ~ fhist, scales = "free", space = "free_x", 
             switch = "y", labeller = "label_parsed") +
  labs(y = "", x = "") +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(-0.01, units = "cm"),
        strip.text.x = element_blank(),
        strip.placement.y = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5,
                                   lineheight = 0.7),
        plot.margin = unit(x = c(0, 3, 3, 3), units = "pt"),
        axis.title = element_blank(),
        legend.key.width = unit(0.5, "lines"),
        legend.key.height = unit(0.5, "lines"),
        legend.title = element_text(size = 7), 
        legend.spacing.y = unit(-20, "pt")) +
  scale_y_continuous(limits = c(-6, 0), breaks = c(0, -3, -6))

p <- plot_grid(
  plot_grid(ggplot() + theme_nothing(), 
            ggplot() + theme_nothing(), 
            ggplot() + theme_nothing(),
            labels = c("(a) fitness function: MSY", 
                       "(b) fitness function: MSY-PA", ""),
            label_size = 9, label_x = c(-0.2, -0.35),
            rel_widths = c(1, 1, 0.2)),
  plot_grid(plot_grid(p_MSY_SSB, p_MSY_catch, p_MSY_ICV, p_MSY_risk,
                      p_MSY_fitness,
                      ncol = 1, align = "v", axis = "lr", 
                      rel_heights = c(1.25, 1, 1, 1, 2.1)),
            plot_grid(p_PA_SSB, p_PA_catch, p_PA_ICV, p_PA_risk, 
                      p_PA_fitness + theme(legend.position = "none"),
                      ncol = 1, align = "v", axis = "lr", 
                      rel_heights = c(1.25, 1, 1, 1, 2.1)),
            plot_grid(NULL, NULL, NULL, NULL, 
                      plot_grid(get_legend(p_PA_fitness), 
                                ggplot() + theme_nothing(), ncol = 2,
                                rel_widths = c(1, 2)),
                      rel_heights = c(1, 1, 1, 1, 2)),
            ncol = 3, rel_widths = c(1, 1, 0.3)),
  ncol = 1, rel_heights = c(0.2, 6))

ggsave(filename = "output/plots/paper/hr_GA_pol_comps.png", type = "cairo",
       width = 17, height = 12, units = "cm", dpi = 600, plot = p)
ggsave(filename = "output/plots/paper/hr_GA_pol_comps.pdf", plot = p,
       width = 17, height = 12, units = "cm")




p_MSY <- pol_GA %>% 
  filter(objective == "MSY") %>%
  mutate(parameters = factor(parameters, levels = rev(levels(parameters)))) %>%
  mutate(fitness_colour = ifelse(risk_Blim >= 0.055 & objective == "MSY-PA", 
                                 NA, fitness)) %>%
  ggplot(aes(x = parameters, y = fitness)) +
  geom_hline(yintercept = 0, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(aes(fill = fitness_colour, 
               alpha = ""), ### dummy to separate legend
           position = "dodge", show.legend = FALSE, width = 0.8, 
           colour = "black", 
           size = 0.1) +
  geom_blank(
    data = data.frame(fhist = c("one-way", "random", "one-way", "random"),
                      objective = c("MSY-PA", "MSY-PA", "MSY", "MSY"),
                      parameters = "all",
                      fitness = -1)) +
  coord_flip(ylim = c(-1.1, 0)) +
  scale_fill_gradient("fitness",
                      high = "white", low = "#2222ff", na.value = "#ff6666",
                      limits = c(-1.1, 0),
                      breaks = c(-1, -0.5, 0),
                      labels = c("-6 (worst)", -3, "0 (best)")) +
  scale_alpha_manual(values = 1, labels = expression("risk" > 5*"%"),
                     guide = guide_legend(order = 2)) +
  guides(alpha = guide_legend("",
                              override.aes = list(colour = "#ff6666",
                                                  fill = "#ff6666"),
                              order = 2)) +
  facet_grid(~ fhist) +
  scale_y_continuous(breaks = c(-1, -0.5, 0)) +
  labs(x = "optimised parameters", title = "(a) fitness function: MSY") +
  theme_bw(base_size = 8) +
  theme(axis.title.x = element_blank(),
        panel.spacing.x = unit(0, units = "pt"),
        plot.title = element_text(face = "bold", size = 9))
p_PA <- pol_GA %>% 
  filter(objective == "MSY-PA") %>%
  mutate(parameters = factor(parameters, levels = rev(levels(parameters)))) %>%
  mutate(fitness_colour = ifelse(risk_Blim >= 0.055 & objective == "MSY-PA", 
                                 NA, fitness)) %>%
  ggplot(aes(x = parameters, y = fitness)) +
  geom_hline(yintercept = 0, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(aes(fill = fitness_colour, 
               alpha = ""), ### dummy to separate legend
           position = "dodge", show.legend = TRUE, width = 0.8, 
           colour = "black", 
           size = 0.1) +
  geom_blank(
    data = data.frame(fhist = c("one-way", "random", "one-way", "random"),
                      objective = c("MSY-PA", "MSY-PA", "MSY", "MSY"),
                      parameters = "all",
                      fitness = -1)) +
  coord_flip(ylim = c(-1.1, 0)) +
  scale_fill_gradient("fitness",
                      high = "white", low = "#2222ff", na.value = "#ff6666",
                      limits = c(-1.1, 0),
                      breaks = c(-1, -0.5, 0),
                      labels = c("-1 (worst)", -0.5, "0 (best)")) +
  scale_alpha_manual(values = 1, labels = expression("risk" > 5*"%"),
                     guide = guide_legend(order = 2)) +
  guides(alpha = guide_legend("",
                              override.aes = list(colour = "#ff6666",
                                                  fill = "#ff6666"),
                              order = 2)) +
  facet_grid(~ fhist) +
  scale_y_continuous(breaks = c(-1, -0.5, 0)) +
  labs(x = "optimised parameters", title = "(b) fitness function: MSY-PA") +
  theme_bw(base_size = 8) +
  theme(axis.title.x = element_blank(),
        panel.spacing.x = unit(0, units = "pt"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(face = "bold", size = 9),
        legend.key.width = unit(0.5, "lines"),
        legend.key.height = unit(0.5, "lines"))

plot_grid(plot_grid(plot_grid(p_MSY, p_PA + theme(legend.position = "none"),
                              ncol = 2, align = "h", rel_widths = c(1.35, 1)),
                    get_legend(p_PA),
                    ncol = 2, rel_widths = c(1, 0.1)),
          ggplot() + theme_nothing(),
          ncol = 1, rel_heights = c(1, 0.04),
          labels = c("", "fitness"), label_size = 8, label_fontface = "plain",
          label_x = 0.48, label_y = 1.1)
ggsave(filename = "output/plots/paper/hr_GA_pol_comps_fitness.png", 
       type = "cairo",
       width = 17, height = 8, units = "cm", dpi = 600)
ggsave(filename = "output/plots/paper/hr_GA_pol_comps_fitness.pdf",
       width = 17, height = 8, units = "cm")

### ------------------------------------------------------------------------ ###
### all stocks compare HR to rfb - table like plot ####
### ------------------------------------------------------------------------ ###

### stats from 2 over 3 rule
# file.copy(from = "../GA_MSE/output/all_stocks_2over_stats.rds", 
#           to = "output/2over3_all_stocks_stats.rds")
stats_2over3 <- readRDS("output/2over3_all_stocks_stats.rds")

### stats from rfb-rule
# file.copy(from = "../GA_MSE/output/all_stocks_GA_optimised_stats.rds", 
#           to = "output/rfb_all_stocks_GA_optimised_stats.rds")
stats_rfb <- readRDS("output/rfb_all_stocks_GA_optimised_stats.rds")

stats_hr <- readRDS("output/hr_all_GA_stats.rds")

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
    filter(capped == TRUE, optimised == "default" & scenario == "PA") %>%
    mutate(catch_rule = "rfb", target = "none", group = "rfb: default"),
  ### rfb-rule, cond. cap, optimisation with multiplier, target PA
  stats_rfb %>% 
    filter(capped == TRUE, optimised == "mult" & scenario == "PA_capped") %>%
    mutate(catch_rule = "rfb", target = "PA", 
           group = "rfb: mult"),
  ### rfb-rule, cond. cap, optimisation with all, target PA
  stats_rfb %>% 
    filter(capped == TRUE, optimised == "all" & scenario == "PA_capped") %>%
    mutate(catch_rule = "rfb", target = "PA", 
           group = "rfb: all"),
  ### hr, cond. cap, default, target PA
  stats_hr %>%
    filter(objective == "MSY-PA" & scenario == "GA_cond_cap" & 
             parameters == "none") %>%
    mutate(catch_rule = "hr", target = "PA", group = "hr: default"),
  ### hr, cond. cap, optimisation with multiplier, target PA
  stats_hr %>%
    filter(objective == "MSY-PA" & scenario == "GA_cond_cap" & 
             parameters == "mult") %>%
    mutate(catch_rule = "hr", target = "PA", group = "hr: mult"),
  ### hr, cond. cap, optimisation with all parameters, target PA
  stats_hr %>%
    filter(objective == "MSY-PA" & scenario == "GA_cond_cap" & 
             parameters == "all") %>%
    mutate(catch_rule = "hr", target = "PA", group = "hr: all")
)
stats_plot <- stats_plot %>%
  select(fhist, stock, catch_rule, target, group, 
         SSB_rel, Catch_rel, ICV, risk_Blim) %>%
  mutate(penalty = penalty(x = risk_Blim, 
                           negative = FALSE, max = 5, 
                           inflection = 0.06, 
                           steepness = 0.5e+3)) %>%
  mutate(fitness = -abs(SSB_rel - 1) - abs(Catch_rel - 1) - ICV - penalty) %>%
  mutate(fitness_colour = ifelse(risk_Blim >= 0.055, NA, fitness)) %>%
  full_join(stocks %>%
    select(stock, k) %>%
    mutate(stock = factor(stock, levels = stock),
           stock_k = paste0(stock, "~(italic(k)==", sprintf(k, fmt =  "%.2f"), 
                            "*year^-1)")) %>%
              mutate(stock_k = factor(stock_k, levels = rev(stock_k)))) %>%
  mutate(
    group = factor(group, 
                   levels = c("zero-fishing", "2 over 3", 
                              "rfb: default",
                              "rfb: mult", "rfb: all",
                              "hr: default",
                              "hr: mult", "hr: all"),
                   labels = c("(a) zero-fishing", "(b) 2 over 3", 
                              "(c) rfb: default",
                              "(d) rfb: mult", "(e) rfb: all",
                              "(f) hr: default",
                              "(g) hr: mult", "(h) hr: all")))

saveRDS(stats_plot, "output/data_all_comparison_table.rds")
stats_plot <- readRDS("output/data_all_comparison_table.rds")

### plot (table style)
p_all_table <- stats_plot %>%
  ggplot(aes(x = group, y = stock_k, fill = fitness_colour)) +
  geom_raster(aes(alpha = "")) +
  geom_text(aes(label = sprintf(fitness, fmt =  "%.1f")),
            colour = "grey40", size = 2) + 
  scale_fill_gradient("fitness",
                      high = "white", low = "#2222ff", na.value = "#ff6666",
                      limits = c(-8.3, 0), 
                      breaks = c(-8, -4, 0), 
                      labels = c("-8 (worst)", -4, "0 (best)")
                      ) +
  scale_alpha_manual(values = 1, labels = expression("risk" > 5*"%"),
                     guide = guide_legend(order = 2)) +
  guides(alpha = guide_legend("", 
                              override.aes = list(colour = "#ff6666",
                                                  fill = "#ff6666"), 
                              order = 2)) +
  facet_grid(~ fhist) +
  scale_y_discrete(labels = scales::parse_format()) +
  theme_bw(base_size = 8, base_family = "sans") +
  theme(axis.title = element_blank(),
        panel.spacing.x = unit(-0.01, units = "cm"),
        legend.key.width = unit(0.7, "lines"),
        legend.key.height = unit(0.5, "lines"),
        legend.title = element_text(size = 7),
        axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5,
                                   lineheight = 0.7),
        axis.text.y = element_text(hjust = 0),
        panel.grid = element_blank()
  )
p_all_table
ggsave(filename = "output/plots/paper/all_comparison_table.png", 
       plot = p_all_table,
       width = 17, height = 16, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/paper/all_comparison_table.pdf",
       plot = p_all_table,
       width = 17, height = 16, units = "cm", dpi = 600)

### ------------------------------------------------------------------------ ###
### all stocks compare HR to rfb - full details ####
### ------------------------------------------------------------------------ ###

### stats from 2 over 3 rule
stats_2over3 <- readRDS("output/2over3_all_stocks_stats.rds")
### stats from rfb-rule
stats_rfb <- readRDS("output/rfb_all_stocks_GA_optimised_stats.rds")
### harvest rate
stats_hr <- readRDS("output/hr_all_GA_stats.rds")

### combine 
stats_full <- bind_rows(
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
  ### rfb-rule, optimisation with multiplier
  stats_rfb %>% 
    filter(capped == FALSE, optimised == "mult" & scenario == "PA") %>%
    mutate(catch_rule = "rfb", target = "PA", group = "rfb: mult"),
  ### rfb-rule, optimisation with all parameters
  stats_rfb %>% 
    filter(capped == FALSE, optimised == "all_cap" & scenario == "PA") %>%
    mutate(catch_rule = "rfb", target = "PA", group = "rfb: all"),
  ### rfb-rule, default with cond. cap
  stats_rfb %>% 
    filter(capped == TRUE, optimised == "default" & scenario == "PA") %>%
    mutate(catch_rule = "rfb", target = "none", group = "rfb: default*"),
  ### rfb-rule, cond. cap, optimisation with multiplier, target PA
  stats_rfb %>% 
    filter(capped == TRUE, optimised == "mult" & scenario == "PA_capped") %>%
    mutate(catch_rule = "rfb", target = "PA", 
           group = "rfb: mult*"),
  ### rfb-rule, cond. cap, optimisation with all, target PA
  stats_rfb %>% 
    filter(capped == TRUE, optimised == "all" & scenario == "PA_capped") %>%
    mutate(catch_rule = "rfb", target = "PA", 
           group = "rfb: all*"),
  ### hr, default, target PA
  stats_hr %>%
    filter(objective == "MSY-PA" & scenario == "GA" & 
             parameters == "none") %>%
    mutate(catch_rule = "hr", target = "PA", group = "hr: default"),
  ### hr, optimisation with multiplier, target PA
  stats_hr %>%
    filter(objective == "MSY-PA" & scenario == "GA" & 
             parameters == "mult") %>%
    mutate(catch_rule = "hr", target = "PA", group = "hr: mult"),
  ### hr, optimisation with all parameters, target PA
  stats_hr %>%
    filter(objective == "MSY-PA" & scenario == "GA" & 
             parameters == "all") %>%
    mutate(catch_rule = "hr", target = "PA", group = "hr: all"),
  ### hr, cond. cap, default, target PA
  stats_hr %>%
    filter(objective == "MSY-PA" & scenario == "GA_cond_cap" & 
             parameters == "none") %>%
    mutate(catch_rule = "hr", target = "PA", group = "hr: default*"),
  ### hr, cond. cap, optimisation with multiplier, target PA
  stats_hr %>%
    filter(objective == "MSY-PA" & scenario == "GA_cond_cap" & 
             parameters == "mult") %>%
    mutate(catch_rule = "hr", target = "PA", group = "hr: mult*"),
  ### hr, cond. cap, optimisation with all parameters, target PA
  stats_hr %>%
    filter(objective == "MSY-PA" & scenario == "GA_cond_cap" & 
             parameters == "all") %>%
    mutate(catch_rule = "hr", target = "PA", group = "hr: all*")
)


breaks <- c(1, 3, 5, 6, 7, 8.5, 9.5, 10.5, 12.5, 13.5, 14.5, 16, 17, 18)
stats_plot_full <- stats_full %>%
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
           mutate(stock_k = factor(stock_k, levels = rev(stock_k)))) %>%
  mutate(
    group = factor(group, 
                   levels = c("zero-fishing", "2 over 3", 
                              "rfb: default", "rfb: mult", "rfb: all",
                              "rfb: default*", "rfb: mult*", "rfb: all*",
                              "hr: default", "hr: mult", "hr: all",
                              "hr: default*", "hr: mult*", "hr: all*"),
                   labels = c("(a) zero-fishing", "(b) 2 over 3", 
                              "(i) rfb: default",
                              "(j) rfb: mult", "(k) rfb: all",
                              "(c) rfb: default*",
                              "(d) rfb: mult*", "(e) rfb: all*",
                              "(l) hr: default",
                              "(m) hr: mult", "(n) hr: all",
                              "(f) hr: default*",
                              "(g) hr: mult*", "(h) hr: all*"))) %>%
  mutate(group_numeric = factor(group, labels = breaks)) %>%
  mutate(group_numeric = as.numeric(as.character(group_numeric)))

stats_plot_dev <- stats_plot_full %>%
  filter(name %in% c("SSB", "Catch")) %>%
  mutate(value_sign = ifelse(name == "SSB",
                             -abs(SSB_rel - 1)/2,
                             -abs(SSB_rel - 1) - abs(Catch_rel - 1)/2)) %>%
  mutate(sign = ifelse(name == "SSB",
                       ifelse(SSB_rel > 1, "+", "-"),
                       ifelse(Catch_rel > 1, "+", "-"))) %>%
  filter(abs(value) >= 0.1)

saveRDS(stats_plot_full, "output/data_all_comparison_full.rds")
stats_plot_full <- readRDS("output/data_all_comparison_full.rds")
saveRDS(stats_plot_dev, "output/data_all_comparison_full_dev.rds")
stats_plot_dev <- readRDS("output/data_all_comparison_full_dev.rds")

plot_comparison <- function(stock, data, data_dev, ylim, breaks) {
  data %>% 
    filter(stock %in% !!stock) %>%
    mutate(value = -abs(value)) %>%
    ggplot(aes(x = group_numeric, y = value, fill = name)) +
    geom_hline(yintercept = 0, linetype = "solid", size = 0.5, colour = "grey") +
    geom_col(position = "stack", width = 1, 
             colour = "black", size = 0.1, show.legend = TRUE) +
    ### for reversing order in legend (not in plot)
    scale_fill_discrete("fitness\nelements",
                        breaks = rev(levels(data$name))) +
    geom_text(data = stats_plot_dev %>%
                filter(stock %in% !!stock),
              aes(x = group_numeric, y = value_sign, label = sign),
              vjust = 0.5, colour = "grey20") +
    facet_grid(fhist ~ stock_k, scales = "free", space = "free_x", 
               labeller = label_parsed) +
    labs(y = "fitness", x = "") +
    scale_x_continuous(breaks = breaks, 
                       labels = levels(data$group),
                       minor_breaks = NULL) +
    coord_cartesian(ylim = ylim) +
    theme_bw(base_size = 8, base_family = "sans") +
    theme(panel.spacing.x = unit(0, units = "cm"),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5,
                                     lineheight = 0.7),
          #plot.margin = unit(x = c(1, 3, 3, 3), units = "pt"),
          legend.key.width = unit(0.65, "lines"),
          legend.key.height = unit(1, "lines"))
}

plot_comparison(stock = stocks$stock[12], data = stats_plot_full, 
                data_dev = data_dev, ylim = c(-8.5, 0), breaks = breaks)

plot_comparison_six <- function(stocks, data, data_dev, ylim, breaks) {
  #browser()
  if (isTRUE(length(stocks) == 5)) {
    p1 <- NULL
    stocks <- c(NA, stocks)
  } else {
    p1 <- plot_comparison(stock = stocks[1], data = data, 
                          data_dev = data_dev, ylim = ylim, breaks) + 
      theme(strip.text.y = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            legend.position = "none",) +
      labs(y = "")
  }
  p2 <- plot_comparison(stock = stocks[2], data = data, 
                        data_dev = data_dev, ylim = ylim, breaks) + 
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank())
  p3 <- plot_comparison(stock = stocks[3], data = data, 
                        data_dev = data_dev, ylim = ylim, breaks) + 
    theme(strip.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none") +
    labs(y = "fitness")
  p4 <- plot_comparison(stock = stocks[4], data = data, 
                        data_dev = data_dev, ylim = ylim, breaks) + 
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none",
          axis.title.y = element_blank())
  p5 <- plot_comparison(stock = stocks[5], data = data, 
                        data_dev = data_dev, ylim = ylim, breaks) + 
    theme(strip.text.y = element_blank(),
          legend.position = "none") +
    labs(y = "")
  p6 <- plot_comparison(stock = stocks[6], data = data, 
                        data_dev = data_dev, ylim = ylim, breaks) + 
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),legend.position = "none",
          axis.title.y = element_blank())
  ### combine all plots
  p <- plot_grid(p1, p2 + theme(legend.position = "none"), get_legend(p2),
                 p3, p4, NULL,
                 p5, p6, NULL,
                 ncol = 3, rel_heights = c(1, 1, 1.25), 
                 rel_widths = c(1, 1, 0.3))
  return(p)
}

### plot all stocks for supplementary material
for (i in seq_along(split(stocks$stock, ceiling(seq_along(stocks$stock)/6)))) {
  #browser()
  plot_comparison_six(stocks = split(stocks$stock, 
                                     ceiling(seq_along(stocks$stock)/6))[[i]],
                      data = stats_plot_full, data_dev = stats_plot_dev, 
                      ylim = c(-8.5, 0), breaks = breaks)
  ggsave(filename = paste0("output/plots/paper/all_comparison_full_", i , 
                           ".png"),
         width = 17, height = 18, units = "cm", dpi = 600, type = "cairo")
  ggsave(filename = paste0("output/plots/paper/all_comparison_full_", i , 
                           ".pdf"),
         width = 17, height = 18, units = "cm", dpi = 600)
}


### ------------------------------------------------------------------------ ###
### tmp: ang3 and cond. cap with multiplier ####
### ------------------------------------------------------------------------ ###

# runs <- readRDS("C:/Users/sf02/OneDrive - CEFAS/data-limited/GA_MSE_HR/output/hr/500_50/GA_cond_cap/random/ang3/multiplier-upper_constraint1.2-lower_constraint0.7--obj_SSB_C_risk_ICV_runs.rds")
# res <- readRDS("C:/Users/sf02/OneDrive - CEFAS/data-limited/GA_MSE_HR/output/hr/500_50/GA_cond_cap/random/ang3/multiplier-upper_constraint1.2-lower_constraint0.7--obj_SSB_C_risk_ICV_res.rds")
# res@solution
# 
# runs <- lapply(runs, function(x) {
#   #browser()
#   tmp <- c(x$pars, x$stats[1:11])
#   names(tmp) <- c(names(x$pars), dimnames(x$stats)[[1]][1:11])
#   tmp <- as.data.frame(tmp)
#   return(tmp)
# })
# runs <- do.call(bind_rows, runs)
# 
# 
# ggplot(data = runs,
#        aes(x = multiplier, y = risk_Blim)) +
#   geom_line() +
#   geom_hline(yintercept = 0.05) +
#   ylim(c(0, NA))

### ------------------------------------------------------------------------ ###
### supplementary tables: optimised parameterisations ####
### ------------------------------------------------------------------------ ###

pol_GA <- readRDS("output/hr_pol_comps_stats.rds")
stats_hr <- readRDS("output/hr_all_GA_stats.rds")

### default fitness
fitness_default <- stats_hr %>% 
  filter(parameters == "none") %>%
  group_by(objective, scenario, stock, fhist) %>%
  select(scenario, objective, fhist, stock, fitness) %>%
  rename(fitness_default = fitness) %>%
  ungroup() %>%
  arrange(scenario, objective, fhist)


### combine 
df <- bind_rows(
  ### pollack explorations
  pol_GA %>%
    filter(parameters != "not optimised") %>%
    mutate(stock = "pol", scenario_no = 0),
  # ### MSY mult
  # stats_hr %>%
  #   filter(objective == "MSY" & scenario == "GA" &
  #            parameters == "mult") %>%
  #   mutate(scenario_no = 0.1),
  # ### MSY all
  # stats_hr %>%
  #   filter(objective == "MSY" & scenario == "GA" &
  #            parameters == "all") %>%
  #   mutate(scenario_no = 0.2),
  ### MSY-PA mult
  stats_hr %>%
    filter(objective == "MSY-PA" & scenario == "GA" &
             parameters == "mult") %>%
    mutate(scenario_no = 1),
  ### MSY-PA all
  stats_hr %>%
    filter(objective == "MSY-PA" & scenario == "GA" &
             parameters == "all") %>%
    mutate(scenario_no = 2),
  ### MSY-PA mult & cond. cap
  stats_hr %>%
    filter(objective == "MSY-PA" & scenario == "GA_cond_cap" &
             parameters == "mult") %>%
    mutate(scenario_no = 3),
    ### MSY-PA all & cond. cap
  stats_hr %>%
  filter(objective == "MSY-PA" & scenario == "GA_cond_cap" &
             parameters == "all") %>%
    mutate(scenario_no = 4)
) %>%
  left_join(fitness_default) %>%
  mutate(fitness_improvement = round((1 - fitness/fitness_default)*100)) %>%
  select(scenario, scenario_no, objective, parameters, fhist, stock, iter, 
         multiplier, idxB_lag, idxB_range_3, comp_b_multiplier, 
         interval, upper_constraint, lower_constraint, 
         SSB_rel, Catch_rel, ICV, risk_Blim, 
         fitness_improvement) %>% 
  arrange(scenario_no, objective, fhist) %>%
  select(-scenario_no) %>%
  mutate(SSB_rel = round(SSB_rel, 2),
         Catch_rel = round(Catch_rel, 2),
         ICV = round(ICV, 2),
         risk_Blim = round(risk_Blim, 3))
# View(df)
write.csv(df, file = "output/hr_summary_table_parameters.csv", 
          row.names = FALSE)

### ------------------------------------------------------------------------ ###
### optimised multiplier  ####
### ------------------------------------------------------------------------ ###

all_mult <- readRDS("output/hr_all_mult_runs_stats.rds")


all_mult_plot <- all_mult %>%
  full_join(stocks %>%
            select(stock, k) %>%
            mutate(stock = factor(stock, levels = stock),
                   stock_k = paste0(stock, "~(italic(k)==", sprintf(k, fmt =  "%.2f"), 
                                    "*year^-1)")) %>%
            mutate(stock_k = factor(stock_k, levels = rev(stock_k))))

all_mult_plot %>%
  #filter(years == 100) %>%
  filter(k < 0.45 & k >= 0.32) %>%
  ggplot(aes(x = multiplier, y = risk_Blim, group = stock, colour = k)) +
  geom_line() +
  facet_grid(fhist ~ years) +
  theme_bw()

all_mult_plot %>%
  filter(optimum == TRUE) %>%
  ggplot(aes(x = k, y = multiplier, colour = fhist)) +
  geom_point() +
  theme_bw() +
  geom_vline(xintercept = 0.32) +
  geom_vline(xintercept = 0.45) +
  geom_smooth() +
  facet_grid(fhist ~ years)

all_mult_plot %>%
  filter(optimum == TRUE & years == 50) %>%
  ggplot(aes(x = k, y = multiplier, colour = fhist, shape = fhist)) +
  geom_line(aes(group = stock), colour = "grey", size = 0.3) +
  geom_point(size = 0.5) +
  #geom_jitter() + 
  geom_vline(xintercept = 0.315, size = 0.3, linetype = "dashed") +
  geom_vline(xintercept = 0.45, size = 0.3, linetype = "dashed") +
  #geom_smooth() +
  scale_colour_brewer("", palette = "Set1") +
  scale_shape_discrete("") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 2)) +
  labs(x = expression(italic(k)~'[year'^{-1}*']'), y = "Multiplier") + 
  theme_bw(base_size = 8) +
  theme(legend.position = c(0.75, 0.82),
        legend.key.width = unit(0.5, "lines"),
        legend.key.height = unit(0.5, "lines"),
        legend.key = element_blank(),
        legend.background = element_blank())
ggsave(filename = "output/plots/paper/all_mult.png", 
       width = 8.5, height = 6, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/paper/all_mult.pdf",
       width = 8.5, height = 6, units = "cm", dpi = 600)

### ------------------------------------------------------------------------ ###
### visualise fishing histories ####
### ------------------------------------------------------------------------ ###

### Fmsy relative to Fcrash for plotting (for pollack)
Fmsy_rel <- c(refpts(brps$pol)["msy", "harvest"]/
                refpts(brps$pol)["crash", "harvest"])
Fcrash <- c(refpts(brps$pol)["crash", "harvest"])
Bmsy_rel <- c(refpts(brps$pol)["msy", "ssb"]/
                refpts(brps$pol)["virgin", "ssb"])

### random scenario
set.seed(2)
n_iter <- 500
yrs_hist <- 100
start <- rep(0, n_iter)
middle <- runif(n = n_iter, min = 0, max = 1)
end <- runif(n = n_iter, min = 0, max = 1)
df <- t(sapply(seq(n_iter), 
  function(x) {
    c(approx(x = c(1, yrs_hist/2), 
             y = c(start[x], middle[x]), 
             n = yrs_hist/2)$y,
      approx(x = c(yrs_hist/2, yrs_hist + 1), 
             y = c(middle[x], end[x]), 
             n = (yrs_hist/2) + 1)$y[-1])
}))
df_rnd <- as.data.frame(t(df)) %>%
  mutate(year = -99:0) %>%
  pivot_longer(-year) %>%
  mutate(iter = gsub(x = name, pattern = "V", replacement = "")) %>%
  mutate(iter = as.numeric(iter)) %>%
  select(-name) %>%
  mutate(fhist = "random")

### one-way scenario
fs <- rep(Fmsy_rel * 0.5, 74)
f0 <- Fmsy_rel * 0.5
fmax <- 0.8
rate <- exp((log(fmax) - log(f0)) / (25))
fs <- c(fs, rate ^ (1:25) * f0)
df_ow <- data.frame(year = -98:0, value = fs, fhist = "one-way")

### roller-coaster
fs <- rep(Fmsy_rel * 0.5, 75)
f0_up <- Fmsy_rel * 0.5
fmax_up <- 0.8
yrs_up <- 15
rate_up <- exp((log(fmax_up) - log(f0_up)) / yrs_up)
yrs_down <- 6
f0_down <- Fmsy_rel
rate_down <- exp((log(fmax_up) - log(f0_down)) / yrs_down)
fs <- c(fs, rate_up ^ seq(yrs_up) * f0_up, rep(fmax_up, 3),
        rev(rate_down ^ seq(yrs_down) * f0_down))
df_rc <- data.frame(year = -98:0, value = c(fs), fhist = "roller-coaster")

### combine all fishing histories
df_all <- bind_rows(df_rnd, df_ow, df_rc) %>%
  mutate(fhist = factor(fhist, 
                        levels = c("one-way", "roller-coaster", "random")))

res_def_colours <- c("one-way" = brewer.pal(n = 4, name = "Set1")[1], 
                     "roller-coaster" = brewer.pal(n = 4, name = "Set1")[4], 
                     "random" = brewer.pal(n = 4, name = "Set1")[2])
p <- ggplot() +
  geom_line(data = df_all %>% filter(fhist == "random"),
            aes(x = year, y = value, group = iter, colour = fhist),
            size = 0.1, show.legend = FALSE, alpha = 0.2) +
  geom_line(data = df_all %>% filter(fhist != "random"),
            aes(x = year, y = value, group = iter, colour = fhist),
            size = 0.4, show.legend = FALSE) +
  facet_wrap(~ fhist) + 
  scale_colour_manual(values = res_def_colours) + 
  scale_x_continuous(breaks = c(-50, 0)) + 
  scale_y_continuous(sec.axis = sec_axis(trans = ~ ./Fmsy_rel,
                                         name = expression(F/F[MSY]))) + 
  coord_cartesian(xlim = c(-100, 0), ylim = c(0, 1), expand = FALSE) + 
  labs(y = expression(F/F[crash]), x = "Year") + 
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "pt"),
        axis.title.y.right = element_text(angle = 90))
ggsave(filename = "output/plots/paper/fhist.png", plot = p,
       width = 8.5, height = 4, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/paper/fhist.pdf", plot = p,
       width = 8.5, height = 4, units = "cm", dpi = 600)
### for paper: remove transparency with raster graphics
### 1st: save as postscript file, then convert to PDF
cairo_ps(filename = "output/plots/paper/fhist_ps_600.ps", 
         width = 8.5/2.54, height = 4/2.54, fallback_resolution = 600)
print(p)
dev.off()

### get SSB and F from pollack OM
pol_fhist <- foreach(fhist = c("one-way", "roller-coaster", "random"), 
                     .combine = bind_rows) %do% {
  #browser()
  om <- readRDS(paste0("input/hr/500_50/OM_1_hist/", fhist, "/pol.rds"))
  tmp <- as.data.frame(FLQuants(SSB = ssb(om$stk)/1000, 
                                Fbar = fbar(om$stk)/Fcrash)) %>%
    select(year, data, qname, iter) %>%
    mutate(qname = as.character(qname),
           fhist = fhist) %>%
    filter(year <= 100) %>%
    mutate(year = year - 100)
  return(tmp)
}
pol_fhist <- pol_fhist %>%
  mutate(fhist = factor(fhist, 
                        levels = c("one-way", "roller-coaster", "random")))

pol_fhist %>%
  ggplot(aes(x = year, y = data, group = iter)) +
  geom_line() +
  facet_grid(qname ~ fhist, scales = "free")

p_fbar <- pol_fhist %>%
  filter(qname == "Fbar") %>%
  ggplot(aes(x = year, y = data, group = iter, colour = fhist)) +
  geom_line(size = 0.1, show.legend = FALSE, alpha = 0.2) +
  facet_wrap(~ fhist) + 
  scale_colour_manual(values = res_def_colours) + 
  scale_x_continuous(breaks = c(-50, 0)) + 
  scale_y_continuous(sec.axis = sec_axis(trans = ~ ./Fmsy_rel,
                                         name = expression(F/F[MSY]))) + 
  coord_cartesian(xlim = c(-100, 0), ylim = c(0, 1), expand = FALSE) + 
  labs(y = expression(F/F[crash])) + 
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "pt"),
        axis.title.y.right = element_text(angle = 90),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_ssb <- pol_fhist %>%
  filter(qname == "SSB") %>%
  ggplot(aes(x = year, y = data, group = iter, colour = fhist)) +
  geom_line(size = 0.1, show.legend = FALSE, alpha = 0.2) +
  facet_wrap(~ fhist) + 
  scale_colour_manual(values = res_def_colours) + 
  scale_x_continuous(breaks = c(-50, 0)) + 
  scale_y_continuous(sec.axis = sec_axis(trans = ~ ./Bmsy_rel,
                                         name = expression(B/B[MSY]))) +
  coord_cartesian(xlim = c(-100, 0), ylim = c(0, 1.5), expand = FALSE) + 
  labs(y = expression(B/B[0]), x = "Year") + 
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "pt"),
        axis.title.y.right = element_text(angle = 90),
        strip.text = element_blank())

p <- p_fbar + p_ssb + plot_layout(ncol = 1)
ggsave(filename = "output/plots/paper/fhist_fbar_ssb.png", plot = p,
       width = 8.5, height = 7, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/paper/fhist_fbar_ssb.pdf", plot = p,
       width = 8.5, height = 7, units = "cm", dpi = 600)
