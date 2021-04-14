### ------------------------------------------------------------------------ ###
### analyse MSY GA runs ####
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

### ------------------------------------------------------------------------ ###
### plot fishing histories - pollack example ####
### ------------------------------------------------------------------------ ###

pol_ow <- readRDS("input/500_50/OM_1_hist/one-way/pol.rds")$stk
pol_rnd <- readRDS("input/500_50/OM_1_hist/random/pol.rds")$stk
pol_brp <- readRDS("input/brps.rds")$pol
refpts <- pol_brp@refpts
F_crash_msy <- c(refpts["crash", "harvest"]/refpts["msy", "harvest"])
B_0_msy <- c(refpts["virgin", "ssb"]/refpts["msy", "ssb"])

df_hist <- rbind(
    cbind(
      as.data.frame(FLQuants(Fbar = fbar(pol_ow)/c(refpts["msy", "harvest"]),
                    SSB = ssb(pol_ow)/c(refpts["msy", "ssb"]))),
      fhist = "one-way"),
    cbind(
      as.data.frame(FLQuants(Fbar = fbar(pol_rnd)/c(refpts["msy", "harvest"]),
                    SSB = ssb(pol_rnd)/c(refpts["msy", "ssb"]))),
      fhist = "random")
)
df_hist <- df_hist %>%
  filter(year <= 100 & year > 0) %>%
  mutate(data = ifelse(year %in% c(0, 1) & fhist == "one-way", NA, data)) %>%
  mutate(year = year - 100) %>%
  mutate(qname = as.character(qname)) %>%
  mutate(qname = factor(qname, levels = c("Fbar", "SSB"),
                        labels = c("F/F[MSY]", "SSB/B[MSY]"))) %>%
  mutate(fhist = as.character(fhist)) %>%
  mutate(fhist = factor(fhist, levels = c("one-way", "random"),
                        labels = c("one-way", "random"))) %>%
  mutate(type = "replicate")
df_hist <- df_hist %>% full_join(
  df_hist %>% 
    group_by(year, qname, fhist) %>%
    summarise(data = median(data, na.rm = TRUE)) %>%
    mutate(type = "median")
) %>%
  mutate(type = factor(type, levels = c("replicate", "median")))

attr(df_hist, "F_crash_msy") <- F_crash_msy
attr(df_hist, "B_0_msy") <- B_0_msy
saveRDS(df_hist, file = "output/plots/data_pol_fhist.rds")

df_hist <- readRDS("output/plots/data_pol_fhist.rds")
F_crash_msy <- attr(df_hist, "F_crash_msy")
B_0_msy <- attr(df_hist, "B_0_msy")

p_f <- df_hist %>%
  #filter(iter %in% c(1:50)) %>%
  filter(qname == "F/F[MSY]") %>%
  ggplot(aes(x = year, y = data, group = iter, colour = type, alpha = type,
             linetype = type, size = type)) +
  geom_line() +
  scale_alpha_manual("", values = c(replicate = 0.075, median = 1), ) +
  scale_colour_manual("", values = c(replicate = "black", median = "red")) +
  scale_linetype_manual("", values = c(replicate = "solid", median = "dashed")) +
  scale_size_manual("", values = c(replicate = 0.1, median = 0.5)) +
  geom_hline(yintercept = 1, size = 0.1, colour = "black", linetype = "dashed") +
  guides(alpha = FALSE) +
  theme_bw(base_size = 8) +
  facet_grid(~ fhist, scales = "free_y") +
  scale_y_continuous(
    sec.axis = dup_axis(trans = ~ . / F_crash_msy,
                        name = expression(F/F[crash])), 
    name = expression(F/F[MSY])) +
  theme(axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = unit(x = c(1, 3, 1, 3), units = "pt"),
        legend.position = c(0.2, 0.7),
        legend.background = element_blank(), 
        legend.key = element_blank(), legend.key.height = unit(0.3, "cm"))

p_ssb <- df_hist %>%
  #filter(iter %in% c(1:50)) %>%
  filter(qname == "SSB/B[MSY]") %>%
  ggplot(aes(x = year, y = data, group = iter, colour = type, alpha = type,
             linetype = type, size = type)) +
  geom_line(show.legend = FALSE) +
  scale_alpha_manual("", values = c(replicate = 0.075, median = 1), ) +
  scale_colour_manual("", values = c(replicate = "black", median = "red")) +
  scale_linetype_manual("", values = c(replicate = "solid", median = "dashed")) +
  scale_size_manual("", values = c(replicate = 0.1, median = 0.5)) +
  geom_hline(yintercept = 1, size = 0.1, colour = "black", linetype = "dashed") +
  theme_bw(base_size = 8) +
  facet_grid(~ fhist, scales = "free_y") +
  scale_y_continuous(
    sec.axis = dup_axis(trans = ~ . / B_0_msy,
                        name = expression(SSB/B[0])), 
    name = expression(SSB/B[MSY])) +
  # scale_x_continuous(breaks = c(-100, -75, -50, -25, 0), 
  #                    #labels = function(x) sub('^-', '\U2212', format(x))) +
  #                    labels = c("\U2212\U0031\U0030\U0030", # -100 
  #                               "\U2212\U0037\U0035",  # -75
  #                               "\U2212\U0035\U0030",  # -50
  #                               "\U2212\U0032\U0035",  # -25
  #                               "\U0030"           # 0
  #                               )) + 
  scale_x_continuous(breaks = c(-100, -75, -50, -25, 0), 
                     #labels = function(x) sub('^-', '\U2212', format(x))) +
                     labels = c("−100", "−75", "−50", "−25", "0"
                     )) + 
  theme(strip.text.x = element_blank(),
        plot.margin = unit(x = c(1, 3, 1, 3), units = "pt"))
p_both <- plot_grid(p_f, p_ssb, ncol = 1, align = "v")
p_both

ggsave(filename = "output/plots/pol_fhist.png", plot = p_both,
       width = 8.5, height = 6, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/pol_fhist.pdf", plot = p_both,
       width = 8.5, height = 6, units = "cm")
### Figure contains transparency
### for proper reproduction in journal, save as eps and convert transparent
### parts into raster graphic but keep rest as vector graphic
ggsave(filename = "output/plots/pol_fhist_600.eps", plot = p_both,
       width = 8.5, height = 6, units = "cm",
       device = cairo_ps, fallback_resolution = 600)
### convert with
###  gs -sOutputFile=pol_fhist_600.pdf -dNOPAUSE -sDEVICE=pdfwrite -dEPSCrop -dSAFER -dColorImageDownsampleType=/Bicubic -dColorImageResolution=600 -dColorImageDepth=4 -dDownsampleColorImages=true -dUseFlatCompression=true -dBATCH pol_fhist_600.eps
### embed fonts in Linux
# gs \
# -sOutputFile=pol_fhist_600_fonts.pdf \
# -dNOPAUSE -sDEVICE=pdfwrite -dEPSCrop \
# -dSAFER \
# -dColorImageDownsampleType=/Bicubic \
# -dColorImageResolution=600 \
# -dColorImageDepth=4 \
# -dDownsampleColorImages=true \
# -dUseFlatCompression=true -dBATCH \
# -c ".setpdfwrite <</NeverEmbed [ ]>> setdistillerparams" \
# -dSubsetFonts=true \
# -dCompatibilityLevel=1.4 \
# -f "pol_fhist_600.eps"

### ------------------------------------------------------------------------ ###
### collate results - pollack objective function explorations ####
### ------------------------------------------------------------------------ ###
n_yrs <- 50
n_iter <- 500

### get optimised parameterisation
pol_obj <- foreach(optimised = c(TRUE, FALSE), .combine = rbind) %:%
  foreach(scenario = "fitness_function", .combine = rbind) %:%
  foreach(stock = "pol", .combine = rbind) %:%
  foreach(fhist = c("one-way", "random"), .combine = rbind) %do% {
    #browser()
    ### find files
    fs <- list.files(path = paste0("output/", n_iter, "_", n_yrs, "/", scenario, 
                                   "/", fhist, "/pol/"), pattern = "*_res.rds")
    if (length(fs) < 1) return(NULL)
    fs <- fs[grep(x = fs, 
                  pattern = paste0("lag_idx-range_idx_1-range_idx_2-exp_r-exp_f",
                                   "-exp_b-interval-multiplier--obj_"))]
    trials <- data.frame(file = fs)
    trials$obj_fun <- gsub(x = trials$file, 
                           pattern = paste0("lag_idx-range_idx_1-range_idx_2-",
                                            "exp_r-exp_f-exp_b-interval-",
                                            "multiplier--obj_|_res\\.rds"), 
                           replacement = "")
    trials$fhist <- fhist
    trials$stock <- stock
    trials$scenario <- scenario
    trials$optimised <- optimised
    ### load GA results
    res_lst <- lapply(trials$file, function(x) {
      readRDS(paste0("output/", n_iter, "_", n_yrs, "/", scenario, "/", fhist, 
                     "/pol/", x))
    })
    res_par <- lapply(res_lst, function(x) {
      tmp <- x@solution[1,]
      tmp[c(1:4, 8)] <- round(tmp[c(1:4, 8)])
      tmp[5:7] <- round(tmp[5:7], 1)
      tmp[9] <- round(tmp[9], 2)
      tmp[10] <- ifelse(is.nan(tmp[10]), Inf, tmp[10])
      return(tmp)
    })
    names(res_par) <- trials$obj_fun
    ### default (non-optimised?)
    if (isFALSE(optimised)) {
      res_par <- lapply(res_par, function(x) {
        x[] <- c(1, 2, 3, 1, 1, 1, 1, 2, 1, Inf, 0)
        x})
    }
    res <- as.data.frame(do.call(rbind, res_par))
    res$file <- trials$file
    res$solution <- sapply(res_par, paste0, collapse = "_")
    res$fitness <- sapply(res_lst, slot, "fitnessValue")
    res$iter <- sapply(res_lst, slot, "iter")
    if (isFALSE(optimised)) {
      res$fitness <- NA
      res$iter <- NA
    }
    ### combine
    res <- merge(trials, res)
    ### load stats of solution
    res_stats <- lapply(split(res, seq(nrow(res))), function(x) {
      res_runs <- readRDS(paste0("output/", n_iter, "_", n_yrs, "/", scenario, 
                                 "/", fhist, "/pol/",
                                 gsub(x = x$file, pattern = "res", 
                                      replacement = "runs")))
      stats_i <- t(res_runs[[x$solution]]$stats)
      as.data.frame(lapply(as.data.frame(stats_i), unlist))
    })
    res_stats <- do.call(rbind, res_stats)
    res <- cbind(res, res_stats)
    if (isFALSE(optimised)) {
      for (i in 1:nrow(res)) {
        res$fitness[i] <- switch(res$obj_fun[i],
          "C"                = -abs(res$Catch_rel[i] - 1),
          "SSB"              = -abs(res$SSB_rel[i] - 1),
          "SSB_risk_ICV"     = -sum(abs(res$Catch_rel[i] - 1),
                                    res$ICV[i], res$risk_Blim[i]),
          "SSB_C_risk_ICV"   = -sum(abs(res$Catch_rel[i] - 1),
                                    abs(res$SSB_rel[i] - 1),
                                    res$ICV[i], res$risk_Blim[i]),
          "SSB_F_C_risk_ICV" = -sum(abs(res$Catch_rel[i] - 1),
                                    abs(res$SSB_rel[i] - 1),
                                    abs(res$Fbar_rel[i] - 1),
                                    res$ICV[i], res$risk_Blim[i]))
      }
    }
    return(res)
}
### fitness improvement
pol_obj_improvement <- pol_obj %>%
  select(obj_fun, fhist, optimised, fitness) %>%
  pivot_wider(names_from = optimised, values_from = fitness) %>%
  mutate(fitness_improvement = (`FALSE` - `TRUE`)/`FALSE` * 100,
         optimised = TRUE, `TRUE` = NULL, `FALSE` = NULL)
pol_obj <- pol_obj %>%
  left_join(pol_obj_improvement)
pol_obj_improvement %>%
  arrange(fhist) %>%
  mutate(fitness_improvement = round(fitness_improvement),
         optimised = NULL) %>%
  print(n = Inf)

saveRDS(pol_obj, file = "output/pol_obj_fun_explorations_stats.rds")
write.csv(pol_obj, file = "output/pol_obj_fun_explorations_stats.csv", 
          row.names = FALSE)

### ------------------------------------------------------------------------ ###
### collate results - all stocks with MSY fitness function ####
### ------------------------------------------------------------------------ ###
n_yrs <- 50
n_iter <- 500

### get optimised parameterisation
all_MSY <- foreach(stock = stocks$stock, .combine = rbind) %:%
  foreach(optimised = c(TRUE, FALSE), .combine = rbind) %:%
  foreach(scenario = "MSY", .combine = rbind) %:%
  foreach(fhist = c("one-way", "random"), .combine = rbind) %do% {#browser()
    #browser()
    ### find files
    fs <- list.files(path = paste0("output/", n_iter, "_", n_yrs, "/", scenario, 
                                   "/", fhist, "/", stock, "/"), 
                     pattern = paste0("lag_idx-range_idx_1-range_idx_2-exp_r-",
                                      "exp_f-exp_b-interval-multiplier--",
                                      "obj_SSB_C_risk_ICV_res.rds"))
    if (length(fs) < 1) return(NULL)
    trials <- data.frame(file = fs)
    trials$obj_fun <- gsub(x = trials$file, 
                           pattern = paste0("lag_idx-range_idx_1-range_idx_2-",
                                            "exp_r-exp_f-exp_b-interval-",
                                            "multiplier--obj_|_res\\.rds"), 
                           replacement = "")
    trials$fhist <- fhist
    trials$scenario <- scenario
    trials$optimised <- optimised
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
      return(tmp)
    })
    names(res_par) <- trials$obj_fun
    ### default (non-optimised?)
    if (isFALSE(optimised)) {
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
                                 gsub(x = x$file, pattern = "res", 
                                      replacement = "runs")))
      stats_i <- t(res_runs[[x$solution]]$stats)
      as.data.frame(lapply(as.data.frame(stats_i), unlist))
    })
    res_stats <- do.call(rbind, res_stats)
    res <- cbind(res, res_stats)
    ### calculate fitness value for non-optimised rule
    if (isFALSE(optimised)) {
      res$fitness <- sapply(split(res, seq(nrow(res))), function(x) {
        -sum(abs(x$SSB_rel - 1), abs(x$Catch_rel - 1), x$risk_Blim, x$ICV)
      })
    }
    
    return(res)
  }
### fitness improvement
all_MSY_improvement <- all_MSY %>%
  select(stock, fhist, optimised, fitness) %>%
  pivot_wider(names_from = optimised, values_from = fitness) %>%
  mutate(fitness_improvement = (`FALSE` - `TRUE`)/`FALSE` * 100,
         optimised = TRUE, `TRUE` = NULL, `FALSE` = NULL)
all_MSY <- all_MSY %>%
  left_join(all_MSY_improvement)
all_MSY_improvement %>%
  arrange(fhist) %>%
  mutate(fitness_improvement = round(fitness_improvement),
         optimised = NULL) %>%
  print(n = Inf)
saveRDS(all_MSY, file = "output/all_stocks_MSY_stats.rds")
write.csv(all_MSY, file = "output/all_stocks_MSY_stats.csv", row.names = FALSE)

### ------------------------------------------------------------------------ ###
### collate results - stock groups with MSY fitness function ####
### ------------------------------------------------------------------------ ###
n_yrs <- 50
n_iter <- 500

### get optimised parameterisation
groups_MSY <- foreach(group = c("low", "medium", "high"), .combine = rbind) %:%
  foreach(optimised = c(TRUE, FALSE), .combine = rbind) %:%
  foreach(scenario = "MSY", .combine = rbind) %:%
  foreach(fhist = c("one-way"), .combine = rbind) %do% {#browser()
    #browser()
    stocks_i <- switch(group,
                       "low" = stocks$stock[1:12],
                       "medium" = stocks$stock[13:20],
                       "high" = stocks$stock[21:29])
    stock <- paste0(stocks_i, collapse = "_")
    ### find files
    fs <- list.files(path = paste0("output/", n_iter, "_", n_yrs, "/", scenario, 
                                   "/", fhist, "/", stock, "/"), 
                     pattern = paste0("lag_idx-range_idx_1-range_idx_2-exp_r-",
                                      "exp_f-exp_b-interval-multiplier--",
                                      "obj_SSB_C_risk_ICV_res.rds"))
    if (length(fs) < 1) return(NULL)
    trials <- data.frame(file = fs)
    trials$obj_fun <- gsub(x = trials$file, 
                           pattern = paste0("lag_idx-range_idx_1-range_idx_2-",
                                            "exp_r-exp_f-exp_b-interval-",
                                            "multiplier--obj_|_res\\.rds"), 
                           replacement = "")
    trials$fhist <- fhist
    trials$scenario <- scenario
    trials$optimised <- optimised
    trials$group <- group
    trials <- trials[rep(1, length(stocks_i)), ]
    trials$stock <- stocks_i
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
      return(tmp)
    })
    names(res_par) <- trials$obj_fun
    ### default (non-optimised?)
    if (isFALSE(optimised)) {
      #trials <- trials[1, ]
      #res_lst <- res_lst[1]
      #res_par <- res_par[1]
      res_par <- lapply(res_par, function(x) {
        x[] <- c(1, 2, 3, 1, 1, 1, 1, 2, 1, Inf, 0)
        return(x)
      })
    }
    res <- as.data.frame(do.call(rbind, res_par))
    #res$file <- trials$file
    res$solution <- sapply(res_par, paste0, collapse = "_")
    res$fitness <- sapply(res_lst, slot, "fitnessValue")
    res$iter <- sapply(res_lst, slot, "iter")
    ### combine
    res <- cbind(trials, res)
    ### load stats of solution
    res_runs <- readRDS(paste0("output/", n_iter, "_", n_yrs, "/", scenario, 
                               "/", fhist, "/", stock, "/",
                               gsub(x = res$file[1], pattern = "res", 
                                    replacement = "runs")))
    stats_i <- t(res_runs[[res$solution[1]]]$stats)
    res_stats <- as.data.frame(lapply(as.data.frame(stats_i), unlist))
    #res_stats <- do.call(rbind, res_stats)
    res <- cbind(res, res_stats)
    ### add separate group results
    res_group <- res[1, ]
    res_group$stock <- "group"
    res_group[, c("risk_Blim", "risk_Bmsy", "risk_halfBmsy", "risk_collapse",
                  "SSB", "Fbar", "Catch", "SSB_rel", "Fbar_rel", "Catch_rel",
                  "ICV")] <- NA
    ### calculate fitness value for non-optimised rule
    if (isFALSE(optimised)) {
      res$fitness <- sapply(split(res, seq(nrow(res))), function(x) {
        -sum(abs(x$SSB_rel - 1), abs(x$Catch_rel - 1), x$risk_Blim, x$ICV)
      })
      res_group$fitness <- sum(res$fitness)
    }
    res$fitness <- NA
    res <- rbind(res_group, res)
    return(res)
  }
### fitness improvement
groups_MSY_improvement <- groups_MSY %>%
  select(stock, group, fhist, optimised, fitness) %>%
  pivot_wider(names_from = optimised, values_from = fitness) %>%
  mutate(fitness_improvement = (`FALSE` - `TRUE`)/`FALSE` * 100,
         optimised = TRUE, `TRUE` = NULL, `FALSE` = NULL)
groups_MSY <- groups_MSY %>%
  left_join(groups_MSY_improvement)
groups_MSY_improvement %>%
  arrange(fhist) %>%
  mutate(fitness_improvement = round(fitness_improvement),
         optimised = NULL) %>%
  print(n = Inf)

saveRDS(groups_MSY, file = "output/groups_MSY_stats.rds")
write.csv(groups_MSY, file = "output/groups_MSY_stats.csv", row.names = FALSE)

### ------------------------------------------------------------------------ ###
### recreate MSE projections for pollack explorations ####
### ------------------------------------------------------------------------ ###

### run optimised solution
library(doParallel)
n_cores <- ifelse(parallel::detectCores() > 20, 20, 10)
cl <- makeCluster(n_cores)
registerDoParallel(cl)
cl_length <- length(cl)
req_pckgs <- c("FLCore", "FLash", "mse", "GA", "doParallel", "doRNG", "FLBRP")
for (i in req_pckgs) library(package = i, character.only = TRUE)
source("funs.R"); source("funs_GA.R")
. <- foreach(i = seq(cl_length)) %dopar% {
  for (i in req_pckgs) library(package = i, character.only = TRUE)
  source("funs.R"); source("funs_GA.R")
}
### load stocks
stocks <- read.csv("input/stocks.csv")
stocks <- "pol"

### optimised parameters
pol_obj <- readRDS("output/pol_obj_fun_explorations_stats.rds")

for (obj_fun in seq(nrow(pol_obj))) {
  
  rm(res_mp)
  
  input <- readRDS(paste0("input/500_50/OM_2_mp_input/", pol_obj$fhist[obj_fun], 
                          "/", "pol.rds"))

  ### prepare input object for MSE
  input$args$nblocks <- cl_length
  input$cut_hist <- FALSE ### retain full history
  ### optimised parameters
  input$ctrl$est@args$idxB_lag     <- pol_obj$lag_idx[obj_fun]
  input$ctrl$est@args$idxB_range_1 <- pol_obj$range_idx_1[obj_fun]
  input$ctrl$est@args$idxB_range_2 <- pol_obj$range_idx_2[obj_fun]
  input$ctrl$est@args$catch_range  <- pol_obj$range_catch[obj_fun]
  input$ctrl$est@args$comp_m <- pol_obj$multiplier[obj_fun]
  input$ctrl$phcr@args$exp_r <- pol_obj$exp_r[obj_fun]
  input$ctrl$phcr@args$exp_f <- pol_obj$exp_f[obj_fun]
  input$ctrl$phcr@args$exp_b <- pol_obj$exp_b[obj_fun]
  input$ctrl$hcr@args$interval <- pol_obj$interval[obj_fun]
  input$ctrl$isys@args$interval <- pol_obj$interval[obj_fun]
  input$ctrl$isys@args$upper_constraint <- pol_obj$upper_constraint[obj_fun]
  input$ctrl$isys@args$lower_constraint <- pol_obj$lower_constraint[obj_fun]
  
  ### run MSE with optimised parameters
  res_mp <- do.call(mp, input)
  path_out <- paste0("output/500_50/fitness_function/", 
                     pol_obj$fhist[obj_fun], "/pol/")
  saveRDS(res_mp, file = paste0(path_out, "mp_", pol_obj[obj_fun,]$solution, ".rds"))
  
}

### ------------------------------------------------------------------------ ###
### pollack: objective function explorations - plots ####
### ------------------------------------------------------------------------ ###
pol_obj <- readRDS("output/pol_obj_fun_explorations_stats.rds")
pol_obj$seq <- seq(nrow(pol_obj))
pol_obj

### get time series
res <- lapply(pol_obj$seq, function(x) {
  res <- readRDS(paste0("output/500_50/fitness_function//", pol_obj$fhist[x], 
                        "/pol/", "mp_", pol_obj$solution[x], ".rds"))
  ### stock metrics
  SSBs <- FLCore::ssb(res@stock)
  Fs <- FLCore::fbar(res@stock)
  Cs <- FLCore::catch(res@stock)
  #yrs <- dim(SSBs)[2]
  its <- dim(SSBs)[6]
  yrs_hist <- ac(50:100)
  yrs_proj <- ac(101:150)
  yrs <- length(yrs_proj)
  ### collapse correction
  if (isTRUE(TRUE)) {
    ### find collapses
    cd <- sapply(seq(its), function(x) {
      min_yr <- min(which(SSBs[, yrs_proj,,,, x] < 1))
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
    SSBs[, yrs_proj]@.Data[cd] <- 0
    Cs[, yrs_proj]@.Data[cd] <- 0
    Fs[, yrs_proj]@.Data[cd] <- 0
  }
  input <- readRDS(paste0("input/500_50/OM_2_mp_input/", pol_obj$fhist[x], 
                          "/pol.rds"))
  Bmsy <- c(input$refpts["msy", "ssb"])
  Fmsy <- c(input$refpts["msy", "harvest"])
  Cmsy <- c(input$refpts["msy", "yield"])
  Blim <- input$Blim
  ### summarise
  qnts <- FLQuants(SSB = SSBs/Bmsy, F = Fs/Fmsy, Catch = Cs/Cmsy)
  qnts <- lapply(qnts, quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
  df_qnts <- as.data.frame(qnts)[, c("year", "iter", "data", "qname")]
  df_qnts <- bind_cols(df_qnts, pol_obj[rep(x, nrow(df_qnts)), ])
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
res_df$obj_fun[res_df$optimised == FALSE] <- "not optimised"
res_df$obj_fun[res_df$obj_fun == "C"] <- "Catch"
res_df$obj_fun[res_df$obj_fun == "SSB"] <- "SSB"
res_df$obj_fun[res_df$obj_fun == "SSB_C_risk_ICV"] <- "SSB+Catch+\nrisk+ICV"
res_df$obj_fun[res_df$obj_fun == "SSB_F_C_risk_ICV"] <- "SSB+F+Catch+\nrisk+ICV"
res_df$obj_fun[res_df$obj_fun == "SSB_risk_ICV"] <- "SSB+risk+ICV"
res_df$obj_fun <- as.factor(res_df$obj_fun)
res_df$obj_fun <- factor(res_df$obj_fun, levels = levels(res_df$obj_fun)[c(2, 1, 3, 6, 4, 5)])
res_df$`50%` <- ifelse(res_df$group == "history" & res_df$optimised == TRUE, 
                       NA, res_df$`50%`)
saveRDS(res_df, file = "output/plots/data_pol_fitness_trajectories.rds")
res_df <- readRDS("output/plots/data_pol_fitness_trajectories.rds")

plot_ssb_rel <- res_df %>% filter(qname == "SSB") %>%
  ggplot(aes(x = year, y = `50%`, colour = obj_fun, linetype = obj_fun)) +
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
  ggplot(aes(x = year, y = `50%`, colour = obj_fun, linetype = obj_fun)) +
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
  ggplot(aes(x = year, y = `50%`, colour = obj_fun, linetype = obj_fun)) +
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
ggsave(filename = "output/plots/pol_trajectories.png", 
       plot = p_pol_trajectories,
       width = 8.5, height = 10, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/pol_trajectories.pdf", 
       plot = p_pol_trajectories,
       width = 8.5, height = 10, units = "cm", dpi = 600)


### add GA progress to same plot
pol_ga <- readRDS(paste0("output/500_50/fitness_function/one-way/pol/",
                         "lag_idx-range_idx_1-range_idx_2-exp_r-exp_f-exp_b",
                         "-interval-multiplier--obj_SSB_C_risk_ICV_res.rds"))
ga_df <- as.data.frame(pol_ga@summary)
ga_df$generation <- seq(nrow(ga_df))
saveRDS(ga_df, file = "output/plots/data_pol_fitness_GA_progress.rds")
ga_df <- readRDS("output/plots/data_pol_fitness_GA_progress.rds")

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
  scale_x_continuous(expand = expand_scale(add = c(1, 1))) +
  scale_y_continuous(expand = expand_scale(add = c(0.4, abs(max(ga_df$max))))) +
  labs(x = "generations", y = "fitness value") +
  theme(legend.position = "right", legend.background = element_blank(),
        legend.key = element_blank(), 
        legend.key.size = unit(0.5, "lines"), legend.key.width = unit(1, "lines"),
        strip.text.x = element_blank(), legend.title = element_blank())
ggsave(filename = "output/plots/pol_ga.png", 
       plot = p_ga,
       width = 8.5, height = 4, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/pol_ga.pdf", 
       plot = p_ga,
       width = 8.5, height = 4, units = "cm", dpi = 600)




### plot stats
pol_obj

### format labels for plotting
pol_obj_plot <- pol_obj
pol_obj_plot$obj_label <- as.character(pol_obj_plot$obj_fun)
pol_obj_plot$obj_label[pol_obj_plot$optimised == FALSE] <- "not optimised"
pol_obj_plot$obj_label[pol_obj_plot$obj_label == "C"] <- "Catch"
pol_obj_plot$obj_label[pol_obj_plot$obj_label == "SSB"] <- "SSB"
pol_obj_plot$obj_label[pol_obj_plot$obj_label == "SSB_C_risk_ICV"] <- "SSB+Catch+\nrisk+ICV"
pol_obj_plot$obj_label[pol_obj_plot$obj_label == "SSB_F_C_risk_ICV"] <- "SSB+F+Catch+\nrisk+ICV"
pol_obj_plot$obj_label[pol_obj_plot$obj_label == "SSB_risk_ICV"] <- "SSB+risk+ICV"
pol_obj_plot$obj_label <- as.factor(pol_obj_plot$obj_label)
pol_obj_plot$obj_label <- factor(pol_obj_plot$obj_label, 
                              levels = levels(pol_obj_plot$obj_label)[c(2, 1, 3, 6, 4, 5)])

pol_obj_plot <- pol_obj_plot %>%
  pivot_longer(c(SSB_rel, Fbar_rel, Catch_rel, risk_Blim, ICV, fitness), 
               names_to = "key", values_to = "value") %>%
  mutate(stat = factor(key, levels = c("SSB_rel", "Fbar_rel", "Catch_rel", "risk_Blim", 
                                       "ICV", "fitness"), 
                       labels = c("SSB/B[MSY]", "F/F[MSY]", "Catch/MSY", 
                                  "B[lim]~risk", "ICV", "fitness~value")))
stats_targets <- data.frame(stat = c("SSB/B[MSY]", "F/F[MSY]", "Catch/MSY", 
                                     "B[lim]~risk", "ICV", "fitness~value"),
                            target = c(1, 1, 1, 0, 0, NA))

saveRDS(pol_obj_plot, file = "output/plots/data_pol_fitness_stats.rds")
pol_obj_plot <- readRDS("output/plots/data_pol_fitness_stats.rds")
saveRDS(stats_targets, file = "output/plots/data_pol_fitness_targets.rds")
stats_targets <- readRDS("output/plots/data_pol_fitness_targets.rds")

### individual plots
y_max <- 1.4
p_pol_stats_SSB <- pol_obj_plot %>% 
  filter(stat %in% c("SSB/B[MSY]")) %>%
  ggplot(aes(x = obj_label, y = value, fill = obj_label,
             colour = obj_label)) +
  geom_hline(yintercept = 1, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = position_dodge2(preserve = "single"), width = 0.8, 
           show.legend = FALSE, colour = "black", size = 0.1) +
  facet_grid(stat ~ fhist, scales = "free", space = "free_x", switch = "y",
             labeller = "label_parsed") +
  labs(y = "", x = "fitness function") +
    scale_fill_manual("fitness\nfunction", 
    values = setNames(c("black", "grey", scales::hue_pal()(4)),
                      c("not optimised", "SSB+Catch+\nrisk+ICV", "Catch", "SSB",
                        "SSB+risk+ICV", "SSB+F+Catch+\nrisk+ICV"))) +
  # scale_colour_manual("fitness\nfunction", 
  #   values = setNames(c("black", "grey", scales::hue_pal()(4)),
  #                     c("not optimised", "SSB+Catch+\nrisk+ICV", "Catch", "SSB",
  #                       "SSB+risk+ICV", "SSB+F+Catch+\nrisk+ICV"))) +
  theme_bw(base_size = 8, base_family = "sans") +
  theme(panel.spacing.x = unit(0, units = "cm"),
        strip.placement.y = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        plot.margin = unit(x = c(1, 3, 0, 3), units = "pt"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) +
  scale_y_continuous(trans = trans_from(), limits = c(0, y_max))
p_pol_stats_F <- pol_obj_plot %>% 
  filter(stat %in% c("F/F[MSY]")) %>%
  ggplot(aes(x = obj_label, y = value, fill = obj_label,
             colour = obj_label)) +
  geom_hline(yintercept = 1, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = "dodge", show.legend = FALSE, width = 0.8, 
           colour = "black", size = 0.1) +
  facet_grid(stat ~ fhist, scales = "free", space = "free_x", switch = "y",
             labeller = "label_parsed") +
  labs(y = "", x = "fitness function") +
    scale_fill_manual("fitness\nfunction", 
    values = setNames(c("black", "grey", scales::hue_pal()(4)),
                      c("not optimised", "SSB+Catch+\nrisk+ICV", "Catch", "SSB",
                        "SSB+risk+ICV", "SSB+F+Catch+\nrisk+ICV"))) +
  # scale_colour_manual("fitness\nfunction", 
  #   values = setNames(c("black", "grey", scales::hue_pal()(4)),
  #                     c("not optimised", "SSB+Catch+\nrisk+ICV", "Catch", "SSB",
  #                       "SSB+risk+ICV", "SSB+F+Catch+\nrisk+ICV"))) +
  theme_bw(base_size = 8, base_family = "sans") +
  theme(panel.spacing.x = unit(0, units = "cm"),
        strip.text.x = element_blank(),
        strip.placement.y = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        plot.margin = unit(x = c(0, 3, 0, 3), units = "pt"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) +
  scale_y_continuous(trans = trans_from(), limits = c(0, y_max))
p_pol_stats_C <- pol_obj_plot %>% 
  filter(stat %in% c("Catch/MSY")) %>%
  ggplot(aes(x = obj_label, y = value, fill = obj_label,
             colour = obj_label)) +
  geom_hline(yintercept = 1, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = "dodge", show.legend = FALSE, width = 0.8, 
           colour = "black", size = 0.1) +
  facet_grid(stat ~ fhist, scales = "free", space = "free_x", switch = "y",
             labeller = "label_parsed") +
  labs(y = "", x = "fitness function") +
    scale_fill_manual("fitness\nfunction", 
    values = setNames(c("black", "grey", scales::hue_pal()(4)),
                      c("not optimised", "SSB+Catch+\nrisk+ICV", "Catch", "SSB",
                        "SSB+risk+ICV", "SSB+F+Catch+\nrisk+ICV"))) +
  # scale_colour_manual("fitness\nfunction", 
  #   values = setNames(c("black", "grey", scales::hue_pal()(4)),
  #                     c("not optimised", "SSB+Catch+\nrisk+ICV", "Catch", "SSB",
  #                       "SSB+risk+ICV", "SSB+F+Catch+\nrisk+ICV"))) +
  theme_bw(base_size = 8, base_family = "sans") +
  theme(panel.spacing.x = unit(0, units = "cm"),
        strip.text.x = element_blank(),
        strip.placement.y = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        plot.margin = unit(x = c(0, 3, 0, 3), units = "pt"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) +
  scale_y_continuous(trans = trans_from(), limits = c(0, y_max))
p_pol_stats_risk <- pol_obj_plot %>% 
  filter(stat %in% c("B[lim]~risk")) %>%
  ggplot(aes(x = obj_label, y = value, fill = obj_label,
             colour = obj_label)) +
  geom_hline(yintercept = 0, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = "dodge", show.legend = FALSE, width = 0.8, 
           colour = "black", size = 0.1) +
  facet_grid(stat ~ fhist, scales = "free", space = "free_x", switch = "y",
             labeller = "label_parsed") +
  labs(y = "", x = "fitness function") +
    scale_fill_manual("fitness\nfunction", 
    values = setNames(c("black", "grey", scales::hue_pal()(4)),
                      c("not optimised", "SSB+Catch+\nrisk+ICV", "Catch", "SSB",
                        "SSB+risk+ICV", "SSB+F+Catch+\nrisk+ICV"))) +
  # scale_colour_manual("fitness\nfunction", 
  #   values = setNames(c("black", "grey", scales::hue_pal()(4)),
  #                     c("not optimised", "SSB+Catch+\nrisk+ICV", "Catch", "SSB",
  #                       "SSB+risk+ICV", "SSB+F+Catch+\nrisk+ICV"))) +
  theme_bw(base_size = 8, base_family = "sans") +
  theme(panel.spacing.x = unit(0, units = "cm"),
        strip.text.x = element_blank(),
        strip.placement.y = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        plot.margin = unit(x = c(0, 3, 0, 3), units = "pt"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) +
  scale_y_continuous(trans = trans_from(0), limits = c(0, y_max))
p_pol_stats_ICV <- pol_obj_plot %>% 
  filter(stat %in% c("ICV")) %>%
  ggplot(aes(x = obj_label, y = value, fill = obj_label,
             colour = obj_label)) +
  geom_hline(yintercept = 0, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = "dodge", show.legend = FALSE, width = 0.8, 
           colour = "black", size = 0.1) +
  facet_grid(stat ~ fhist, scales = "free", space = "free_x", switch = "y",
             labeller = "label_parsed") +
  labs(y = "", x = "fitness function") +
    scale_fill_manual("fitness\nfunction", 
    values = setNames(c("black", "grey", scales::hue_pal()(4)),
                      c("not optimised", "SSB+Catch+\nrisk+ICV", "Catch", "SSB",
                        "SSB+risk+ICV", "SSB+F+Catch+\nrisk+ICV"))) +
  # scale_colour_manual("fitness\nfunction", 
  #   values = setNames(c("black", "grey", scales::hue_pal()(4)),
  #                     c("not optimised", "SSB+Catch+\nrisk+ICV", "Catch", "SSB",
  #                       "SSB+risk+ICV", "SSB+F+Catch+\nrisk+ICV"))) +
  theme_bw(base_size = 8, base_family = "sans") +
  theme(panel.spacing.x = unit(0, units = "cm"),
        strip.text.x = element_blank(),
        strip.placement.y = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        plot.margin = unit(x = c(0, 3, 0, 3), units = "pt"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) +
  scale_y_continuous(trans = trans_from(0), limits = c(0, y_max))
p_pol_stats_fitness <- pol_obj_plot %>% 
  filter(stat %in% c("fitness~value")) %>%
  ggplot(aes(x = obj_label, y = value, fill = obj_label,
             colour = obj_label)) +
  geom_hline(yintercept = 0, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = "dodge", show.legend = FALSE, width = 0.8, 
           colour = "black", size = 0.1) +
  facet_grid(stat ~ fhist, scales = "free", space = "free_x", switch = "y",
             labeller = "label_parsed") +
  labs(y = "", x = "fitness function") +
    scale_fill_manual("fitness\nfunction", 
    values = setNames(c("black", "grey", scales::hue_pal()(4)),
                      c("not optimised", "SSB+Catch+\nrisk+ICV", "Catch", "SSB",
                        "SSB+risk+ICV", "SSB+F+Catch+\nrisk+ICV"))) +
  # scale_colour_manual("fitness\nfunction", 
  #   values = setNames(c("black", "grey", scales::hue_pal()(4)),
  #                     c("not optimised", "SSB+Catch+\nrisk+ICV", "Catch", "SSB",
  #                       "SSB+risk+ICV", "SSB+F+Catch+\nrisk+ICV"))) +
  theme_bw(base_size = 8, base_family = "sans") +
  theme(panel.spacing.x = unit(0, units = "cm"),
        strip.text.x = element_blank(),
        strip.placement.y = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        plot.margin = unit(x = c(0, 3, 3, 3), units = "pt")) +
  scale_y_continuous(trans = trans_from(0), limits = c(-y_max, NA),
                     breaks = c(0, -0.5, -1), 
                     minor_breaks = c(-0.25, -0.75, -1.25))
p_pol_stats_comb <- plot_grid(p_pol_stats_SSB, p_pol_stats_F, p_pol_stats_C,
                              p_pol_stats_risk, p_pol_stats_ICV,
                              p_pol_stats_fitness,
                              ncol = 1, align = "v",
                              rel_heights = c(1.25, 1, 1, 1, 1, 2.1))
ggsave(filename = "output/plots/pol_stats.png", 
       plot = p_pol_stats_comb,
       width = 8.5, height = 13, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/pol_stats.pdf", 
       plot = p_pol_stats_comb,
       width = 8.5, height = 13, units = "cm", dpi = 600)



### combine all plots
p <- plot_grid(
  plot_grid(p_pol_trajectories, 
            plot_grid(plot_grid(NULL, p_ga + theme(legend.position = "none"),
                                rel_widths = c(0.07, 1), nrow = 1), 
                      get_legend(p_ga), 
                      nrow = 1, rel_widths = c(1, 0.39)), 
            ncol = 1, rel_heights = c(4, 1), 
            labels = c("(a)", "(b)"), label_size = 10, hjust = -0.1, vjust = 1.1), 
  p_pol_stats_comb, ncol = 2, rel_widths = c(1.2, 1),
  labels = c("", "(c)"), label_size = 10, hjust = -0.1, vjust = 1.1)
ggsave(filename = "output/plots/pol_combined.png", plot = p,
       width = 17, height = 13, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/pol_combined.pdf", plot = p,
       width = 17, height = 13, units = "cm", dpi = 600)
### Figure contains transparency
### for proper reproduction in journal, save as eps and convert transparent
### parts into raster graphic but keep rest as vector graphic
ggsave(filename = "output/plots/pol_combined_2400.eps", plot = p,
       width = 17, height = 13, units = "cm",
       device = cairo_ps, fallback_resolution = 2400)
### convert with
###  gs -sOutputFile=pol_combined_2400.pdf -dNOPAUSE -sDEVICE=pdfwrite -dEPSCrop -dSAFER -dColorImageDownsampleType=/Bicubic -dColorImageResolution=2400 -dColorImageDepth=4 -dDownsampleColorImages=true -dUseFlatCompression=true -dBATCH pol_combined_2400.eps

### ------------------------------------------------------------------------ ###
### pollack: check impact of fixing TAC interval ####
### ------------------------------------------------------------------------ ###

### get optimised parameterisation
pol_interval <- foreach(stock = "pol", .combine = rbind) %:%
  foreach(optimised = c(TRUE, FALSE), .combine = rbind) %:%
  foreach(scenario = "fitness_function", .combine = rbind) %:%
  foreach(interval = 1:3, .combine = rbind) %:%
  foreach(fhist = c("one-way"), .combine = rbind) %do% {#browser()
    res_ga <- readRDS(paste0("output/500_50/", scenario, "/", fhist, "/", stock, 
                          "/lag_idx-range_idx_1-range_idx_2-exp_r-exp_f-exp_b",
                          "-interval", interval, "-multiplier--obj_SSB_C_risk",
                          "_ICV_res.rds"))
    trials <- data.frame(file = NA)
    trials$obj_fun <- "SSB_Catch_risk_ICV"
    trials$fhist <- fhist
    trials$scenario <- scenario
    trials$optimised <- optimised
    trials$stock <- stock
    ### load GA results
    res_par <- res_ga@solution[1,]
    res_par[c(1:4, 8)] <- round(res_par[c(1:4, 8)])
    res_par[5:7] <- round(res_par[5:7], 1)
    res_par[9] <- round(res_par[9], 2)
    res_par[10] <- ifelse(is.nan(res_par[10]), Inf, res_par[10])
    ### default (non-optimised?)
    if (isFALSE(optimised)) {
      res_par[] <- c(1, 2, 3, 1, 1, 1, 1, interval, 1, Inf, 0)
    }
    res <- cbind(trials, as.data.frame(t(res_par)))
    res$solution <- paste0(res_par, collapse = "_")
    res$fitness <- res_ga@fitnessValue
    res$iter <- res_ga@iter
    ### load stats of solution
    runs <- readRDS(paste0("output/500_50/", scenario, "/", fhist, "/", stock, 
                          "/lag_idx-range_idx_1-range_idx_2-exp_r-exp_f-exp_b",
                          "-interval", interval, "-multiplier--obj_SSB_C_risk",
                          "_ICV_runs.rds"))
    stats <- t(runs[[res$solution]]$stats)
    stats <-  as.data.frame(lapply(as.data.frame(stats), unlist))
    res <- cbind(res, stats)
    ### calculate fitness value for non-optimised rule
    if (isFALSE(optimised)) {
      res$fitness <- -sum(abs(res$SSB_rel - 1), 
                          abs(res$Catch_rel - 1), 
                          res$risk_Blim, res$ICV)
    }
    return(res)
}
saveRDS(pol_interval, file = "output/pol_interval_MSY_stats.rds")
write.csv(pol_interval, file = "output/pol_interval_MSY_stats.csv", 
          row.names = FALSE)
pol_interval <- readRDS("output/pol_interval_MSY_stats.rds")

### check values
pol_interval %>% 
  filter(optimised == TRUE) %>%
  select(stock, interval, fitness) %>% 
  filter(interval == 2) %>%
  transmute(fitness2 = fitness, stock = stock) %>%
  full_join(pol_interval %>% 
    filter(optimised == TRUE) %>%
    select(stock, interval, fitness)) %>%
  mutate(change = (1 - fitness/fitness2)*100)

### ------------------------------------------------------------------------ ###
### all stocks stats: default vs. optimised ####
### ------------------------------------------------------------------------ ###

### load data
all_MSY <- readRDS("output/all_stocks_MSY_stats.rds")
stocks_sorted <- stocks %>%
  select(stock, k) %>%
  arrange(k)

### format for plotting
stats_plot <- all_MSY %>%
  filter(fhist == "one-way") %>%
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
         rule = ifelse(optimised == TRUE, "optimised", "default"))
saveRDS(stats_plot, file = "output/plots/data_stocks_stats.rds")
stats_plot <- readRDS("output/plots/data_stocks_stats.rds")

### individual plots
p_stats_SSB <- stats_plot %>%
  filter(stat %in% c("SSB/B[MSY]")) %>%
  ggplot(aes(x = stock, y = value, fill = rule, colour = rule)) +
  geom_hline(yintercept = 1, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = position_dodge2(), width = 0.8, 
           show.legend = FALSE, colour = "black", size = 0.1) +
  scale_fill_manual("rfb rule", 
                    values = c(default = "black", optimised = "grey")) + 
  facet_grid(stat ~ "stock~specific~optimisation", 
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
           show.legend = FALSE, colour = "black", size = 0.1) +
  scale_fill_manual("rfb rule", 
                    values = c(default = "black", optimised = "grey")) + 
  facet_grid(stat ~ "stock~specific~optimisation", 
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
           show.legend = FALSE, colour = "black", size = 0.1) +
  scale_fill_manual("rfb rule", 
                    values = c(default = "black", optimised = "grey")) + 
  facet_grid(stat ~ "stock~specific~optimisation", 
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
  geom_col(position = position_dodge2(), width = 0.8, 
           show.legend = FALSE, colour = "black", size = 0.1) +
  scale_fill_manual("rfb rule", 
                    values = c(default = "black", optimised = "grey")) + 
  facet_grid(stat ~ "stock~specific~optimisation", 
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
           show.legend = FALSE, colour = "black", size = 0.1) +
  scale_fill_manual("rfb rule", 
                    values = c(default = "black", optimised = "grey")) + 
  facet_grid(stat ~ "stock~specific~optimisation", 
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
           colour = "black", size = 0.1, show.legend = FALSE) +
  scale_fill_manual("rfb rule", 
                    values = c(default = "black", optimised = "grey")) + 
  facet_grid(stat ~ "stock~specific~optimisation", 
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
        axis.title.y = element_blank()) +
  scale_y_continuous(trans = trans_from(0), limits = c(NA, NA))
p_stats_combined <- plot_grid(p_stats_SSB, p_stats_F, p_stats_C,
                              p_stats_risk, p_stats_ICV,
                              p_stats_fitness,
                              ncol = 1, align = "v",
                              rel_heights = c(1.25, 1, 1, 1, 1, 1.5))
ggsave(filename = "output/plots/stats_single.png", 
       plot = p_stats_combined,
       width = 15, height = 13, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/stats_single.pdf", 
       plot = p_stats_combined,
       width = 15, height = 13, units = "cm", dpi = 600)


### ------------------------------------------------------------------------ ###
### k-groups stats: default vs. optimised ####
### ------------------------------------------------------------------------ ###

groups_MSY <- readRDS("output/groups_MSY_stats.rds")
stocks_sorted <- stocks %>%
  select(stock, k) %>%
  arrange(k) %>%
  bind_rows(data.frame(stock = "group", k = Inf))

### format for plotting
stats_group_plot <- groups_MSY %>%
  filter(fhist == "one-way") %>%
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
         rule = ifelse(optimised == TRUE, "optimised", "default"),
         k_group = paste0(group, "-italic(k)~group"),
         k_group = factor(k_group, levels = unique(k_group)[c(1, 2, 3)]))
saveRDS(stats_group_plot, file = "output/plots/data_groups_stats.rds")
stats_group_plot <- readRDS("output/plots/data_groups_stats.rds")

### individual plots
p_stats_ms_SSB <- stats_group_plot %>%
  filter(stat %in% c("SSB/B[MSY]")) %>%
  ggplot(aes(x = stock, y = value, fill = rule, colour = rule)) +
  geom_hline(yintercept = 1, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = position_dodge2(), width = 0.8, 
           show.legend = FALSE, colour = "black", size = 0.1) +
  scale_fill_manual("rfb rule", 
                    values = c(default = "black", optimised = "grey")) + 
  facet_grid(stat ~ k_group, space = "free_x", switch = "y",
             scales = "free",# strip.position = "left", 
             labeller = "label_parsed") +
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
p_stats_ms_F <- stats_group_plot %>%
  filter(stat %in% c("F/F[MSY]")) %>%
  ggplot(aes(x = stock, y = value, fill = rule, colour = rule)) +
  geom_hline(yintercept = 1, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = position_dodge2(), width = 0.8, 
           show.legend = FALSE, colour = "black", size = 0.1) +
  scale_fill_manual("rfb rule", 
                    values = c(default = "black", optimised = "grey")) + 
  facet_grid(stat ~ k_group, space = "free_x", switch = "y",
             scales = "free",# strip.position = "left", 
             labeller = "label_parsed") +
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
p_stats_ms_C <- stats_group_plot %>%
  filter(stat %in% c("Catch/MSY")) %>%
  ggplot(aes(x = stock, y = value, fill = rule, colour = rule)) +
  geom_hline(yintercept = 1, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = position_dodge2(), width = 0.8, 
           show.legend = FALSE, colour = "black", size = 0.1) +
  scale_fill_manual("rfb rule", 
                    values = c(default = "black", optimised = "grey")) + 
  facet_grid(stat ~ k_group, space = "free_x", switch = "y",
             scales = "free",# strip.position = "left", 
             labeller = "label_parsed") +
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
p_stats_ms_risk <- stats_group_plot %>%
  filter(stat %in% c("B[lim]~risk")) %>%
  ggplot(aes(x = stock, y = value, fill = rule, colour = rule)) +
  geom_hline(yintercept = 0, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = position_dodge2(), width = 0.8, 
           show.legend = FALSE, colour = "black", size = 0.1) +
  scale_fill_manual("rfb rule", 
                    values = c(default = "black", optimised = "grey")) + 
  facet_grid(stat ~ k_group, space = "free_x", switch = "y",
             scales = "free",# strip.position = "left", 
             labeller = "label_parsed") +
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
p_stats_ms_ICV <- stats_group_plot %>%
  filter(stat %in% c("ICV")) %>%
  ggplot(aes(x = stock, y = value, fill = rule, colour = rule)) +
  geom_hline(yintercept = 0, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = position_dodge2(), width = 0.8, 
           show.legend = FALSE, colour = "black", size = 0.1) +
  scale_fill_manual("rfb rule", 
                    values = c(default = "black", optimised = "grey")) + 
  facet_grid(stat ~ k_group, space = "free_x", switch = "y",
             scales = "free",# strip.position = "left", 
             labeller = "label_parsed") +
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
p_stats_ms_fitness <- stats_group_plot %>%
  filter(stat %in% c("fitness~value")) %>%
  ggplot(aes(x = stock, y = value, fill = rule, colour = rule)) +
  geom_hline(yintercept = 0, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = position_dodge2(), width = 0.8, 
           colour = "black", size = 0.1) +
  scale_fill_manual("rfb-rule", 
                    values = c(default = "black", optimised = "grey")) + 
  facet_grid(stat ~ k_group, space = "free_x", switch = "y",
             scales = "free",# strip.position = "left", 
             labeller = "label_parsed") +
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
        legend.key.size = unit(0.2, "lines"),
        legend.position = c(0.15, 0.5),
        legend.background = element_blank()) +
  scale_y_continuous(trans = trans_from(0), limits = c(NA, NA))
p_stats_ms_combined <- plot_grid(p_stats_ms_SSB, p_stats_ms_F, p_stats_ms_C,
                              p_stats_ms_risk, p_stats_ms_ICV,
                              p_stats_ms_fitness,
                              ncol = 1, align = "v",
                              rel_heights = c(1.25, 1, 1, 1, 1, 1.5))

### combine single and multi-species stats plots ####
plot_grid(plot_grid(NULL, p_stats_combined + theme(legend.position = "none"),
                    ncol = 1, rel_heights = c(0, 1)),#c(0.0395, 1)), 
          plot_grid(p_stats_ms_combined,# + theme(legend.position = "none"),
                    NULL,#get_legend(p_stats_ms),
                    ncol = 1, rel_heights = c(1, 0)),#0.17)),
          labels = c("(a)", "(b)"), label_size = 10, rel_widths = c(0.95, 1.05))
ggsave(filename = "output/plots/stats_ms_and_single.png",
       width = 17, height = 12, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/stats_ms_and_single.pdf",
       width = 17, height = 12, units = "cm")

### compare stock fitness: stock-specific vs groups ####
# stats <- readRDS("output/all_stocks_MSY_stats.rds")
# stats_ms <- readRDS("output/groups_MSY_stats.rds")
# stats_tmp <- stats_ms %>%
#   rename(fitness_ms = fitness,
#          iter_ms = iter) %>%
#   select(fhist, optimised, stock, group, iter_ms, fitness_ms) %>%
#   left_join(stats) %>%
#   select(optimised, stock, group, iter, iter_ms, fitness, 
#          fitness_ms) %>%
#   mutate(fitness_ratio = fitness_ms/fitness)
# stats_tmp
# ### check fitness per group
# stats_tmp %>% group_by(optimised, group) %>%
#   summarise(fitness = sum(fitness, na.rm = TRUE),
#             fitness_ms = sum(fitness_ms, na.rm = TRUE))

### ------------------------------------------------------------------------ ###
### ICES default 2 over 3 rule ####
### ------------------------------------------------------------------------ ###

### stats from rfb-rule: default and optimised
stats_rfb <- readRDS("output/all_stocks_MSY_stats.rds")
stats_rfb <- stats_rfb %>%
  mutate(catch_rule = ifelse(optimised, "optimised", "default"))
### stats from rfb-rule with 20% constraint
stats_rfb_constrained <- foreach(stock = stocks$stock, .combine = "rbind") %:%
  foreach(fhist = c("one-way", "random"), .combine = "rbind") %do% {
    res_i <- readRDS(paste0("output/500_50/MSY/", fhist, "/", stock, 
                              "/1_2_3_1_1_1_1_2_1_1.2_0.8.rds"))
    stats_i <- data.frame(file = NA, obj_fun = NA, fhist = fhist, 
                          scenario = "MSY", optimised = FALSE, stock = stock,
                          lag_idx = 1, range_idx_1 = 2, range_idx_2 = 3,
                          range_catch = 1, exp_r = 1, exp_f = 1, exp_b = 1,
                          interval = 2, multiplier = 1, upper_constraint = 1.2,
                          lower_constraint = 0.8, solution = NA, fitness = NA,
                          iter = NA,
                          stringsAsFactors = FALSE)
    res_i <- as.data.frame(lapply(as.data.frame(t(res_i)), unlist))
    res_i <- res_i[, 1:11]
    stats_i <- cbind(stats_i, res_i)
    stats_i$fhist <- fhist
    stats_i$stock = stock
    stats_i$catch_rule = "constrained"
    stats_i$fitness <- -sum(abs(stats_i$SSB_rel - 1), 
                            abs(stats_i$Catch_rel - 1), 
                            stats_i$risk_Blim, stats_i$ICV)
    return(stats_i)
}
### stats from 2 over 3 rule
stats_2over3 <- foreach(stock = stocks$stock, .combine = "rbind") %:%
  foreach(fhist = c("one-way", "random"), .combine = "rbind") %do% {
  
    stats_i <- readRDS(paste0("output/500_50/2over3/", fhist, "/", stock, 
                              "_stats.rds"))
    stats_i <- as.data.frame(lapply(as.data.frame(t(stats_i)), unlist))
    stats_i$fhist <- fhist
    stats_i$stock = stock
    stats_i$catch_rule = "2 over 3"
    stats_i$fitness <- -sum(abs(stats_i$SSB_rel - 1), 
                            abs(stats_i$Catch_rel - 1), 
                            stats_i$risk_Blim, stats_i$ICV)
    return(stats_i)
}

### combine 
stats_plot <- bind_rows(stats_rfb %>%
                          select(fhist, catch_rule, fitness, stock, 
                                 risk_Blim:ICV), 
                        stats_2over3, stats_rfb_constrained)

stats_plot <- stats_plot %>% 
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
         fhist = factor(fhist, levels = c("one-way", "random")),
         stock = factor(stock, level = stocks$stock),
         catch_rule = factor(catch_rule, 
                             levels = c("2 over 3", "default", "constrained", 
                                        "optimised"), 
                             labels = c("2 over 3", "rfb\n(default)",
                                        "rfb\n(constrained)", 
                                        "rfb\n(optimised)")))
saveRDS(stats_plot, file = "output/plots/data_2over3_stats.rds")
stats_plot <- readRDS("output/plots/data_2over3_stats.rds")
stats_plot_full <- stats_plot
stats_plot <- stats_plot %>%
  filter(!catch_rule %in% "rfb\n(constrained)")

### individual plots
p_stats_23_SSB <- stats_plot %>%
  filter(stat %in% c("SSB/B[MSY]")) %>%
  ggplot(aes(x = stock, y = value, fill = catch_rule)) +
  geom_hline(yintercept = 1, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = position_dodge2(preserve = "single"), width = 0.8, 
           show.legend = FALSE, colour = "black", size = 0.1) +
  scale_fill_manual("catch\nrule",
                    values = c("rfb\n(default)" = "black", 
                               "rfb\n(optimised)" = "grey",
                               "2 over 3" = "white")) +
  facet_grid(stat ~ fhist, space = "free_x", switch = "y",
             scales = "free_y", labeller = "label_parsed") +
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
p_stats_23_F <- stats_plot %>%
  filter(stat %in% c("F/F[MSY]")) %>%
  ggplot(aes(x = stock, y = value, fill = catch_rule)) +
  geom_hline(yintercept = 1, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = position_dodge2(preserve = "single"), width = 0.8, 
           show.legend = FALSE, colour = "black", size = 0.1) +
  scale_fill_manual("catch\nrule",
                    values = c("rfb\n(default)" = "black", 
                               "rfb\n(optimised)" = "grey",
                               "2 over 3" = "white")) +
  facet_grid(stat ~ fhist, space = "free_x", switch = "y",
             scales = "free_y", labeller = "label_parsed") +
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
p_stats_23_C <- stats_plot %>%
  filter(stat %in% c("Catch/MSY")) %>%
  ggplot(aes(x = stock, y = value, fill = catch_rule)) +
  geom_hline(yintercept = 1, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = position_dodge2(preserve = "single"), width = 0.8, 
           show.legend = FALSE, colour = "black", size = 0.1) +
  scale_fill_manual("catch\nrule",
                    values = c("rfb\n(default)" = "black", 
                               "rfb\n(optimised)" = "grey",
                               "2 over 3" = "white")) +
  facet_grid(stat ~ fhist, space = "free_x", switch = "y",
             scales = "free_y", labeller = "label_parsed") +
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
p_stats_23_risk <- stats_plot %>%
  filter(stat %in% c("B[lim]~risk")) %>%
  ggplot(aes(x = stock, y = value, fill = catch_rule)) +
  geom_hline(yintercept = 0, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = position_dodge2(preserve = "single"), width = 0.8, 
           show.legend = FALSE, colour = "black", size = 0.1) +
  scale_fill_manual("catch\nrule",
                    values = c("rfb\n(default)" = "black", 
                               "rfb\n(optimised)" = "grey",
                               "2 over 3" = "white")) +
  facet_grid(stat ~ fhist, space = "free_x", switch = "y",
             scales = "free_y", labeller = "label_parsed") +
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
p_stats_23_ICV <- stats_plot %>%
  filter(stat %in% c("ICV")) %>%
  ggplot(aes(x = stock, y = value, fill = catch_rule)) +
  geom_hline(yintercept = 0, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = position_dodge2(preserve = "single"), width = 0.8, 
           show.legend = FALSE, colour = "black", size = 0.1) +
  scale_fill_manual("catch\nrule",
                    values = c("rfb\n(default)" = "black", 
                               "rfb\n(optimised)" = "grey",
                               "2 over 3" = "white")) +
  facet_grid(stat ~ fhist, space = "free_x", switch = "y",
             scales = "free_y", labeller = "label_parsed") +
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
p_stats_23_fitness <- stats_plot %>%
  filter(stat %in% c("fitness~value")) %>%
  ggplot(aes(x = stock, y = value, fill = catch_rule)) +
  geom_hline(yintercept = 0, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = position_dodge2(preserve = "single"), width = 0.8, 
           colour = "black", size = 0.1) +
  scale_fill_manual("catch\nrule",
                    values = c("rfb\n(default)" = "black", 
                               "rfb\n(optimised)" = "grey",
                               "2 over 3" = "white")) +
  facet_grid(stat ~ fhist, space = "free_x", switch = "y",
             scales = "free_y", labeller = "label_parsed") +
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
        legend.key.width = unit(0.15, "lines"),
        legend.position = "right",
        legend.background = element_blank()) +
  scale_y_continuous(trans = trans_from(0), limits = c(NA, NA))
p_stats_23_combined <- plot_grid(
  plot_grid(p_stats_23_SSB, p_stats_23_F, p_stats_23_C,p_stats_23_risk, 
            p_stats_23_ICV, p_stats_23_fitness + theme(legend.position = "none"), 
            ncol = 1, align = "v", rel_heights = c(1.25, 1, 1, 1, 1, 1.5)),
  get_legend(p_stats_23_fitness),
  ncol = 2, rel_widths = c(1, 0.1))
ggsave(filename = "output/plots/2over3_stats.png", plot = p_stats_23_combined,
       width = 17, height = 12, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/2over3_stats.pdf", plot = p_stats_23_combined,
       width = 17, height = 12, units = "cm")



### ------------------------------------------------------------------------ ###
### ICES default 2 over 3 vs rfb with and without constraint ####
### ------------------------------------------------------------------------ ###
stats_plot <- readRDS("output/plots/data_2over3_stats.rds")
stats_plot <- stats_plot %>%
  filter(!catch_rule %in% "rfb\n(optimised)")

### individual plots
p_stats_23_SSB <- stats_plot %>%
  filter(stat %in% c("SSB/B[MSY]")) %>%
  ggplot(aes(x = stock, y = value, fill = catch_rule)) +
  geom_hline(yintercept = 1, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = position_dodge2(preserve = "single"), width = 0.8, 
           show.legend = FALSE, colour = "black", size = 0.1) +
  scale_fill_manual("catch\nrule",
                    values = c("rfb\n(default)" = "black", 
                               "rfb\n(constrained)" = "grey",
                               "2 over 3" = "white")) +
  facet_grid(stat ~ fhist, space = "free_x", switch = "y",
             scales = "free_y", labeller = "label_parsed") +
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
p_stats_23_F <- stats_plot %>%
  filter(stat %in% c("F/F[MSY]")) %>%
  ggplot(aes(x = stock, y = value, fill = catch_rule)) +
  geom_hline(yintercept = 1, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = position_dodge2(preserve = "single"), width = 0.8, 
           show.legend = FALSE, colour = "black", size = 0.1) +
  scale_fill_manual("catch\nrule",
                    values = c("rfb\n(default)" = "black", 
                               "rfb\n(constrained)" = "grey",
                               "2 over 3" = "white")) +
  facet_grid(stat ~ fhist, space = "free_x", switch = "y",
             scales = "free_y", labeller = "label_parsed") +
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
p_stats_23_C <- stats_plot %>%
  filter(stat %in% c("Catch/MSY")) %>%
  ggplot(aes(x = stock, y = value, fill = catch_rule)) +
  geom_hline(yintercept = 1, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = position_dodge2(preserve = "single"), width = 0.8, 
           show.legend = FALSE, colour = "black", size = 0.1) +
  scale_fill_manual("catch\nrule",
                    values = c("rfb\n(default)" = "black", 
                               "rfb\n(constrained)" = "grey",
                               "2 over 3" = "white")) +
  facet_grid(stat ~ fhist, space = "free_x", switch = "y",
             scales = "free_y", labeller = "label_parsed") +
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
p_stats_23_risk <- stats_plot %>%
  filter(stat %in% c("B[lim]~risk")) %>%
  ggplot(aes(x = stock, y = value, fill = catch_rule)) +
  geom_hline(yintercept = 0, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = position_dodge2(preserve = "single"), width = 0.8, 
           show.legend = FALSE, colour = "black", size = 0.1) +
  scale_fill_manual("catch\nrule",
                    values = c("rfb\n(default)" = "black", 
                               "rfb\n(constrained)" = "grey",
                               "2 over 3" = "white")) +
  facet_grid(stat ~ fhist, space = "free_x", switch = "y",
             scales = "free_y", labeller = "label_parsed") +
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
p_stats_23_ICV <- stats_plot %>%
  filter(stat %in% c("ICV")) %>%
  ggplot(aes(x = stock, y = value, fill = catch_rule)) +
  geom_hline(yintercept = 0, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = position_dodge2(preserve = "single"), width = 0.8, 
           show.legend = FALSE, colour = "black", size = 0.1) +
  scale_fill_manual("catch\nrule",
                    values = c("rfb\n(default)" = "black", 
                               "rfb\n(constrained)" = "grey",
                               "2 over 3" = "white")) +
  facet_grid(stat ~ fhist, space = "free_x", switch = "y",
             scales = "free_y", labeller = "label_parsed") +
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
p_stats_23_fitness <- stats_plot %>%
  filter(stat %in% c("fitness~value")) %>%
  ggplot(aes(x = stock, y = value, fill = catch_rule)) +
  geom_hline(yintercept = 0, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = position_dodge2(preserve = "single"), width = 0.8, 
           colour = "black", size = 0.1) +
  scale_fill_manual("catch\nrule",
                    values = c("rfb\n(default)" = "black", 
                               "rfb\n(constrained)" = "grey",
                               "2 over 3" = "white")) +
  facet_grid(stat ~ fhist, space = "free_x", switch = "y",
             scales = "free_y", labeller = "label_parsed") +
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
        legend.key.width = unit(0.15, "lines"),
        legend.position = "right",
        legend.background = element_blank()) +
  scale_y_continuous(trans = trans_from(0), limits = c(NA, NA))
p_stats_23_combined <- plot_grid(
  plot_grid(p_stats_23_SSB, p_stats_23_F, p_stats_23_C,p_stats_23_risk, 
            p_stats_23_ICV, p_stats_23_fitness + theme(legend.position = "none"), 
            ncol = 1, align = "v", rel_heights = c(1.25, 1, 1, 1, 1, 1.5)),
  get_legend(p_stats_23_fitness),
  ncol = 2, rel_widths = c(1, 0.1))
ggsave(filename = "output/plots/2over3_stats_constrained.png", 
       plot = p_stats_23_combined,
       width = 17, height = 12, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/2over3_stats_constrained.pdf", 
       plot = p_stats_23_combined,
       width = 17, height = 12, units = "cm")

