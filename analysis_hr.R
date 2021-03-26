### ------------------------------------------------------------------------ ###
### analysis of constant harvest rate rule ####
### ------------------------------------------------------------------------ ###

library(mse)
library(tidyverse)
library(doParallel)
library(scales)
library(cowplot)
library(RColorBrewer)
library(akima)
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
### ... ####
### ------------------------------------------------------------------------ ###


### "correct" time series
### once a replicate collapsed, keep it at zero
#trgts <- c("Fmsy", "LFeM", "uniform", "length")
trgts <- "length"
res <- foreach(trgt = trgts, .errorhandling = "pass") %:%
  foreach(stock = stocks_subset, .errorhandling = "pass") %dopar% {
  tmp <- readRDS(paste0("output/10000_100/test/random/", stock, 
                        "/", trgt, "_1_", 
                        ifelse(trgt == "length", "TRUE", "FALSE"),
                        "_1_1_1_Inf_0.rds"))
  stk <- tmp@stock
  rm(tmp); gc()
  ### stock metrics
  SSBs <- FLCore::ssb(stk)
  Fs <- FLCore::fbar(stk)
  Cs <- FLCore::catch(stk)
  yrs_check <- ac(100:200)
  yrs <- dim(stk[, yrs_check])[2]
  its <- dim(stk[, yrs_check])[6]
  ### collapse correction
  if (isTRUE(TRUE)) {
    ### find collapses
    cd <- sapply(seq(its), function(x) {
      min_yr <- min(which(SSBs[, yrs_check,,,, x] < 1))
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
    SSBs[, yrs_check]@.Data[cd] <- 0
    Cs[, yrs_check]@.Data[cd] <- 0
    Fs[, yrs_check]@.Data[cd] <- 0
  }
  brp <- brps[[stock]]
  Bmsy <- c(refpts(brp)["msy", "ssb"])
  Fmsy <- c(refpts(brp)["msy", "harvest"])
  Cmsy <- c(refpts(brp)["msy", "yield"])
  ### summarise
  qnts <- FLQuants(SSB = SSBs/Bmsy, F = Fs/Fmsy, Catch = Cs/Cmsy)
  rm(stk); gc()
  return(qnts)
}
names(res) <- trgts
res <- lapply(res, function(x) {
  names(x) <- stocks_subset
  return(x)
})
saveRDS(res, paste0("output/10000_100/test/", paste0(trgts, collapse = "_")
                    , "_corrected_21-29.rds"))
# quit(save = "no")
# exit
# exit


### ------------------------------------------------------------------------ ###
### harvest rate vs. depletion (from uniform HR testing) ####
### ------------------------------------------------------------------------ ###
res <- readRDS("output/10000_100/test/corrected_21-29.rds")
### harvest rates
set.seed(33)
hr_rates <- runif(n = 10000, min = 0, max = 1)

### uniform: extract stats
stats <- lapply(stocks_subset, function(x) {
  stats_yrs <- seq(from = 10, to = 100, by = 10)
  names(stats_yrs) <- stats_yrs
  stats_tmp <- lapply(stats_yrs, function(y) {
    data.frame(
      SSB = c(apply(res$uniform[[x]]$SSB[, ac(seq(from = 101, to = 100 + y))], 6, median)),
      catch = c(apply(res$uniform[[x]]$Catch[, ac(seq(from = 101, to = 100 + y))], 6,
                      median)),
      years = y,
      HR = hr_rates,
      status = c(res$uniform[[x]]$SSB[, ac(100)]),
      stock = x
    )
  })
  stats_tmp <- do.call(rbind, stats_tmp)
})
stats <- do.call(rbind, stats)
saveRDS(stats, "output/10000_100/test/corrected_21-29_stats.rds")
stats <- readRDS("output/10000_100/test/corrected_21-29_stats.rds")

### summarise into discrete cells
cut_mean <- function(x, bins = NULL, interval = NULL) {
  if (!is.null(bins)) {
    tmp <- cut(x, bins)
  } else if (!is.null(interval)) {
    bins_tmp <- seq(from = -1e-10, to = max(x) + interval, by = interval)
    tmp <- cut(x, bins_tmp)
    
  }
  tmp <- as.character(tmp)
  tmp <- gsub(x = tmp, pattern = "\\(|\\]", replacement = "")
  tmp <- strsplit(x = tmp, split = ",")
  tmp <- sapply(tmp, as.numeric, simplify = FALSE)
  tmp <- sapply(tmp, median)
  return(tmp)
  
}
stats2 <- stats %>%
  group_by(stock, years) %>%
  mutate(HR = cut_mean(HR, interval = 0.05),
         status = cut_mean(status, interval = 0.05)) %>%
  group_by(stock, years, HR, status) %>%
  summarise(SSB = median(SSB),
            catch = median(catch)) %>%
  mutate(period = paste0(years, "~years"))
stats2$period <- factor(stats2$period, levels = unique(stats2$period))
### add k
stats3 <- stats2 %>% 
  left_join(stocks[, c("stock", "k")]) %>%
  mutate(stock_k = paste0(stock, "~~(italic(k)==", k, ")"))
stats3$stock_k <- factor(stats3$stock_k, levels = unique(stats3$stock_k))

saveRDS(stats3, "output/10000_100/test/corrected_21-29_stats_plot.rds")
stats3 <- readRDS("output/10000_100/test/corrected_21-29_stats_plot.rds")

# df_brill <- stats3 %>%
#   filter(stock == "bll")
max_y <- max(stats3$catch, na.rm = TRUE)
stats3 %>%
  ggplot(aes(x = status, y = HR, fill = catch)) +
  geom_raster() + 
  # geom_contour(aes(z = catch), breaks = c(0.75), colour = "black",
  #              alpha = 0.5, size = 0.1) +
  facet_grid(stock_k ~ period, labeller = "label_parsed") +
  scale_fill_gradientn(expression(catch/MSY), 
                       colours = c("red", "orange", "yellow", "green",
                                   "darkgreen"),
                       values = c(0, 0.33, 0.66, 1, max_y)/max_y
                       ) +
  theme_bw(base_size = 7) +
  labs(x = expression(initial~stock~status~(SSB/italic(B)[MSY])), 
       y = "harvest rate") +
  #xlim(c(0, 5)) +
  coord_cartesian(xlim = c(0, 5), y = c(0, 1)) +
  theme(panel.spacing.x = unit(0, units = "pt"),
        panel.spacing.y = unit(0, units = "pt"))

ggsave(filename = "output/plots/HR_catch.pdf",
       width = 25, height = 15, units = "cm", dpi = 600)
ggsave(filename = "output/plots/HR_catch2.png", type = "cairo",
       width = 25, height = 15, units = "cm", dpi = 3840/(25/2.54))

### only 10 & 100 years
stats3 %>%
  filter(period %in% c("10~years", "100~years")) %>%
  ggplot(aes(x = status, y = HR, fill = catch)) +
  geom_raster() + 
  # geom_contour(aes(z = catch), breaks = c(0.75), colour = "black",
  #              alpha = 0.5, size = 0.1) +
  facet_grid(period ~ stock_k, labeller = "label_parsed") +
  scale_fill_gradientn(expression(catch/MSY), 
                       colours = c("red", "orange", "yellow", "green",
                                   "darkgreen"),
                       values = c(0, 0.33, 0.66, 1, max_y)/max_y
                       ) +
  theme_bw(base_size = 7) +
  labs(x = expression(initial~stock~status~(SSB/italic(B)[MSY])), 
       y = "harvest rate") +
  #xlim(c(0, 5)) +
  coord_cartesian(xlim = c(0, 5), y = c(0, 1)) +
  theme(panel.spacing.x = unit(0, units = "pt"),
        panel.spacing.y = unit(0, units = "pt"))

ggsave(filename = "output/plots/HR_catch_10_100.pdf",
       width = 17, height = 5, units = "cm", dpi = 600)
ggsave(filename = "output/plots/HR_catch_10_100.png", type = "cairo",
       width = 17, height = 5, units = "cm", dpi = 300)


### ------------------------------------------------------------------------ ###
### compare stats vs. harvest rate ####
### ------------------------------------------------------------------------ ###

res <- readRDS("output/10000_100/test/corrected_21-29.rds")
### harvest rates
# set.seed(33)
# hr_rates <- runif(n = 10000, min = 0, max = 1)
n_blocks <- 50
hr_blocks <- cut(hr_rates, n_blocks)
hr_blocks_val <- seq(from = 0, to = 1 - 1/n_blocks, length.out = n_blocks) + ((1/n_blocks)/2)

### uniform: extract stats
stats_hr <- lapply(stocks_subset, function(stock) {#browser()
  stats_tmp <- lapply(seq_along(hr_blocks_val), function(block) {#browser()
    sub_i <- which(hr_blocks == sort(unique(hr_blocks))[[block]])
    data.frame(
      HR = hr_blocks_val[block],
      SSB_rel = c(apply(res$uniform[[stock]]$SSB[, ac(101:200),,,, sub_i], 1, median)),
      Catch_rel = c(apply(res$uniform[[stock]]$Catch[, ac(101:200),,,, sub_i], 1, median)),
      F_rel = c(apply(res$uniform[[stock]]$Catch[, ac(101:200),,,, sub_i], 1, median)),
      ICV = iav(res$uniform[[stock]]$Catch[, ac(100:200),,,, sub_i], from = 100, period = 1,
                  summary_all = median),
      risk_Blim = mean(c(res$uniform[[stock]]$SSB[, ac(100:200),,,, sub_i] < 
                     brps[[stock]]@Blim/c(brps[[stock]]@refpts["msy", "ssb"])), 
                     na.rm = TRUE),
      stock = stock
    )
  })
  stats_tmp <- do.call(rbind, stats_tmp)
})
stats_hr <- do.call(rbind, stats_hr)

### format for plotting
stats_hr <- stats_hr %>%
  left_join(stocks[, c("stock", "k")]) %>%
  mutate(stock_k = paste0(stock, "~~(italic(k)==", k, ")")) %>%
  mutate(stock_k = factor(stock_k, levels = unique(stock_k))) %>%
  pivot_longer(c(SSB_rel, Catch_rel, F_rel, ICV, risk_Blim)) %>%
  mutate(name = factor(name, 
                       levels = c("SSB_rel", "Catch_rel", "F_rel", "ICV", 
                                  "risk_Blim"), 
                       labels = c("SSB/italic(B)[MSY]", "Catch/MSY",
                                  "F/italic(F)[MSY]", "ICV", 
                                  "italic(B)[lim]~risk")))

### all stats
stats_hr %>% 
  ggplot(aes(x = HR, y = value)) +
  geom_smooth(method = "loess", span = 0.2, n = 100, colour = "grey", 
              level = 0, size = 0.3) +
  geom_point(size = 0.5, stroke = 0) +
  facet_grid(name ~ stock_k, labeller = "label_parsed", switch = "y",
             scales = "free_y") +
  theme_bw(base_size = 8) +
  theme(strip.placement.y = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8)) +
  labs(x = "harvest rate", y = "") +
  scale_x_continuous(breaks = c(0, 0.5, 1),
                     expand = expansion(mult = c(0.1, 0.1)))
  # scale_y_continuous(#expand = expansion(mult = c(0, 0), add = 0),
  #                    #limits = c(0, 1)
  #                    breaks = c(0, 0.25, 0.5, 0.75, 1))
ggsave(filename = "output/plots/HR_stats.pdf",
       width = 17, height = 10, units = "cm", dpi = 600)
ggsave(filename = "output/plots/HR_stats.png", type = "cairo",
       width = 17, height = 10, units = "cm", dpi = 3840/(17/2.54))

### catch peak vs. k
stats_hr_k <- stats_hr %>%
  filter(name == "Catch/MSY") %>%
  group_by(stock) %>%
  filter(value == max(value))
cor.test(x = stats_hr_k$HR, y = stats_hr_k$value)
stats_hr_k %>% 
  ggplot(aes(x = k, y = value)) +
  geom_smooth(method = "lm", formula = y~x) +
  geom_point(size = 0.5, stroke = 0) +
  geom_text(aes(label = stock), nudge_x = 0.03, size = 2) +
  theme_bw(base_size = 8) +
  labs(x = expression(italic(k)), y = "optimum harvest rate") +
  xlim(0, NA) + ylim(0, NA)
ggsave(filename = "output/plots/HR_vs_k.pdf",
       width = 8.5, height = 6, units = "cm", dpi = 600)
ggsave(filename = "output/plots/HR_vs_k.png", type = "cairo",
       width = 8.5, height = 6, units = "cm", dpi = 3840/(8.5/2.54))

### ------------------------------------------------------------------------ ###
### stats of uniform, Fmsy and LFeM proxy HRs ####
### ------------------------------------------------------------------------ ###
res <- readRDS("output/10000_100/test/corrected_21-29.rds")
stats <- foreach(trgt = trgts) %:% foreach(stock = stocks_subset) %do% {
  
  tmp <- res[[trgt]][[stock]]
  refpts <- brps[[stock]]@refpts
  
  list(
  SSB_rel = median(c(tmp$SSB[, ac(101:150)]), na.rm = TRUE),
  F_rel = median(c(tmp$F[, ac(101:150)]), na.rm = TRUE),
  Catch_rel = median(c(tmp$Catch[, ac(101:150)]), na.rm = TRUE),
  ICV = iav(tmp$Catch[, ac(100:150)], from = 100, period = 1,
                  summary_all = median),
  risk_Blim = mean(c(tmp$SSB[, ac(101:150)] < 
                       brps[[stock]]@Blim/c(refpts["msy", "ssb"])), 
                   na.rm = TRUE),
  risk_Bmsy = mean(c(tmp$SSB[, ac(101:150)] < 1), na.rm = TRUE),
  risk_halfBmsy = mean(c(tmp$SSB[, ac(101:150)]) < 0.5, na.rm = TRUE),
  risk_collapse = mean(c(tmp$SSB[, ac(101:150)] < 1/c(refpts["msy", "ssb"])), 
                       na.rm = TRUE),
  stock = stock, target = trgt
  )
  
}
stats <- do.call(rbind, lapply(stats, function(x) {do.call(rbind, x)}))



      list(
        risk_Blim = mean(c(SSBs < Blim), na.rm = TRUE),
        risk_Bmsy = mean(c(SSBs < Bmsy), na.rm = TRUE),
        risk_halfBmsy = mean(c(SSBs < Bmsy/2), na.rm = TRUE),
        risk_collapse = mean(c(SSBs < 1), na.rm = TRUE),
        SSB = median(c(SSBs), na.rm = TRUE), Fbar = median(c(Fs), na.rm = TRUE),
        Catch = median(c(Cs), na.rm = TRUE),
        SSB_rel = median(c(SSBs/Bmsy), na.rm = TRUE),
        Fbar_rel = median(c(Fs/Fmsy), na.rm = TRUE),
        Catch_rel = median(c(Cs/Cmsy), na.rm = TRUE),
        ICV = iav(Cs_long, from = 100, period = TAC_intvl,
                  summary_all = median)
      )





### ------------------------------------------------------------------------ ###
### catch rate in fishing history ####
### ------------------------------------------------------------------------ ###

inp_bll <- readRDS("input/10000_100/OM_2_mp_input/random/bll.rds")
plot(catch(inp_bll$om@stock), iter = 1:10)
plot(inp_bll$oem@observations$idx$idxB, iter = 1:10)

df <- FLQuants(catch = catch(inp_bll$om@stock)[, ac(50:100)],
               idx = inp_bll$oem@observations$idx$idxB[, ac(50:100)],
               rate = catch(inp_bll$om@stock)[, ac(50:100)] %=% NA_real_)
df$rate <- df$catch/df$idx
plot(df)

df <- as.data.frame(df)
df %>% 
  filter(iter %in% 1:5) %>%
  ggplot(aes(x = year, y = data, colour = iter)) +
  geom_line() +
  facet_wrap(~ qname, scales = "free_y", ncol = 1) +
  theme_bw()









### ------------------------------------------------------------------------ ###
### "length" scenario ####
### ------------------------------------------------------------------------ ###
stats <- lapply(stocks_subset, function(x) {
  readRDS(paste0("output/500_50/length/random/", x, 
                 "/collated_stats_length_0-1_1_1_1_1_Inf_0.rds"))
})
stats0 <- lapply(stocks_subset, function(x) {
  readRDS(paste0("output/500_50/length/random/", x, 
                 "/collated_stats_length_0-1_1_0_1_1_Inf_0.rds"))
})
stats <- do.call(rbind, stats)
stats0 <- do.call(rbind, stats0)
stats <- data.frame(lapply(stats, unlist))
stats0 <- data.frame(lapply(stats0, unlist))
stats <- rbind(stats, stats0)
### add k
stats <- stats %>%
  left_join(stocks[, c("stock", "k")]) %>%
  mutate(stock_k = paste0(stock, "~(italic(k)==", k, ")")) %>%
  mutate(stock_k = factor(stock_k, levels = unique(stock_k)))

stats %>% group_by(stock, idxB_lag) %>%
  filter(multiplier == 1) %>%
  print(n = Inf, width = Inf)

stats %>% group_by(stock, idxB_lag) %>%
  filter(Catch_rel == max(Catch_rel)) %>%
  print(n = Inf, width = Inf)

stats %>% group_by(stock) %>%
  filter(multiplier == 0.2) %>%
  print(n = Inf, width = Inf)

### plot catch vs. multiplier
### for idxB_lag = 1
stats %>%
  filter(idxB_lag == 1) %>%
  ggplot(aes(x = multiplier, y = Catch_rel)) +
  geom_line(size = 0.3) +
  facet_wrap(~ stock_k, labeller = "label_parsed") +
  theme_bw(base_size = 8) +
  labs(y = "Catch/MSY", x = "Multiplier")
ggsave(filename = "output/plots/HR_length_catch_vs_mult.pdf",
       width = 17, height = 8, units = "cm", dpi = 600)

### for idxB_lag = 0 & 1
stats %>%
  ggplot(aes(x = multiplier, y = Catch_rel, 
             colour = as.factor(idxB_lag), linetype = as.factor(idxB_lag))) +
  geom_line(size = 0.3) +
  scale_linetype_discrete("idx lag\n[years]") +
  scale_colour_discrete("idx lag\n[years]") +
  facet_wrap(~ stock_k, labeller = "label_parsed") +
  theme_bw(base_size = 8) +
  labs(y = "Catch/MSY", x = "Multiplier")
ggsave(filename = "output/plots/HR_length_catch_vs_mult_idxB_lag.pdf",
       width = 17, height = 8, units = "cm", dpi = 600)

### full stats vs. multiplier
stats_full <- stats %>%
  mutate(ICV = ifelse(multiplier == 0, NA, ICV)) %>%
  pivot_longer(c(SSB_rel, Catch_rel, Fbar_rel, ICV, risk_Blim)) %>%
    mutate(name = factor(name, 
                         levels = c("SSB_rel", "Catch_rel", "Fbar_rel", "ICV", 
                                    "risk_Blim"), 
                         labels = c("SSB/B[MSY]", "Catch/MSY",
                                    "F/F[MSY]", "ICV", 
                                    "B[lim]~risk")))
p <- stats_full %>% 
  ggplot(aes(x = multiplier, y = value,
             colour = as.factor(idxB_lag), linetype = as.factor(idxB_lag))) +
  geom_line(size = 0.3) +
  facet_grid(name ~ stock_k, labeller = "label_parsed", switch = "y",
             scales = "free_y") +
  scale_linetype_discrete("idx lag\n[years]") +
  scale_colour_discrete("idx lag\n[years]") +
  theme_bw(base_size = 8) +
  theme(strip.placement.y = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 6)) +
  labs(x = "multiplier", y = "") +
  ylim(c(0, NA)) +
  scale_x_continuous(breaks = c(0, 0.5, 1),
                     expand = expansion(mult = c(0.1, 0.1)))
p
ggsave(filename = "output/plots/HR_length_stats_1.pdf",
       width = 17, height = 10, units = "cm", dpi = 600)
ggsave(filename = "output/plots/HR_length_stats_1.png", type = "cairo",
       width = 17, height = 10, units = "cm", dpi = 3840/(17/2.54))

### multiplier = 1 only
stats_full %>% 
  ggplot(aes(x = multiplier, y = value,
             colour = as.factor(idxB_lag), linetype = as.factor(idxB_lag))) +
  # geom_line(size = 0.3, alpha = 0) +
  # geom_point(data = stats_full %>% filter(multiplier == 1),
  #            size = 1, stroke = 0) +
  geom_col(data = stats_full %>% filter(multiplier == 1),
           aes(fill = as.factor(idxB_lag)), 
           position = position_dodge(), 
           linetype = "solid") +
  scale_shape_discrete("idx lag\n[years]") +
  facet_grid(name ~ stock_k, labeller = "label_parsed", switch = "y",
             scales = "free_y") +
  scale_linetype_discrete("idx lag\n[years]") +
  scale_colour_discrete("idx lag\n[years]") +
  scale_fill_discrete("idx lag\n[years]") +
  theme_bw(base_size = 8) +
  theme(strip.placement.y = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 6)) +
  labs(x = "multiplier", y = "") +
  ylim(c(0, NA)) +
  scale_x_continuous(breaks = c(1),
                     expand = expansion(mult = c(0.1, 0.1)))
ggsave(filename = "output/plots/HR_length_stats_0.pdf",
       width = 17, height = 10, units = "cm", dpi = 600)
ggsave(filename = "output/plots/HR_length_stats_0.png", type = "cairo",
       width = 17, height = 10, units = "cm", dpi = 3840/(17/2.54))

### add vertical lines
### catch optimum by stock
### (first attempt, vertical lines, too distracting)
# p + geom_vline(data = stats_full %>%
#                  group_by(stock, idxB_lag) %>%
#                  filter(name == "Catch/MSY") %>%
#                  filter(value == max(value)) %>%
#                  select(-name),
#                aes(xintercept = multiplier, linetype = as.factor(idxB_lag)), 
#                show.legend = FALSE)
### use points
stats_tmp <- stats_full %>%
  group_by(stock, idxB_lag) %>%
  filter(name == "Catch/MSY") %>%
  filter(value == max(value)) %>%
  select(-risk_Bmsy, -risk_halfBmsy, -risk_collapse, -SSB, -Fbar, -Catch,
         -name, -value) %>%
  mutate(selection = "max_catch") %>%
  full_join(stats_full)
p + geom_point(data = stats_tmp %>% 
                 filter(selection == "max_catch"),
               aes(shape = as.factor(idxB_lag)),
               size = 1, stroke = 0) +
  scale_shape_discrete("idx lag\n[years]")
ggsave(filename = "output/plots/HR_length_stats_2_max_by_stock.pdf",
       width = 17, height = 10, units = "cm", dpi = 600)
ggsave(filename = "output/plots/HR_length_stats_2_max_by_stock.png", 
       type = "cairo",
       width = 17, height = 10, units = "cm", dpi = 3840/(17/2.54))

### lowest multiplier for catch optimum
stats_tmp <- stats_full %>%
  group_by(stock, idxB_lag) %>%
  filter(name == "Catch/MSY") %>%
  filter(value == max(value)) %>%
  ungroup() %>% group_by(idxB_lag) %>%
  filter(multiplier == min(multiplier)) %>%
  select(-risk_Bmsy, -risk_halfBmsy, -risk_collapse, -SSB, -Fbar, -Catch,
         -name, -value, -stock, -stock_k, -k) %>%
  mutate(selection = "lowest_multilier") %>%
  full_join(stats_full)
p + geom_point(data = stats_tmp %>% 
                 filter(selection == "lowest_multilier"),
               aes(shape = as.factor(idxB_lag)),
               size = 1, stroke = 0) +
  scale_shape_discrete("idx lag\n[years]")
ggsave(filename = "output/plots/HR_length_stats_3_min_mult_for_max_by_stock.pdf",
       width = 17, height = 10, units = "cm", dpi = 600)
ggsave(filename = "output/plots/HR_length_stats_3_min_mult_for_max_by_stock.png", 
       type = "cairo",
       width = 17, height = 10, units = "cm", dpi = 3840/(17/2.54))

### risk below 5% by stock
stats_tmp <- stats_full %>%
  group_by(stock, idxB_lag) %>%
  filter(name == "B[lim]~risk") %>%
  filter(value <= 0.05) %>%
  filter(value == max(value)) %>%
  select(-risk_Bmsy, -risk_halfBmsy, -risk_collapse, -SSB, -Fbar, -Catch,
         -name, -value) %>%
  mutate(selection = "5%") %>%
  full_join(stats_full)
p + geom_point(data = stats_tmp %>% 
                 filter(selection == "5%"),
               aes(shape = as.factor(idxB_lag)),
               size = 1, stroke = 0) +
  scale_shape_discrete("idx lag\n[years]")
ggsave(filename = "output/plots/HR_length_stats_4_risk_by_stock.pdf",
       width = 17, height = 10, units = "cm", dpi = 600)
ggsave(filename = "output/plots/HR_length_stats_4_risk_by_stock.png", 
       type = "cairo",
       width = 17, height = 10, units = "cm", dpi = 3840/(17/2.54))

### all stock below 5%
stats_tmp <- stats_full %>%
  group_by(stock, idxB_lag) %>%
  filter(name == "B[lim]~risk") %>%
  filter(value <= 0.05) %>%
  filter(value == max(value)) %>%
  ungroup() %>% group_by(idxB_lag) %>%
  filter(multiplier == min(multiplier)) %>%
  select(-risk_Bmsy, -risk_halfBmsy, -risk_collapse, -SSB, -Fbar, -Catch,
         -name, -value, -stock, -stock_k, -k) %>%
  mutate(selection = "lowest_multilier_5%") %>%
  full_join(stats_full)
p + geom_point(data = stats_tmp %>% 
                 filter(selection == "lowest_multilier_5%"),
               aes(shape = as.factor(idxB_lag)),
               size = 1, stroke = 0) +
  scale_shape_discrete("idx lag\n[years]")
ggsave(filename = "output/plots/HR_length_stats_5_risk_all_stocks.pdf",
       width = 17, height = 10, units = "cm", dpi = 600)
ggsave(filename = "output/plots/HR_length_stats_5_risk_all_stocks.png", 
       type = "cairo",
       width = 17, height = 10, units = "cm", dpi = 3840/(17/2.54))



### plot scenarios
### lag=1 mult = 1

def_with_lag <- stats %>% group_by(stock, idxB_lag) %>%
  filter(multiplier == 1 & idxB_lag == 1) %>%
  mutate(scenario = "default (with lag)")
def_no_lag <- stats %>% group_by(stock, idxB_lag) %>%
  filter(multiplier == 1 & idxB_lag == 0) %>%
  mutate(scenario = "default (no lag)")
max_catch_lag <- stats %>% group_by(stock, idxB_lag) %>%
  filter(Catch_rel == max(Catch_rel) & idxB_lag == 1) %>%
  mutate(scenario = "max catch (with lag)")
max_catch_no_lag <- stats %>% group_by(stock, idxB_lag) %>%
  filter(Catch_rel == max(Catch_rel) & idxB_lag == 0) %>%
  mutate(scenario = "max catch (no lag)")
min_mult_max_catch_lag <- stats %>% group_by(stock, idxB_lag) %>%
  filter(multiplier == min(max_catch_lag$multiplier) & idxB_lag == 1) %>%
  mutate(scenario = "min mult for max catch (with lag)")
min_mult_max_catch_no_lag <- stats %>% group_by(stock, idxB_lag) %>%
  filter(multiplier == min(max_catch_no_lag$multiplier) & idxB_lag == 1) %>%
  mutate(scenario = "min mult for max catch (no lag)")
stats_df <- do.call(rbind, list(def_with_lag, def_no_lag, max_catch_lag,
                                max_catch_no_lag, min_mult_max_catch_lag,
                                min_mult_max_catch_no_lag))



trans_from <- function(from = 1) {
  trans <- function(x) x - from
  inv <- function(x) x + from
  trans_new("from", trans, inv, 
            domain = c(from, Inf))
}

p_theme <- theme_bw(base_size = 8, base_family = "sans") +
  theme(panel.spacing.x = unit(0, units = "cm"),
        strip.placement.y = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        plot.margin = unit(x = c(1, 3, 0, 3), units = "pt"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())
  
p_SSB <- stats_df %>% 
  ggplot(aes(x = scenario, y = SSB_rel, fill = scenario,
             colour = scenario)) +
  geom_hline(yintercept = 1, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = position_dodge2(preserve = "single"), width = 0.8, 
           show.legend = FALSE, colour = "black", size = 0.1) +
  facet_grid("SSB/B[MSY]" ~ stock_k, scales = "free", space = "free_x", switch = "y",
             labeller = "label_parsed") +
  labs(y = "", x = "fitness function") +
  p_theme +
  scale_y_continuous(trans = trans_from(), limits = c(0, NA))
p_F <- stats_df %>% 
  ggplot(aes(x = scenario, y = Fbar_rel, fill = scenario,
             colour = scenario)) +
  geom_hline(yintercept = 1, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = position_dodge2(preserve = "single"), width = 0.8, 
           show.legend = FALSE, colour = "black", size = 0.1) +
  facet_grid("F/F[MSY]" ~ stock_k, scales = "free", space = "free_x", switch = "y",
             labeller = "label_parsed") +
  labs(y = "", x = "fitness function") +
  p_theme +
  theme(strip.text.x = element_blank(),
          plot.margin = unit(x = c(0, 3, 0, 3), units = "pt")) + 
  scale_y_continuous(trans = trans_from(), limits = c(0, 1.25))
p_Catch <- stats_df %>% 
  ggplot(aes(x = scenario, y = Catch_rel, fill = scenario,
             colour = scenario)) +
  geom_hline(yintercept = 1, linetype = "solid", size = 0.5, colour = "grey") +
  geom_col(position = position_dodge2(preserve = "single"), width = 0.8, 
           show.legend = FALSE, colour = "black", size = 0.1) +
  facet_grid("Catch/MSY" ~ stock_k, scales = "free", space = "free_x", switch = "y",
             labeller = "label_parsed") +
  labs(y = "", x = "fitness function") +
  p_theme +
  theme(strip.text.x = element_blank(),
          plot.margin = unit(x = c(0, 3, 0, 3), units = "pt")) + 
  scale_y_continuous(trans = trans_from(), limits = c(0, 1.25))
p_risk <- stats_df %>% 
    ggplot(aes(x = scenario, y = risk_Blim, fill = scenario,
               colour = scenario)) +
    geom_hline(yintercept = 0.05, 
               linetype = "solid", size = 0.5, 
               colour = "red") +
    geom_col(position = "dodge", show.legend = FALSE, width = 0.8, 
             colour = "black", size = 0.1) +
    facet_grid("B[lim]~risk" ~ stock_k, scales = "free", space = "free_x", switch = "y",
               labeller = "label_parsed") +
    labs(y = "", x = "fitness function") +
    p_theme +
    theme(strip.text.x = element_blank(),
          plot.margin = unit(x = c(0, 3, 0, 3), units = "pt")) + 
    scale_y_continuous(trans = trans_from(0), limits = c(0, 1))
p_ICV <- stats_df %>% 
    ggplot(aes(x = scenario, y = ICV, fill = scenario,
               colour = scenario)) +
    geom_hline(yintercept = 0, linetype = "solid", size = 0.5, colour = "grey") +
    geom_col(position = "dodge", show.legend = FALSE, width = 0.8, 
             colour = "black", size = 0.1) +
    facet_grid("ICV" ~ stock_k, scales = "free", space = "free_x", switch = "y",
               labeller = "label_parsed") +
    labs(y = "", x = "") +
    scale_y_continuous(trans = trans_from(0), limits = c(0, 1)) +
    theme_bw(base_size = 8, base_family = "sans") +
    theme(panel.spacing.x = unit(0, units = "cm"),
          strip.text.x = element_blank(),
          strip.placement.y = "outside",
          strip.background.y = element_blank(),
          strip.text.y = element_text(size = 8),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          plot.margin = unit(x = c(0, 3, 3, 3), units = "pt"))
p_all <- plot_grid(p_SSB, p_F, p_Catch, p_risk, p_ICV,
                   ncol = 1, align = "v",
                   rel_heights = c(1.25, 1, 1, 1, 3))
p_all
ggsave(filename = "output/plots/HR_length_stats.pdf",
       width = 17, height = 12, units = "cm", dpi = 600)


### ------------------------------------------------------------------------ ###
### visualise calculation of target HR ####
### ------------------------------------------------------------------------ ###

stock <- "gut"
input <- readRDS(paste0("input/500_50/OM_2_mp_input/random/", stock, ".rds"))
lhist <- stocks[stocks$stock == stock, ]

stk <- input$om@stock[, ac(50:100)]
idxL <- input$oem@observations$idx$idxL[, ac(50:100)]
idxL_dev <- input$oem@deviances$idx$idxL[, ac(50:100)]
idxL <- idxL * idxL_dev
Lc <- calc_lc(stk = stk[, ac(50:100)], a = lhist$a, b = lhist$b)
LFeM <- (lhist$linf + 2*1.5*c(Lc)) / (1 + 2*1.5)
idxB <- input$oem@observations$idx$idxB[, ac(50:100)]
idxB_dev <- input$oem@deviances$idx$idxB[, ac(50:100)]
idxB <- idxB * idxB_dev
i <- 13
yrs_above <- dimnames(idxL)$year[which((idxL/LFeM)[,,,,, i] >= 1)]

### mean length
as.data.frame(idxL[,,,,, i]) %>%
  ggplot(aes(x = year - 100, y = data)) +
  geom_line() +
  theme_bw(base_size = 8) +
  labs(x = "year", y = "mean catch length [cm]")
ggsave(filename = "output/plots/length_procedure_lmean_1.png",
       width = 10, height = 6, units = "cm", dpi = 600,
       type = "cairo")
### add LFeM
as.data.frame(idxL[,,,,, i]) %>%
  ggplot(aes(x = year - 100, y = data)) +
  geom_line() +
  geom_hline(yintercept = LFeM, colour = "black", linetype = "dashed") +
  theme_bw(base_size = 8) +
  labs(x = "year", y = "mean catch length [cm]")
ggsave(filename = "output/plots/length_procedure_lmean_2.png",
       width = 10, height = 6, units = "cm", dpi = 600,
       type = "cairo")
### mark length above LFeM
as.data.frame(idxL[,,,,, i]) %>%
  ggplot(aes(x = year - 100, y = data)) +
  geom_line() +
  geom_hline(yintercept = LFeM, colour = "black", linetype = "dashed") +
  geom_point(data = as.data.frame(idxL[,,,,, i]) %>% 
               filter(year %in% yrs_above), colour = "red") +
  theme_bw(base_size = 8) +
  labs(x = "year", y = "mean catch length [cm]")
ggsave(filename = "output/plots/length_procedure_lmean_3.png",
       width = 10, height = 6, units = "cm", dpi = 600,
       type = "cairo")
### catch
as.data.frame(catch(stk)[,,,,, i]) %>%
  ggplot(aes(x = year - 100, y = data)) +
  geom_line() +
  theme_bw(base_size = 8) +
  labs(x = "year", y = "catch")
ggsave(filename = "output/plots/length_procedure_catch.png",
       width = 10, height = 6, units = "cm", dpi = 600,
       type = "cairo")
### biomass index
as.data.frame(idxB[,,,,, i]) %>%
  ggplot(aes(x = year - 100, y = data)) +
  geom_line() +
  theme_bw(base_size = 8) +
  labs(x = "year", y = "biomass index")
ggsave(filename = "output/plots/length_procedure_idxB.png",
       width = 10, height = 6, units = "cm", dpi = 600,
       type = "cairo")
### catch/biomass index
as.data.frame(idxB[,,,,, i]/catch(stk)[,,,,, i]) %>%
  ggplot(aes(x = year - 100, y = data)) +
  geom_line() +
  theme_bw(base_size = 8) +
  labs(x = "year", y = "catch/biomass index")
ggsave(filename = "output/plots/length_procedure_cr_1.png",
       width = 10, height = 6, units = "cm", dpi = 600,
       type = "cairo")
### add points
as.data.frame(idxB[,,,,, i]/catch(stk)[,,,,, i]) %>%
  ggplot(aes(x = year - 100, y = data)) +
  geom_line() +
  geom_point(data = as.data.frame(idxB[,,,,, i]/catch(stk)[,,,,, i]) %>% 
               filter(year %in% yrs_above), colour = "red") +
  theme_bw(base_size = 8) +
  labs(x = "year", y = "catch/biomass index")
ggsave(filename = "output/plots/length_procedure_cr_2.png",
       width = 10, height = 6, units = "cm", dpi = 600,
       type = "cairo")
### add mean
as.data.frame(idxB[,,,,, i]/catch(stk)[,,,,, i]) %>%
  ggplot(aes(x = year - 100, y = data)) +
  geom_line() +
  geom_point(data = as.data.frame(idxB[,,,,, i]/catch(stk)[,,,,, i]) %>% 
               filter(year %in% yrs_above), colour = "red") +
  geom_hline(data = as.data.frame(idxB[,,,,, i]/catch(stk)[,,,,, i]) %>% 
               filter(year %in% yrs_above) %>%
               summarise(data = mean(data)),
             aes(yintercept = data), colour = "red") + 
  theme_bw(base_size = 8) +
  labs(x = "year", y = "catch/biomass index")
ggsave(filename = "output/plots/length_procedure_cr_3.png",
       width = 10, height = 6, units = "cm", dpi = 600,
       type = "cairo")



### ------------------------------------------------------------------------ ###
### HR with +20 -30% uncertainty cap - multipliers ####
### ------------------------------------------------------------------------ ###

res_cap <- foreach(stock = stocks$stock[1:29], .combine = bind_rows) %:%
  foreach(fhist = c("one-way", "random"), .combine = bind_rows) %do% {#browser()
    ### load data
    path <- paste0("output/500_50/length/", fhist, "/", stock, "/")
    path_runs <- paste0(path, "collated_stats_length_0-1_1_1_1_1_1.2_0.7.rds")
    if (!file.exists(path_runs)) return(NULL)
    print("found something")
    runs <- readRDS(path_runs)
    runs <- as.data.frame(lapply(runs, unlist))
    runs$fhist <- fhist
    return(runs)
}

res_cap <- res_cap %>%
  left_join(stocks[, c("stock", "k")]) %>%
  mutate(stock_k = paste0(stock, "~(italic(k)==", k, ")")) %>%
  mutate(stock_k = factor(stock_k, levels = unique(stock_k))) %>%
  mutate(stock = factor(stock, levels = stocks$stock))
  


saveRDS(res_cap, file = "output/500_50/length/cap2030.rds")
res_cap <- readRDS("output/500_50/length/cap2030.rds")


### plot
res_cap <- res_cap %>%
  #filter(stat_yrs == "all") %>%
  select(multiplier, risk_Blim, SSB_rel, Fbar_rel, Catch_rel, ICV,
         stock, fhist, stock_k) %>%
  pivot_longer(c(SSB_rel, Fbar_rel, Catch_rel, risk_Blim, ICV), 
               names_to = "key", values_to = "value") %>%
  mutate(stat = factor(key, levels = c("SSB_rel", "Fbar_rel", "Catch_rel",
                                       "risk_Blim", "ICV", "fitness"), 
                       labels = c("SSB/B[MSY]", "F/F[MSY]", "Catch/MSY", 
                                  "B[lim]~risk", "ICV", "fitness~value")))
stats_targets <- data.frame(stat = c("SSB/B[MSY]", "F/F[MSY]", "Catch/MSY", 
                                     "B[lim]~risk", "ICV"),
                            target = c(1, 1, 1, 0, 0))
### plot full period rfb-rule
p <- res_cap %>% 
  ggplot(aes(x = multiplier, y = value,
             colour = as.factor(fhist), linetype = as.factor(fhist))) +
  geom_line(size = 0.3) +
  geom_hline(data = data.frame(stat = "B[lim]~risk", y = 0.05),
             aes(yintercept = y), colour = "red") +
  facet_grid(stat ~ stock_k, labeller = "label_parsed", switch = "y",
             scales = "free_y") +
  scale_linetype_discrete("fishing\nhistory") +
  scale_colour_discrete("fishing\nhistory") +
  theme_bw(base_size = 8) +
  theme(strip.placement.y = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 6)) +
  labs(x = "multiplier", y = "") +
  ylim(c(0, NA)) +
  scale_x_continuous(breaks = c(0, 0.5, 1)#,
                     #expand = expansion(mult = c(0.1, 0.1))
                     )
p
ggsave(filename = "output/plots/length_cap2030.pdf",
       width = 50, height = 10, units = "cm")
ggsave(filename = "output/plots/length_cap2030.png", type = "cairo",
       width = 50, height = 10, units = "cm", dpi = 600)



### ------------------------------------------------------------------------ ###
### HR with +20 -30% uncertainty cap if index above trigger - multipliers ####
### ------------------------------------------------------------------------ ###

stocks_subset <- stocks$stock
res_cap <- foreach(stock = stocks_subset, .combine = bind_rows) %:%
  foreach(fhist = c("one-way", "random"), .combine = bind_rows) %do% {#browser()
    ### load data
    path <- paste0("output/500_100/length_cap_b/", fhist, "/", stock, "/")
    path_runs <- paste0(path, "collated_stats_length_0-1_1_1_1_1_1.2_0.7.rds")
    if (!file.exists(path_runs)) return(NULL)
    print("found something")
    runs <- readRDS(path_runs)
    runs <- as.data.frame(lapply(runs, unlist))
    runs$fhist <- fhist
    return(runs)
}

### split reporting periods
smry_stats <- c("risk_Blim", "risk_Bmsy", "risk_halfBmsy", "risk_collapse", 
                "SSB", "Fbar", "Catch", "SSB_rel", "Fbar_rel", "Catch_rel", "ICV")
periods <- c("first10", "41to50", "last10", "firsthalf", "lastfhalf", "11to50")
res_cap <- lapply(split(res_cap, seq(nrow(res_cap))), function(x) {
  stats_more <- t(sapply(periods, function(y) {
    as.numeric(x[grep(x = names(x), pattern = y)])
  }))
  stats <- rbind(all = as.numeric(x[smry_stats]), stats_more)
  stats <- as.data.frame(stats)
  colnames(stats) <- smry_stats
  stats$stat_yrs <- rownames(stats)
  rownames(stats) <- NULL
  ### add scenario definition
  stats <- cbind(x[c("stock", "multiplier", "comp_b", "idxB_lag", "idxB_range_3", 
                     "interval", "upper_constraint", "lower_constraint", 
                     "fhist")],
                 stats)
  return(stats)
})
res_cap <- do.call(rbind, res_cap)

### add k
res_cap <- res_cap %>%
  left_join(stocks[, c("stock", "k")]) %>%
  mutate(stock_k = paste0(stock, "~(italic(k)==", k, ")")) %>%
  mutate(stock_k = factor(stock_k, levels = unique(stock_k))) %>%
  mutate(stock = factor(stock, levels = stocks$stock))

saveRDS(res_cap, file = "output/500_100/length_cap_b/collated.rds")
res_cap <- readRDS("output/500_100/length_cap_b/collated.rds")


### prepare data for plotting
res_plot <- res_cap %>%
  filter(stat_yrs == "all") %>%
  select(multiplier, risk_Blim, SSB_rel, Fbar_rel, Catch_rel, ICV,
         stock, fhist, stock_k) %>%
  pivot_longer(c(SSB_rel, Fbar_rel, Catch_rel, risk_Blim, ICV), 
               names_to = "key", values_to = "value") %>%
  mutate(stat = factor(key, levels = c("SSB_rel", "Fbar_rel", "Catch_rel",
                                       "risk_Blim", "ICV", "fitness"), 
                       labels = c("SSB/B[MSY]", "F/F[MSY]", "Catch/MSY", 
                                  "B[lim]~risk", "ICV", "fitness~value")))
### plot summary stats for all stocks
p <- res_plot %>% 
  ggplot(aes(x = multiplier, y = value,
             colour = as.factor(fhist), linetype = as.factor(fhist))) +
  geom_line(size = 0.3) +
  geom_hline(data = data.frame(stat = "B[lim]~risk", y = 0.05),
             aes(yintercept = y), colour = "red") +
  facet_grid(stat ~ stock_k, labeller = "label_parsed", switch = "y",
             scales = "free_y") +
  scale_linetype_discrete("fishing\nhistory") +
  scale_colour_discrete("fishing\nhistory") +
  theme_bw(base_size = 8) +
  theme(strip.placement.y = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 6)) +
  labs(x = "multiplier", y = "") +
  ylim(c(0, NA)) +
  scale_x_continuous(breaks = c(0, 0.5, 1)#,
                     #expand = expansion(mult = c(0.1, 0.1))
                     )
p
ggsave(filename = "output/plots/length_cap2030_b/all_stocks_stats.pdf",
       width = 50, height = 10, units = "cm")
ggsave(filename = "output/plots/length_cap2030_b/all_stocks_stats.png", 
       type = "cairo", width = 50, height = 10, units = "cm", dpi = 600)

### medians
res_plot <- res_cap %>%
  filter(stat_yrs == "all") %>%
  mutate(ICV = ifelse(multiplier == 0 & ICV == 1, NA, ICV)) %>%
  filter(stock %in% c("gut", "whg", "bll", "lem", "ane")) %>%
  select(multiplier, risk_Blim, SSB_rel, Fbar_rel, Catch_rel, ICV,
         fhist, stock, stock_k, k)
res_mult <- res_plot %>% 
  group_by(multiplier) %>%
  summarise(risk_Blim = median(risk_Blim),
            SSB_rel = median(SSB_rel),
            Fbar_rel = median(Fbar_rel),
            Catch_rel = median(Catch_rel),
            ICV = median(ICV)) %>%
  mutate(stock = "median")
res_mult %>%
  filter(risk_Blim >= 0.04 & risk_Blim <= 0.07)
res_plot <- res_plot %>%
  full_join(res_mult)
res_plot <- res_plot %>%
  pivot_longer(c(SSB_rel, Fbar_rel, Catch_rel, risk_Blim, ICV), 
               names_to = "key", values_to = "value") %>%
  mutate(stat = factor(key, levels = c("SSB_rel", "Fbar_rel", "Catch_rel",
                                       "risk_Blim", "ICV", "fitness"), 
                       labels = c("SSB/B[MSY]", "F/F[MSY]", "Catch/MSY", 
                                  "B[lim]~risk", "ICV", "fitness~value"))) %>%
  mutate(group = ifelse(stock == "median", "median", NA),
         group = ifelse(stock != "median" & fhist == "one-way", "one-way", group),
         group = ifelse(stock != "median" & fhist == "random", "random", group))

res_plot %>% 
  ggplot(aes(x = multiplier, y = value,
             colour = group, 
             alpha = group,
             group = interaction(stock, fhist))) +
  geom_line(size = 0.3) +
  geom_hline(data = data.frame(stat = "B[lim]~risk", y = 0.05),
             aes(yintercept = y), colour = "red") +
  facet_wrap(~ stat, labeller = "label_parsed", switch = "y",
             scales = "free_y") +
  #scale_colour_manual("", values = c("TRUE" = "black", "FALSE" = "blue")) +
  scale_colour_manual("", values = c("median" = "black", "one-way" = "blue",
                                     "random" = "green")) +
  scale_alpha_manual("", values = c("median" = 1, "one-way" = 0.2,
                                     "random" = 0.2)) +
  theme_bw(base_size = 8) +
  theme(strip.placement.y = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 6)) +
  labs(x = "multiplier", y = "") +
  ylim(c(0, NA)) +
  scale_x_continuous(breaks = c(0, 0.5, 1))
ggsave(filename = "output/plots/length_cap2030_b/multiplier.pdf",
       width = 17, height = 10, units = "cm")
ggsave(filename = "output/plots/length_cap2030_b/multiplier.png", 
       type = "cairo", width = 17, height = 10, units = "cm", dpi = 600)



### ------------------------------------------------------------------------ ###
### pure HR - multipliers - 50 years ####
### ------------------------------------------------------------------------ ###

res_def <- foreach(stock = stocks$stock[1:29], .combine = bind_rows) %:%
  foreach(fhist = c("one-way", "random"), .combine = bind_rows) %do% {#browser()
    ### load data
    path <- paste0("output/500_50/length/", fhist, "/", stock, "/")
    path_runs <- paste0(path, "collated_stats_length_0-2_1_1_1_1_Inf_0.rds")
    if (!file.exists(path_runs)) return(NULL)
    print("found something")
    runs <- readRDS(path_runs)
    runs <- as.data.frame(lapply(runs, unlist))
    runs$fhist <- fhist
    return(runs)
}

res_def <- res_def %>%
  left_join(stocks[, c("stock", "k")]) %>%
  mutate(stock_k = paste0(stock, "~(italic(k)==", k, ")")) %>%
  mutate(stock_k = factor(stock_k, levels = unique(stock_k))) %>%
  mutate(stock = factor(stock, levels = stocks$stock))
  


saveRDS(res_def, file = "output/500_50/length/all_def_mult.rds")
res_def <- readRDS("output/500_50/length/all_def_mult.rds")


### plot
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
### plot all stocks
p <- res_def_p %>% 
  ggplot(aes(x = multiplier, y = value,
             colour = as.factor(fhist), linetype = as.factor(fhist))) +
  geom_line(size = 0.3) +
  # geom_hline(data = data.frame(stat = "B[lim]~risk", y = 0.05),
  #            aes(yintercept = y), colour = "red") +
  facet_grid(stat ~ stock_k, labeller = "label_parsed", switch = "y",
             scales = "free_y") +
  scale_linetype_discrete("fishing\nhistory") +
  scale_colour_discrete("fishing\nhistory") +
  theme_bw(base_size = 8) +
  theme(strip.placement.y = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 6)) +
  labs(x = "multiplier", y = "") +
  ylim(c(0, NA))# +
  # scale_x_continuous(breaks = c(0, 0.5, 1)#,
  #                    #expand = expansion(mult = c(0.1, 0.1))
  #                    )
p
ggsave(filename = "output/plots/length_all_def_mult.pdf",
       width = 50, height = 10, units = "cm")
ggsave(filename = "output/plots/length_all_def_mult.png", type = "cairo",
       width = 50, height = 10, units = "cm", dpi = 600)

### plot subset and stats individually
p_SSB <- res_def_p %>% 
  filter(stock %in% c("ang3", "pol", "bll", "san")) %>%
  filter(key == "SSB_rel") %>%
  ggplot(aes(x = multiplier, y = value,
             colour = as.factor(fhist), linetype = as.factor(fhist))) +
  geom_line(size = 0.3) +
  facet_grid(stat ~ stock_k, labeller = "label_parsed", switch = "y",
             scales = "free_y") +
  scale_linetype_discrete("fishing history") +
  scale_colour_brewer("fishing history", palette = "Dark2") +
  labs(x = "multiplier", y = "") +
  ylim(c(0, NA)) +
  theme_bw(base_size = 8) +
  theme(strip.placement.y = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 8),
        legend.key.height = unit(0.6, "lines"),
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
  geom_line(size = 0.3, show.legend = FALSE) +
  scale_colour_brewer("fishing history", palette = "Dark2") +
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
  geom_line(size = 0.3, show.legend = FALSE) +
  scale_colour_brewer("fishing history", palette = "Dark2") +
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
ggsave(filename = "output/plots/length_all_subset_def_mult.pdf",
       width = 17, height = 10, units = "cm")
ggsave(filename = "output/plots/length_all_subset_def_mult.png", type = "cairo",
       width = 17, height = 10, units = "cm", dpi = 600)

### find max catch
max_catch <- res_def %>%
  group_by(stock, fhist, k) %>%
  filter(Catch_rel == max(Catch_rel))
summary(max_catch)
table(max_catch$multiplier)

### estimate correlation
max_catch_cor <- foreach(fhist = c("one-way", "random")) %:%
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
#  	Pearson's product-moment correlation 
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
# Correlation for fhist: one-way; stat: multiplier
#   
#  	Pearson's product-moment correlation 
#   
#  data:  data_i$k and data_i$stat 
#  t = -10.358, df = 27, p-value = 6.636e-11 
#  alternative hypothesis: true correlation is not equal to 0 
#  95 percent confidence interval: 
#   -0.9493434 -0.7842543 
#  sample estimates: 
#         cor  
#  -0.8938402  
#   
# 
# Correlation for fhist: random; stat: Catch_rel
#   
#  	Pearson's product-moment correlation 
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
# Correlation for fhist: random; stat: multiplier
#   
#  	Pearson's product-moment correlation 
#   
#  data:  data_i$k and data_i$stat 
#  t = -8.5716, df = 27, p-value = 3.473e-09 
#  alternative hypothesis: true correlation is not equal to 0 
#  95 percent confidence interval: 
#   -0.9301314 -0.7116912 
#  sample estimates: 
#         cor  
#  -0.8551425


### prepare for plotting
max_catch_p <- max_catch %>%
  select(multiplier, Catch_rel, risk_Blim, stock, fhist, k) %>%
  pivot_longer(c(multiplier, Catch_rel, risk_Blim), 
               names_to = "key", values_to = "value") %>%
  mutate(stat = factor(key, levels = c("multiplier", "Catch_rel", "risk_Blim"), 
                       labels = c("multiplier~(x)", "catch/MSY", "B[lim]~risk")))
### plot
p_max_catch <- max_catch_p %>%
  filter(stat %in% c("multiplier~(x)", "catch/MSY")) %>%
  ggplot(aes(x = k, y = value, colour = fhist, shape = fhist, linetype = fhist)) +
  geom_smooth(method = lm, se = FALSE, size = 0.3) +
  geom_point(size = 0.4) +
  scale_colour_brewer("fishing history", palette = "Dark2") +
  scale_shape("fishing history") +
  scale_linetype("fishing history") +
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
        legend.key.height = unit(0.6, "lines"),
        legend.position = c(0.77, 0.92))
p_max_catch

### combine plots
plot_grid(p_stats, p_max_catch,
          ncol = 2, rel_widths = c(0.67, 0.33),
          labels = c("(a)", "(b)"), label_size = 10)
ggsave(filename = "output/plots/HR_mult_max_catch_cor.png", 
       type = "cairo",
       width = 17, height = 10, units = "cm", dpi = 600)
ggsave(filename = "output/plots/HR_mult_max_catch_cor.pdf",
       width = 17, height = 10, units = "cm")



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

ggplot() +
  annotate(geom = "segment", x = 1, xend = 1, y = 0, yend = 1,
           linetype = "dotted") +
  annotate(geom = "segment", x = 0, xend = 1, y = 1, yend = 1,
           linetype = "dotted") +
  geom_line(data = data.frame(x = seq(0, 1, length.out = 1000), 
                              y = seq(0, 1, length.out = 1000),
                              z = seq(0, 1, length.out = 1000)),
            aes(x = x, y = y, colour = z),
            size = 1, show.legend = FALSE) +
  scale_colour_gradient(low = "red", high = "yellow") +
  geom_line(data = data.frame(x = c(1, 2), 
                              y = c(1, 1)),
            aes(x = x, y = y), colour = "green3", size = 1) +
  scale_x_continuous("index I", expand = c(0, 0),
                     breaks = c(0, 1), labels = c(0, expression(I[trigger]))) +
  scale_y_continuous("target harvest rate", limits = c(0, 1.2), expand = c(0, 0),
                     breaks = c(0, 1), labels = c(0, "H")) +
  theme_classic()
ggsave(filename = "output/plots/HR_principle.png",
       width = 8.5, height = 5, units = "cm", dpi = 600,
       type = "cairo")
ggsave(filename = "output/plots/HR_principle.pdf",
       width = 8.5, height = 5, units = "cm")

data.frame(x = c(0, 1, 2),
           y = c(0, 1, 1)) %>%
  ggplot(aes(x = x, y = y)) +
  annotation_raster(matrix(colorRampPalette(c("red", "orange", "yellow"))(255),
                           nrow = 1),
                    xmin = 0, xmax = 1, ymin = 0, ymax = 1, interpolate = TRUE) +
  annotate("polygon", x = c(0, 1, 0), y = c(0, 1.01, 1.01), fill = "white") +
  annotate("polygon", x = c(1, 2, 2, 1), y = c(0, 0, 1, 1), fill = "green") +
  geom_line() +
  scale_x_continuous("index I", expand = c(0, 0),
                     breaks = c(0, 1), labels = c(0, expression(I[trigger]))) +
  scale_y_continuous("target harvest rate", limits = c(0, 1.2), expand = c(0, 0),
                     breaks = c(0, 1), labels = c(0, "H")) +
  annotate(geom = "segment", x = 1, xend = 1, y = 0, yend = 1,
           linetype = "dotted") +
  annotate(geom = "segment", x = 0, xend = 1, y = 1, yend = 1,
           linetype = "dotted") +
  #theme(panel.background = element_blank())
  theme_classic()
ggsave(filename = "output/plots/HR_principle_area.png",
       width = 8.5, height = 5, units = "cm", dpi = 600,
       type = "cairo")
ggsave(filename = "output/plots/HR_principle_area.pdf",
       width = 8.5, height = 5, units = "cm")


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


### ------------------------------------------------------------------------ ###
### sensitivity - stats over time ####
### ------------------------------------------------------------------------ ###

stats_time <- foreach(fhist = c("random", "one-way"),
                     .combine = rbind) %do% {
  
  res <- readRDS(paste0("output/hr/500_100/sensitivity/", fhist, 
                        "/pol/mp_length_1_TRUE_1_1_1_Inf_0.rds"))
  
  ### collapse correction
  res_corrected <- collapse_correction(stk = res@stock, yrs = 101:200)
  
  ### template
  df_i <- data.frame(year = 1:100)
  
  ### Blim risk
  df_i$risk_cumulative <- sapply(1:100, function(x) {
    mean(c(res_corrected$ssb[, ac(seq(from = 101, length.out = x))] < Blim), 
         na.rm = TRUE)
  })
  df_i$risk_annual <- c(apply(res_corrected$ssb < Blim, 2, mean, na.rm = TRUE))
  
  ### SSB
  df_i$SSB_annual <- sapply(1:100, function(x) {
    median(c(res_corrected$ssb[, x]/Bmsy), na.rm = TRUE)
  })
  df_i$SSB_cumulative <-  sapply(1:100, function(x) {
    median(c(res_corrected$ssb[, ac(seq(from = 101, length.out = x))]/Bmsy), 
           na.rm = TRUE)
  })
  
  ### Catch
  df_i$Catch_annual <- sapply(1:100, function(x) {
    median(c(res_corrected$catch[, x]/MSY), na.rm = TRUE)
  })
  df_i$Catch_cumulative <-  sapply(1:100, function(x) {
    median(c(res_corrected$catch[, ac(seq(from = 101, length.out = x))]/MSY), 
           na.rm = TRUE)
  })
  
  df_i <- df_i %>%
    pivot_longer(2:7, names_to = c(".value", "period"), names_sep = "_") %>%
    mutate(fhist = !!fhist)
  return(df_i)
}
saveRDS(stats_time, file = "output/hr_pol_sensitivity_time.rds")
write.csv(stats_time, row.names = FALSE, 
          file = "output/hr_pol_sensitivity_time.csv")
stats_time <- readRDS("output/hr_pol_sensitivity_time.rds")

p_time_ssb <- stats_time %>%
  ggplot(aes(x = year, y = SSB, colour = fhist, linetype = period)) +
  geom_line(show.legend = FALSE) +
  scale_colour_brewer(palette = "Dark2") +
  xlim(0, NA) + ylim(0, NA) +
  labs(x = "year", y = expression(SSB/B[MSY])) +
  theme_bw(base_size = 8)
p_time_catch <- stats_time %>%
  ggplot(aes(x = year, y = Catch, colour = fhist, linetype = period)) +
  geom_line(show.legend = FALSE) +
  scale_colour_brewer(palette = "Dark2") +
  xlim(0, NA) + ylim(0, NA) +
  labs(x = "year", y = expression(Catch/MSY)) +
  theme_bw(base_size = 8)
p_time_risk <- stats_time %>%
  ggplot(aes(x = year, y = risk, colour = fhist, linetype = period)) +
  geom_line(show.legend = TRUE) +
  scale_colour_brewer("fishing history", palette = "Dark2") +
  scale_linetype("") +
  xlim(0, NA) + ylim(0, NA) +
  labs(x = "year", y = expression(B[lim]~risk)) +
  theme_bw(base_size = 8) +
  theme(legend.position = c(0.7, 0.8),
        legend.key.height = unit(0.8, "line"),
        legend.key = element_blank(),
        legend.background = element_blank())

p_time <- plot_grid(p_time_ssb, p_time_catch, p_time_risk,
                    ncol = 3, align = "hv")
p_time

### ------------------------------------------------------------------------ ###
### sensitivity - stock status ####
### ------------------------------------------------------------------------ ###

stats_status <- foreach(fhist = c("random", "one-way"),
                     .combine = rbind) %do% {
  
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
  
  ### number of replicates per group:
  group_n <- sapply(SSB_levels, function(x) {
    length(which(SSB_groups %in% x))
  })
  group_n[sort(names(group_n))]
  
  ### Blim risk per group
  ### SSB is on absolut scale 
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
  
  df_i <- data.frame(SSB0 = unlist(SSB_levels)[-pos_remove],
                     SSB = unlist(SSB_group)[-pos_remove],
                     Catch = unlist(Catch_group)[-pos_remove],
                     risk = unlist(risk_group)[-pos_remove],
                     fhist = fhist,
                     iter = unlist(group_n)[-pos_remove])
  row.names(df_i) <- NULL
  
  return(df_i)
}
write.csv(stats_status, row.names = FALSE, 
          file = "output/hr_pol_sensitivity_status.csv")
saveRDS(stats_status, file = "output/hr_pol_sensitivity_status.rds")
stats_status <- readRDS("output/hr_pol_sensitivity_status.rds")


p_status_ssb <- stats_status %>%
  filter(iter >= 200) %>%
  ggplot(aes(x = SSB0, y = SSB, colour = fhist)) +
  geom_line(show.legend = FALSE) +
  scale_colour_brewer(palette = "Dark2") +
  xlim(0, 1.5) + ylim(0, NA) +
  labs(x = expression(SSB[y == 0]/B[MSY]), y = expression(SSB/B[MSY])) +
  theme_bw(base_size = 8)
p_status_catch <- stats_status %>%
  filter(iter >= 200) %>%
  ggplot(aes(x = SSB0, y = Catch, colour = fhist)) +
  geom_line(show.legend = FALSE) +
  scale_colour_brewer(palette = "Dark2") +
  xlim(0, 1.5) + ylim(0, NA) +
  labs(x = expression(SSB[y == 0]/B[MSY]), y = expression(Catch/MSY)) +
  theme_bw(base_size = 8)
p_status_risk <- stats_status %>%
  filter(iter >= 200) %>%
  ggplot(aes(x = SSB0, y = risk, colour = fhist)) +
  geom_line(show.legend = FALSE) +
  scale_colour_brewer(palette = "Dark2") +
  xlim(0, 1.5) + ylim(0, NA) +
  labs(x = expression(SSB[y == 0]/B[MSY]), y = expression(B[lim]~risk)) +
  theme_bw(base_size = 8)

p_status <- plot_grid(p_status_ssb, p_status_catch, p_status_risk,
                      ncol = 3, align = "hv")
p_status

### ------------------------------------------------------------------------ ###
### sensitivity - recruitment variability ####
### ------------------------------------------------------------------------ ###

stats_rec <- foreach(fhist = c("one-way", "random"), 
                     .combine = bind_rows) %do% {
  tmp <- readRDS(
    paste0("output/hr/500_50/sensitivity/", fhist, 
           "/pol/collated_stats_length_1_1_1_1_1_Inf_0_0.2_0.2_0-1_0.rds"))
  tmp <- as.data.frame(lapply(tmp, unlist))
  tmp$fhist <- fhist
  tmp$n_yrs <- 50
  tmp$n_iter <- 500
  return(tmp)
}
write.csv(stats_rec, row.names = FALSE, 
          file = "output/hr_pol_sensitivity_rec.csv")
saveRDS(stats_rec, file = "output/hr_pol_sensitivity_rec.rds")
stats_rec <- readRDS("output/hr_pol_sensitivity_rec.rds")

p_rec_ssb <- stats_rec %>%
  ggplot(aes(x = sigmaR, y = SSB_rel, colour = fhist)) +
  geom_line(show.legend = FALSE) +
  scale_colour_brewer(palette = "Dark2") +
  xlim(0, NA) + ylim(0, NA) +
  labs(x = expression(sigma[R]), y = expression(SSB/B[MSY])) +
  theme_bw(base_size = 8)
p_rec_catch <- stats_rec %>%
  ggplot(aes(x = sigmaR, y = Catch_rel, colour = fhist)) +
  geom_line(show.legend = FALSE) +
  scale_colour_brewer(palette = "Dark2") +
  xlim(0, NA) + ylim(0, NA) +
  labs(x = expression(sigma[R]), y = expression(Catch/MSY)) +
  theme_bw(base_size = 8)
p_rec_risk <- stats_rec %>%
  ggplot(aes(x = sigmaR, y = risk_Blim, colour = fhist)) +
  geom_line(show.legend = FALSE) +
  scale_colour_brewer(palette = "Dark2") +
  xlim(0, NA) + ylim(0, NA) +
  labs(x = expression(sigma[R]), y = expression(B[lim]~risk)) +
  theme_bw(base_size = 8)

p_rec <- plot_grid(p_rec_ssb, p_rec_catch, p_rec_risk,
                   ncol = 3, align = "hv")
p_rec

### ------------------------------------------------------------------------ ###
### sensitivity - biomass index ####
### ------------------------------------------------------------------------ ###

stats_idx <- foreach(fhist = c("one-way", "random"), 
                     .combine = bind_rows) %do% {
  tmp <- readRDS(
    paste0("output/hr/500_50/sensitivity/", fhist, 
           "/pol/collated_stats_length_1_1_1_1_1_Inf_0_0.2_0-1_0.6_0.rds"))
  tmp <- as.data.frame(lapply(tmp, unlist))
  tmp$fhist <- fhist
  tmp$n_yrs <- 50
  tmp$n_iter <- 500
  return(tmp)
}
write.csv(stats_idx, row.names = FALSE, 
          file = "output/hr_pol_sensitivity_idx.csv")
saveRDS(stats_idx, file = "output/hr_pol_sensitivity_idx.rds")
stats_idx <- readRDS("output/hr_pol_sensitivity_idx.rds")



p_idx_ssb <- stats_idx %>%
  ggplot(aes(x = sigmaB, y = SSB_rel, colour = fhist)) +
  geom_line(show.legend = FALSE) +
  scale_colour_brewer(palette = "Dark2") +
  xlim(0, NA) + ylim(0, NA) +
  labs(x = expression(sigma[B]), y = expression(SSB/B[MSY])) +
  theme_bw(base_size = 8)
p_idx_catch <- stats_idx %>%
  ggplot(aes(x = sigmaB, y = Catch_rel, colour = fhist)) +
  geom_line(show.legend = FALSE) +
  scale_colour_brewer(palette = "Dark2") +
  xlim(0, NA) + ylim(0, NA) +
  labs(x = expression(sigma[B]), y = expression(Catch/MSY)) +
  theme_bw(base_size = 8)
p_idx_risk <- stats_idx %>%
  ggplot(aes(x = sigmaB, y = risk_Blim, colour = fhist)) +
  geom_line(show.legend = FALSE) +
  scale_colour_brewer(palette = "Dark2") +
  xlim(0, NA) + ylim(0, NA) +
  labs(x = expression(sigma[B]), y = expression(B[lim]~risk)) +
  theme_bw(base_size = 8)

p_idx <- plot_grid(p_idx_ssb, p_idx_catch, p_idx_risk,
                   ncol = 3, align = "hv")
p_idx

### ------------------------------------------------------------------------ ###
### sensitivity - length index ####
### ------------------------------------------------------------------------ ###

stats_lngth <- foreach(fhist = c("one-way", "random"), 
                     .combine = bind_rows) %do% {
  tmp <- readRDS(
    paste0("output/hr/500_50/sensitivity/", fhist, 
           "/pol/collated_stats_length_1_1_1_1_1_Inf_0_0-1_0.2_0.6_0.rds"))
  tmp <- as.data.frame(lapply(tmp, unlist))
  tmp$fhist <- fhist
  tmp$n_yrs <- 50
  tmp$n_iter <- 500
  return(tmp)
}
write.csv(stats_lngth, row.names = FALSE, 
          file = "output/hr_pol_sensitivity_lngth.csv")
saveRDS(stats_lngth, file = "output/hr_pol_sensitivity_lngth.rds")
stats_lngth <- readRDS("output/hr_pol_sensitivity_lngth.rds")



p_lngth_ssb <- stats_lngth %>%
  ggplot(aes(x = sigmaL, y = SSB_rel, colour = fhist)) +
  geom_line(show.legend = FALSE) +
  scale_colour_brewer(palette = "Dark2") +
  xlim(0, NA) + ylim(0, NA) +
  labs(x = expression(sigma[L]), y = expression(SSB/B[MSY])) +
  theme_bw(base_size = 8)
p_lngth_catch <- stats_lngth %>%
  ggplot(aes(x = sigmaL, y = Catch_rel, colour = fhist)) +
  geom_line(show.legend = FALSE) +
  scale_colour_brewer(palette = "Dark2") +
  xlim(0, NA) + ylim(0, NA) +
  labs(x = expression(sigma[L]), y = expression(Catch/MSY)) +
  theme_bw(base_size = 8)
p_lngth_risk <- stats_lngth %>%
  ggplot(aes(x = sigmaL, y = risk_Blim, colour = fhist)) +
  geom_line(show.legend = FALSE) +
  scale_colour_brewer(palette = "Dark2") +
  xlim(0, NA) + ylim(0, NA) +
  labs(x = expression(sigma[L]), y = expression(B[lim]~risk)) +
  theme_bw(base_size = 8)

p_lngth <- plot_grid(p_lngth_ssb, p_lngth_catch, p_lngth_risk,
                   ncol = 3, align = "hv")
p_lngth

### ------------------------------------------------------------------------ ###
### sensitivity - length & biomass index ####
### ------------------------------------------------------------------------ ###

stats_idx_lngth <- foreach(fhist = c("one-way", "random"), 
                       .combine = bind_rows) %do% {
  tmp <- readRDS(
    paste0("output/hr/500_50/sensitivity/", fhist, 
           "/pol/collated_stats_length_1_1_1_1_1_Inf_0_0-1_0-1_0.6_0.rds"))
  tmp <- as.data.frame(lapply(tmp, unlist))
  tmp$fhist <- fhist
  tmp$n_yrs <- 50
  tmp$n_iter <- 500
  return(tmp)
}
write.csv(stats_idx_lngth, row.names = FALSE, 
          file = "output/hr_pol_sensitivity_lngth.csv")
saveRDS(stats_idx_lngth, file = "output/hr_pol_sensitivity_idx_lngth.rds")
stats_idx_lngth <- readRDS("output/hr_pol_sensitivity_idx_lngth.rds")

p_idx_lngth_ssb <- stats_idx_lngth %>%
  ggplot(aes(x = sigmaB, y = SSB_rel, colour = fhist)) +
  geom_line(show.legend = FALSE) +
  scale_colour_brewer(palette = "Dark2") +
  xlim(0, NA) + ylim(0, NA) +
  labs(x = expression(sigma[L]~"&"~sigma[B]), y = expression(SSB/B[MSY])) +
  theme_bw(base_size = 8)
p_idx_lngth_catch <- stats_idx_lngth %>%
  ggplot(aes(x = sigmaB, y = Catch_rel, colour = fhist)) +
  geom_line(show.legend = FALSE) +
  scale_colour_brewer(palette = "Dark2") +
  xlim(0, NA) + ylim(0, NA) +
  labs(x = expression(sigma[L]~"&"~sigma[B]), y = expression(Catch/MSY)) +
  theme_bw(base_size = 8)
p_idx_lngth_risk <- stats_idx_lngth %>%
  ggplot(aes(x = sigmaB, y = risk_Blim, colour = fhist)) +
  geom_line(show.legend = FALSE) +
  scale_colour_brewer(palette = "Dark2") +
  xlim(0, NA) + ylim(0, NA) +
  labs(x = expression(sigma[L]~"&"~sigma[B]), y = expression(B[lim]~risk)) +
  theme_bw(base_size = 8)

p_idx_lngth <- plot_grid(p_idx_lngth_ssb, p_idx_lngth_catch, p_idx_lngth_risk,
                     ncol = 3, align = "hv")
p_idx_lngth

### ------------------------------------------------------------------------ ###
### combine sensitivity plots ####
### ------------------------------------------------------------------------ ###
df_stats_idx_lngth <- stats_idx_lngth %>% 
  select(fhist, sigmaB, SSB_rel, Catch_rel, risk_Blim) %>%
  rename(x = sigmaB, SSB = SSB_rel, Catch = Catch_rel, risk = risk_Blim) %>%
  mutate(sensitivity = "obs") %>%
  pivot_longer(c(SSB, Catch, risk))
df_stats_rec <- stats_rec %>% 
  select(fhist, sigmaR, SSB_rel, Catch_rel, risk_Blim) %>%
  rename(x = sigmaR, SSB = SSB_rel, Catch = Catch_rel, risk = risk_Blim) %>%
  mutate(sensitivity = "rec") %>%
  pivot_longer(c(SSB, Catch, risk))
df_stats_status <- stats_status %>%
  filter(SSB0 <= 2) %>% ### remove SSB/Bmsy>2
  rename(x = SSB0) %>%
  mutate(sensitivity = "stock_status") %>%
  select(-iter) %>%
  pivot_longer(c(SSB, Catch, risk))
df_stats_time <- stats_time %>%
  filter(period == "cumulative") %>%
  rename(x = year) %>%
  mutate(sensitivity = "period") %>%
  pivot_longer(c(SSB, Catch, risk))

df_sens <- bind_rows(df_stats_idx_lngth, df_stats_rec, df_stats_status,
                     df_stats_time) %>%
  mutate(period = ifelse(is.na(period), "annual", period)) %>%
  mutate(label_top = factor(sensitivity, 
                            levels = c("rec", "obs", "stock_status", "period"), 
                            labels = c("recruitment~variability", 
                                       "observation~uncertainty",
                                       "initial~stock~status",
                                       "projection~period"))) %>%
  mutate(label_bottom = factor(sensitivity, 
                            levels = c("rec", "obs", "stock_status", "period"), 
                            labels = c("sigma[R]", 
                                       "sigma[obs]",
                                       "SSB[y==0]/B[MSY]",
                                       "years"))) %>%
  mutate(name = factor(name, levels = c("SSB", "Catch", "risk"),
                       labels = c("SSB/B[MSY]", "Catch/MSY", "B[lim]~risk"))) %>%
  mutate(value = ifelse(sensitivity == "stock_status" & fhist == "one-way",
                        NA, value))

df_blank <- data.frame(name = rep(c("SSB/B[MSY]", "Catch/MSY", "B[lim]~risk"),
                                  each = 2),
                       x = c(0, 1, 0, 1, 0, 1),
                       value = c(0, 1.8, 0, 1.4, 0, 1))

p_sens_rec <- df_sens %>%
  filter(sensitivity == "rec") %>%
  ggplot(aes(x = x, y = value, fill = fhist, colour = fhist, linetype = fhist)) +
  stat_smooth(n = 50, span = 0.2, se = FALSE, geom = "line", size = 0.4,
              show.legend = FALSE) + 
  geom_point(size = 0.3, stroke = 0, shape = 21, show.legend = FALSE) +
  geom_blank(data = df_blank) +
  facet_grid(name ~ label_top, scales = "free", labeller = "label_parsed",
             switch = "y") +
  scale_colour_brewer("", palette = "Dark2") +
  scale_fill_brewer("", palette = "Dark2") +
  labs(x = expression(sigma[R])) +
  theme_bw(base_size = 8) +
  theme(strip.placement = "outside",
        strip.background.y = element_blank(),
        axis.title.y = element_blank(), 
        strip.switch.pad.grid = unit(0, "pt"),
        plot.margin = unit(c(4, 2, 4, 4), "pt"))
p_sens_obs <- df_sens %>%
  filter(sensitivity == "obs") %>%
  ggplot(aes(x = x, y = value, fill = fhist, colour = fhist, linetype = fhist)) +
  stat_smooth(n = 50, span = 0.2, se = FALSE, geom = "line", size = 0.4,
              show.legend = FALSE) + 
  geom_point(size = 0.3, stroke = 0, shape = 21, show.legend = FALSE) +
  geom_blank(data = df_blank) +
  facet_grid(name ~ label_top, scales = "free", labeller = "label_parsed",
             switch = "y") +
  scale_colour_brewer("", palette = "Dark2") +
  scale_fill_brewer("", palette = "Dark2") +
  labs(x = expression(sigma[obs])) +
  theme_bw(base_size = 8) +
  theme(strip.placement = "outside",
        strip.text.y = element_blank(),
        strip.background.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        strip.switch.pad.grid = unit(0, "pt"),
        plot.margin = unit(c(4, 2, 4, 0), "pt"))
p_sens_status <- df_sens %>%
  filter(sensitivity == "stock_status") %>%
  ggplot(aes(x = x, y = value, fill = fhist, colour = fhist, linetype = fhist)) +
  stat_smooth(n = 50, span = 0.4, se = FALSE, geom = "line", size = 0.4,
              show.legend = FALSE) + 
  geom_point(size = 0.3, stroke = 0, shape = 21, show.legend = FALSE) +
  geom_blank(data = df_blank) +
  facet_grid(name ~ label_top, scales = "free", labeller = "label_parsed",
             switch = "y") +
  scale_colour_brewer("", palette = "Dark2") +
  scale_fill_brewer("", palette = "Dark2") +
  labs(x = expression(SSB[y == 0]/B[MSY])) +
  theme_bw(base_size = 8) +
  theme(strip.placement = "outside",
        strip.text.y = element_blank(),
        ### manual margins because no letter goes below base
        strip.text.x = element_text(margin = unit(c(3.8, 0, 3.8, 0), "pt")),
        strip.background.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        strip.switch.pad.grid = unit(0, "pt"),
        plot.margin = unit(c(4, 2, 4, 0), "pt"))
p_sens_period <- df_sens %>%
  filter(sensitivity == "period") %>%
  ggplot(aes(x = x, y = value, fill = fhist, colour = fhist, linetype = fhist)) +
  stat_smooth(n = 50, span = 0.1, se = FALSE, geom = "line", size = 0.4) + 
  geom_point(size = 0.3, stroke = 0, shape = 21) +
  geom_blank(data = df_blank) +
  facet_grid(name ~ label_top, scales = "free", labeller = "label_parsed",
             switch = "y") +
  scale_colour_brewer("fishing history", palette = "Dark2") +
  scale_fill_brewer("fishing history", palette = "Dark2") +
  scale_linetype("fishing history") +
  labs(x = "years") +
  theme_bw(base_size = 8) +
  theme(strip.placement = "outside",
        strip.text.y = element_blank(),
        strip.background.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        strip.switch.pad.grid = unit(0, "pt"),
        plot.margin = unit(c(4, 4, 4, 0), "pt"),
        legend.position = c(0.7, 0.45),
        legend.background = element_blank(),
        legend.key.height = unit(0.6, "lines"),
        legend.key = element_blank())

plot_grid(p_sens_rec, p_sens_obs, p_sens_status, p_sens_period, 
          nrow = 1, rel_widths = c(1.2, 1, 1, 1), align = "h")

ggsave(filename = "output/plots/hr_sensitivity_stats.png", type = "cairo",
       width = 17, height = 12, units = "cm", dpi = 600)
ggsave(filename = "output/plots/hr_sensitivity_stats.pdf",
       width = 17, height = 12, units = "cm")

### ------------------------------------------------------------------------ ###
### GA: pollack explorations ####
### ------------------------------------------------------------------------ ###

n_yrs <- 50
n_iter <- 500

### get optimised parameterisation
pol_GA <- foreach(stock = "pol", .combine = rbind) %:%
  foreach(optimised = c(TRUE, FALSE), .combine = rbind) %:%
  foreach(parameters = c("all", "caps", "upper_cap", "lower_cap", "multiplier",
                         "interval", "idx_range", "all_no_caps", "idx_lag",
                         "b_multiplier"), .combine = rbind) %:%
  foreach(objective = c("MSY", "MSY-PA"), .combine = rbind) %:%
  foreach(scenario = "GA", .combine = rbind) %:%
  foreach(catch_rule = "hr", .combine = rbind) %:%
  foreach(stat_yrs = "all", .combine = rbind) %:%
  foreach(fhist = c("one-way", "random"), .combine = rbind) %do% {#browser()
    
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
      "b_multiplier" = "comp_b_multiplier"
      )
    file_suffix <- switch(objective,
      "MSY" = paste0("obj_SSB_C_risk_ICV"),
      "MSY-PA" = "obj_ICES_MSYPA")
    file_res <- paste0(file_prefix, "--", file_suffix, "_res.rds")
    file_runs <- paste0(file_prefix, "--", file_suffix, "_runs.rds")
    
    if (!file.exists(paste0("output/", catch_rule, "/", n_iter, "_", n_yrs, "/", 
                            scenario, "/", fhist, "/", stock, "/", file_res))) { 
      return(NULL)
    }
    
    df_res <- data.frame(optimised = optimised, stock = stock, 
                         parameters = parameters, objective = objective, 
                         scenario = scenario, catch_rule = catch_rule, 
                         stat_yrs = stat_yrs, fhist = fhist)
    ### load GA results
    res <- readRDS(paste0("output/", catch_rule, "/", n_iter, "_", n_yrs, "/", 
                          scenario, "/", fhist, "/", stock, "/", file_res))
    runs <- readRDS(paste0("output/", catch_rule, "/", n_iter, "_", n_yrs, "/", 
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
    if (isFALSE(optimised)) {
      solution[] <- c(1, 1, 1, 1.4, 1, 1, Inf, 0)
    }
    df_res <- cbind(df_res, solution)
    solution_c <- paste0(solution, collapse = "_")
    ### replace NaN for upper_constraint with Inf
    solution_c <- gsub(x = solution_c, pattern = "NaN", replacement = "Inf")
    stats <- runs[[solution_c]]$stats
    stats <- as.data.frame(t(data.frame(unlist(stats), 
                                        row.names = rownames(stats))))
    row.names(stats) <- NULL
    df_res <- cbind(df_res, stats)
    ### prepare risk penalty
    penalty_tmp <- penalty(x = stats$risk_Blim, 
                           negative = FALSE, max = 5, 
                           inflection = 0.05 + 0.01, 
                           steepness = 0.5e+3)
    if (isTRUE(optimised)) {
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
pol_GA %>% View()

# pol_GA %>% 
#   filter(optimised == TRUE) %>%
#   bind_rows(pol_GA %>% 
#               filter(optimised == FALSE & 
#                        parameters == "all") %>%
#               mutate(parameters = "not_optimised")) %>%
#   mutate(parameters = factor(parameters,
#     levels = c("not_optimised", "idx_lag", "idx_range", "b_multiplier", 
#                "multiplier", "interval", "upper_cap", "lower_cap", "caps",
#                "all_no_caps", "all"))) %>%
#   ggplot(aes(x = as.factor(parameters), y = fitness, 
#              fill = as.factor(parameters))) +
#   geom_col(colour = "black") +
#   #scale_fill_brewer("parameters", palette = "Dark2") +
#   #brewer.pal(n = 8, name = "Dark2")
#   scale_fill_manual("parameters", 
#                     values = c("#D95F02", 
#                                colorRampPalette(c("grey30", "grey90"))(9), 
#                                "#1B9E77")) +
#   facet_wrap(~ fhist) +
#   labs(x = "optimised parameters", y = "fitness") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


### recreate fitness components
pol_GA_plot <- pol_GA %>%
  mutate(comp_Catch = Catch_rel - 1,
         comp_SSB = SSB_rel - 1,
         comp_ICV = ICV,
         comp_risk = risk_Blim,
         comp_risk_penalty = risk_penalty) %>%
  mutate(comp_risk = ifelse(objective == "MSY-PA", NA, comp_risk)) %>%
  pivot_longer(c(comp_Catch, comp_SSB, comp_ICV, comp_risk, comp_risk_penalty),
               names_prefix = "comp_") 
pol_GA_plot <- pol_GA_plot %>%
  filter(optimised == TRUE) %>%
  bind_rows(pol_GA_plot %>% 
              filter(optimised == FALSE & 
                       parameters == "all") %>%
              mutate(parameters = "not_optimised")) %>%
  mutate(parameters = factor(parameters,
    levels = c("not_optimised", "idx_lag", "idx_range", "b_multiplier", 
               "multiplier", "interval", "upper_cap", "lower_cap", "caps",
               "all_no_caps", "all"),
    labels = c("not optimised", "idx lag", "idx range", "b multiplier", 
               "multiplier", "interval", "upper cap", "lower cap", "caps",
               "all no caps", "all")
    )) %>%
  mutate(name = factor(name,
                       levels = rev(c("SSB", "Catch", "ICV", "risk",
                                      "risk_penalty")),
                       labels = rev(c("SSB", "Catch", "ICV", "risk",
                                   "risk penalty")))) %>%
  mutate(objective = factor(objective, levels = c("MSY", "MSY-PA"))) 

pol_GA_plot_dev <- pol_GA_plot %>%
  filter(name %in% c("SSB", "Catch")) %>%
  mutate(value_sign = ifelse(name == "SSB",
                             -abs(SSB_rel - 1)/2,
                             -abs(SSB_rel - 1) - abs(Catch_rel - 1)/2)) %>%
  mutate(sign = ifelse(name == "SSB",
                       ifelse(SSB_rel > 1, "+", "-"),
                       ifelse(Catch_rel > 1, "+", "-"))) %>%
  filter(abs(value) >= 0.02)


pol_GA_plot %>%
  mutate(value = -abs(value)) %>%
  ggplot(aes(x = parameters, y = value, fill = name)) +
  geom_col(colour = "black", size = 0.2) +
  # scale_fill_brewer("fitness\nelements", palette = "Dark2") +
  scale_fill_manual("fitness\nelements", 
                    values = brewer.pal(n = 5, name = "Dark2"),
                    breaks = rev(levels(pol_GA_plot$name))) +
  geom_text(data = pol_GA_plot_dev,
            aes(x = parameters, y = value_sign, label = sign), 
            vjust = 0.5, colour = "grey20") +
  facet_grid(objective ~ fhist) +
  coord_cartesian(ylim = c(-1.25, 0)) + 
  scale_y_continuous(breaks = c(0, -0.5, -1), 
                     minor_breaks = c(-0.25, -0.75, -1.25)) +
  labs(x = "optimised parameters", y = "fitness") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.spacing.x = unit(-0.5, "pt"),
        legend.key.width = unit(1, "lines"))
ggsave(filename = "output/plots/hr_GA_pol_comps.png", type = "cairo",
       width = 17, height = 12, units = "cm", dpi = 600)
ggsave(filename = "output/plots/hr_GA_pol_comps.pdf",
       width = 17, height = 12, units = "cm")


### ------------------------------------------------------------------------ ###
### GA: load results for all stocks ####
### ------------------------------------------------------------------------ ###
n_yrs <- 50
n_iter <- 500

### get optimised parameterisation
stocks_GA <- foreach(stock = stocks$stock, .combine = rbind) %:%
  foreach(optimised = c(TRUE, FALSE), .combine = rbind) %:%
  foreach(parameters = c("all", "multiplier"), .combine = rbind) %:%
  foreach(objective = c("MSY", "MSY-PA"), .combine = rbind) %:%
  foreach(scenario = "GA", .combine = rbind) %:%
  foreach(catch_rule = "hr", .combine = rbind) %:%
  foreach(stat_yrs = "all", .combine = rbind) %:%
  foreach(fhist = c("one-way", "random"), .combine = rbind) %do% {#browser()
    
    ### find files
    file_prefix <- switch(parameters,
      "all" = paste0("idxB_lag-idxB_range_3-comp_b_multiplier-interval-",
                     "multiplier-upper_constraint-lower_constraint"),
      "multiplier" = "multiplier",
    )
    file_suffix <- switch(objective,
                          "MSY" = paste0("obj_SSB_C_risk_ICV"),
                          "MSY-PA" = "obj_ICES_MSYPA")
    file_res <- paste0(file_prefix, "--", file_suffix, "_res.rds")
    file_runs <- paste0(file_prefix, "--", file_suffix, "_runs.rds")
    
    if (!file.exists(paste0("output/", catch_rule, "/", n_iter, "_", n_yrs, "/", 
                            scenario, "/", fhist, "/", stock, "/", file_res))) { 
      return(NULL)
    }
    
    df_res <- data.frame(optimised = optimised, stock = stock, 
                         parameters = parameters, objective = objective, 
                         scenario = scenario, catch_rule = catch_rule, 
                         stat_yrs = stat_yrs, fhist = fhist)
    ### load GA results
    res <- readRDS(paste0("output/", catch_rule, "/", n_iter, "_", n_yrs, "/", 
                          scenario, "/", fhist, "/", stock, "/", file_res))
    runs <- readRDS(paste0("output/", catch_rule, "/", n_iter, "_", n_yrs, "/", 
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
    if (isFALSE(optimised)) {
      solution[] <- c(1, 1, 1, 1.4, 1, 1, Inf, 0)
    }
    df_res <- cbind(df_res, solution)
    solution_c <- paste0(solution, collapse = "_")
    ### replace NaN for upper_constraint with Inf
    solution_c <- gsub(x = solution_c, pattern = "NaN", replacement = "Inf")
    stats <- runs[[solution_c]]$stats
    stats <- as.data.frame(t(data.frame(unlist(stats), 
                                        row.names = rownames(stats))))
    row.names(stats) <- NULL
    df_res <- cbind(df_res, stats)
    ### prepare risk penalty
    penalty_tmp <- penalty(x = stats$risk_Blim, 
                           negative = FALSE, max = 5, 
                           inflection = 0.05 + 0.01, 
                           steepness = 0.5e+3)
    if (isTRUE(optimised)) {
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
saveRDS(stocks_GA, "output/hr_GA_all_stocks_stats.rds")
stocks_GA <- readRDS("output/hr_GA_all_stocks_stats.rds")
stocks_GA %>% View()

### ------------------------------------------------------------------------ ###
### default vs. optimised: 4 example stocks ####
### ------------------------------------------------------------------------ ###

stocks_GA <- readRDS("output/hr_GA_all_stocks_stats.rds")

stocks_GA_plot <- stocks_GA %>%
  filter(stock %in% c("ang3", "pol", "tur", "san")) %>%
  mutate(parameters = as.character(parameters)) %>%
  mutate(parameters = ifelse(optimised == TRUE, parameters, "not optimised")) %>%
  unique() %>%
  mutate(stock = factor(stock, levels = c("ang3", "pol", "tur", "san"))) %>%
  mutate(objective = factor(objective, levels = c("MSY", "MSY-PA"))) %>%
  mutate(fhist = factor(fhist, levels = c("one-way", "random"))) %>%
  mutate(parameters = factor(parameters, 
                             levels = c("not optimised", "multiplier", "all"),
                             labels = c("not\noptimised", "multiplier", "all")))

stocks_GA_plot %>%
  ggplot(aes(x = stock, y = fitness, fill = parameters)) +
  geom_col(colour = "black", size = 0.2,
           position = "dodge") +
  #scale_fill_brewer("optimised\nparameters", palette = "Dark2") +
  scale_fill_manual("optimised\nparameters",
                    values = brewer.pal(n = 4, name = "Set1")[c(1, 2, 4)]
                    #breaks = rev(levels(pol_GA_plot$name))
                    ) +
  facet_grid(objective ~ fhist) +
  #coord_cartesian(ylim = c(-1.25, 0)) + 
  # scale_y_continuous(breaks = c(0, -0.5, -1), 
  #                    minor_breaks = c(-0.25, -0.75, -1.25)) +
  labs(x = "stocks", y = "fitness") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.spacing.x = unit(-0.5, "pt"),
        legend.key.width = unit(0.9, "lines"))
ggsave(filename = "output/plots/hr_GA_stocks_optimised.png", type = "cairo",
       width = 17, height = 10, units = "cm", dpi = 600)
ggsave(filename = "output/plots/hr_GA_stocks_optimised.pdf",
       width = 17, height = 10, units = "cm")

### ------------------------------------------------------------------------ ###
### comparison: HR vs. rfb vs. 2over3 ####
### ------------------------------------------------------------------------ ###

### load results from previous simulations
stats_rfb_PA <- readRDS("../wklife9_GA_tmp/output/all_stocks_PA_stats.rds")
stats_rfb_MSY <- readRDS("../wklife9_GA_tmp/output/all_stocks_MSY_stats.rds")
stats_2over3 <- readRDS("../wklife9_GA_tmp/output/all_stocks_MSY_stats.rds")
### load HR results
stats_hr <- readRDS("output/hr_GA_all_stocks_stats.rds")

### combine all results
stats_all_rules <- bind_rows(
  stats_rfb_PA %>%
    select(obj_fun, fhist, scenario, optimised, stock,
           lag_idx, range_idx_1, range_idx_2, range_catch, exp_r, exp_f, exp_b,
           interval, multiplier, upper_constraint, lower_constraint, 
           fitness, iter, risk_Blim, SSB_rel, Fbar_rel, Catch_rel, ICV) %>% View()
  
  stats_rfb_MSY %>%
    select(obj_fun, fhist, scenario, optimised, stock,
           lag_idx, range_idx_1, range_idx_2, range_catch, exp_r, exp_f, exp_b,
           interval, multiplier, upper_constraint, lower_constraint, 
           fitness, iter, risk_Blim, SSB_rel, Fbar_rel, Catch_rel, ICV) %>%
    mutate(rule = "rfb", obj_fun = "MSY")
  
  mutate()
  
)


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
  max_y_local <- max(stats3$catch, na.rm = TRUE)
  stats3 %>%
    ggplot(aes(x = Status, y = HR, fill = catch)) +
    facet_grid(years_label ~ stock_k, labeller = label_parsed) +
    geom_raster(interpolate = FALSE) +
    scale_fill_gradientn(expression(catch/MSY),
                         colours = c("red", "orange", "yellow", "green",
                                     "darkgreen"),
                         values = c(0, 0.33, 0.66, 1, max_y_global)/max_y_local
    ) +
    theme_bw(base_size = 7) +
    labs(x = expression(initial~stock~status~(SSB/B[MSY])),
         y = "harvest rate") +
    xlim(c(0, 3))
  ggsave(filename = paste0("output/plots/HR_principle_catch_examples",
                           names(groups)[i], ".png"), 
         type = "cairo", width = 17, height = 8, units = "cm", dpi = 600)
  ggsave(filename = paste0("output/plots/HR_principle_catch_examples", 
                           names(groups)[i], ".pdf"), 
         width = 17, height = 8, units = "cm")
}






