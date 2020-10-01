library(mse)
library(tidyverse)
library(doParallel)
library(scales)
library(cowplot)
source("funs.R")

cl <- makeCluster(3)
registerDoParallel(cl)
clusterEvalQ(cl = cl, expr = {library(mse)})

### stock list
stocks <- read.csv("input/stocks.csv", stringsAsFactors = FALSE)
stocks_subset <- stocks$stock[21:29]
names(stocks_subset) <- stocks_subset
### brps
brps <- readRDS("input/brps.rds")

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


# df_brill <- stats3 %>%
#   filter(stock == "bll")
max_y <- max(stats3$catch, na.rm = TRUE)
stats3 %>%
  ggplot(aes(x = status, y = HR, fill = catch)) +
  geom_raster() +
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
ggsave(filename = "output/plots/HR_catch.png", type = "cairo",
       width = 25, height = 15, units = "cm", dpi = 3840/(25/2.54))


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
       width = 10, height = 6, units = "cm", dpi = 3840/(10/2.54),
       type = "cairo")
### add LFeM
as.data.frame(idxL[,,,,, i]) %>%
  ggplot(aes(x = year - 100, y = data)) +
  geom_line() +
  geom_hline(yintercept = LFeM, colour = "black", linetype = "dashed") +
  theme_bw(base_size = 8) +
  labs(x = "year", y = "mean catch length [cm]")
ggsave(filename = "output/plots/length_procedure_lmean_2.png",
       width = 10, height = 6, units = "cm", dpi = 3840/(10/2.54),
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
       width = 10, height = 6, units = "cm", dpi = 3840/(10/2.54),
       type = "cairo")
### catch
as.data.frame(catch(stk)[,,,,, i]) %>%
  ggplot(aes(x = year - 100, y = data)) +
  geom_line() +
  theme_bw(base_size = 8) +
  labs(x = "year", y = "catch")
ggsave(filename = "output/plots/length_procedure_catch.png",
       width = 10, height = 6, units = "cm", dpi = 3840/(10/2.54),
       type = "cairo")
### biomass index
as.data.frame(idxB[,,,,, i]) %>%
  ggplot(aes(x = year - 100, y = data)) +
  geom_line() +
  theme_bw(base_size = 8) +
  labs(x = "year", y = "biomass index")
ggsave(filename = "output/plots/length_procedure_idxB.png",
       width = 10, height = 6, units = "cm", dpi = 3840/(10/2.54),
       type = "cairo")
### catch/biomass index
as.data.frame(idxB[,,,,, i]/catch(stk)[,,,,, i]) %>%
  ggplot(aes(x = year - 100, y = data)) +
  geom_line() +
  theme_bw(base_size = 8) +
  labs(x = "year", y = "catch/biomass index")
ggsave(filename = "output/plots/length_procedure_cr_1.png",
       width = 10, height = 6, units = "cm", dpi = 3840/(10/2.54),
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
       width = 10, height = 6, units = "cm", dpi = 3840/(10/2.54),
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
       width = 10, height = 6, units = "cm", dpi = 3840/(10/2.54),
       type = "cairo")


