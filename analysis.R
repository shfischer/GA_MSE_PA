library(mse)
library(ggplot2)
library(tidyr)
library(dplyr)


Blim <- readRDS("input/brps.rds")[["pol"]]@Blim
refpts <- refpts(readRDS("input/brps.rds")[["pol"]])

### ------------------------------------------------------------------------ ###
### pollack: multiplier vs risk ####
### ------------------------------------------------------------------------ ###

runs <- readRDS("output/500_50/risk/random/pol/all_runs.rds")

runs2df <- function(x) {
  res <- lapply(x, function(y) {
    stats <- as.data.frame(t(y$stats))
    stats$stock <- colnames(y$stats)
    stats <- lapply(stats, unlist)
    pars <- y$pars
    df <- c(pars, stats)
    return(df)
  })
  res <- do.call(rbind, res)
  res <- as.data.frame(res)
  res$file <- rownames(res)
  rownames(res) <- NULL
  res <- lapply(res, unlist)
  res <- as.data.frame(res)
  return(res)
}

df_mult <- runs2df(runs)
df_mult <- df_mult %>%
  filter(lag_idx == 1 & range_idx_1 == 2 & range_idx_2 == 3 & range_catch == 1 &
         exp_r == 1 & exp_f == 1 & exp_b == 1 & interval == 2 &
         upper_constraint == Inf & lower_constraint == 0 &
           sigmaL %in% c(0.2, 0.4, 0.6) & sigmaB %in% c(0.2, 0.4, 0.6))

### plot Blim risk vs multiplier
# plot(df_mult$risk_Blim ~ df_mult$multiplier)

df_mult %>%
  filter(sigmaL == 0.2) %>%
  ggplot(aes(x = multiplier, y = risk_Blim)) +
  geom_smooth(method = "loess", span = 0.1, n = 10000, colour = "grey", 
              level = 0, size = 0.5) +
  geom_point(size = 0.5, stroke = 0) +
  geom_hline(yintercept = 0.05, colour = "red", size = 0.4) +
  geom_vline(xintercept = 0.75, colour = "black", linetype = "dashed",
             size = 0.4) +
  theme_bw() +
  labs(x = "catch rule multiplier", y = expression(italic(B)[lim]~risk)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0), add = 0),
                     limits = c(0, 1)) +
  scale_x_continuous(expand = expansion(mult = 0, add = 0))
ggsave(filename = "output/plots/pol_risk_multiplier.pdf",
       width = 8.5, height = 6, units = "cm")
ggsave(filename = "output/plots/pol_risk_multiplier.png",
       width = 8.5, height = 6, units = "cm", dpi = 600, type = "cairo")

### for different levels of observation uncertainty
df_mult %>%
  filter(sigmaB %in% c(0.2, 0.4, 0.6)) %>%
  ggplot(aes(x = multiplier, y = risk_Blim)) +
  geom_smooth(method = "loess", span = 0.1, n = 10000, colour = "grey", 
              level = 0, size = 0.5) +
  geom_point(size = 0.5, stroke = 0) +
  geom_hline(yintercept = 0.05, colour = "red", size = 0.4) +
  geom_vline(xintercept = 0.705, colour = "black", linetype = "dashed",
             size = 0.4) +
  theme_bw() +
  facet_wrap(~ paste0("CV = ", sigmaB), ncol = 1) +
  labs(x = "catch rule multiplier", y = expression(italic(B)[lim]~risk)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0), add = 0),
                     limits = c(0, 1)) +
  scale_x_continuous(expand = expansion(mult = 0, add = 0))
ggsave(filename = "output/plots/pol_risk_multiplier_obs_levels.pdf",
       width = 8.5, height = 10, units = "cm")

### ------------------------------------------------------------------------ ###
### numerical approximation of risk slope ####
### ------------------------------------------------------------------------ ###

NumDiffTable <- function(par, val, method = "central") {
  ### forward difference approximation
  if (identical(method, "forward")) {
    diff <- sapply(seq_along(par)[1:(length(par) - 1)], function(x) {
      (val[x + 1] - val[x])/(par[x + 1] - par[x])
    })
    diff <- c(diff, NA)
  ### backward difference approximation
  } else if (identical(method, "backward")) {
    diff <- sapply(seq_along(par)[2:(length(par))], function(x) {
      (val[x] - val[x - 1])/(par[x] - par[x - 1])
    })
    diff <- c(NA, diff)
  ### use central difference approximation
  } else if (identical(method, "central")) {
    diff <- sapply(seq_along(par)[2:(length(par) - 1)], function(x) {
      (val[x + 1] - val[x - 1])/(par[x + 1] - par[x - 1])
    })
    ### first/last values
    diff0 <- (val[2] - val[1])/(par[2] - par[1])
    diff1 <- (val[length(val)] - val[length(val) - 1]) /
      (par[length(val)] - par[length(val) - 1])
    diff <- c(diff0, diff, diff1)
  }
  return(diff)
}

df_diff2 <- df_mult %>%
  group_by(sigmaB) %>%
  mutate(risk_Blim_diff = NumDiffTable(par = multiplier, val = risk_Blim,
                                       method = "central")) %>%
  mutate(risk_Blim_diff2 = NumDiffTable(par = multiplier, val = risk_Blim_diff,
                                       method = "central"))
df_plot2 <- df_diff2 %>%
  select(multiplier, risk_Blim, risk_Blim_diff, sigmaB) %>%
  pivot_longer(c(risk_Blim, risk_Blim_diff)) %>%
  mutate(name = factor(name, levels = c("risk_Blim", "risk_Blim_diff"),
                       labels = c("italic(B)[lim]~risk",
                                  "B[lim]~risk~bold('\\'')")))

p <- df_plot2 %>%
  ggplot(aes(x = multiplier, y = value)) +
  geom_smooth(method = "loess", span = 0.2, n = 10000, colour = "grey",
              level = 0, size = 0.5) +
  geom_point(size = 0.5, stroke = 0) +
  facet_grid(name ~ paste0("CV==", sigmaB), scales = "free", labeller = label_parsed) +
  geom_hline(yintercept = 0.05, colour = "red", size = 0.4) +
  geom_vline(xintercept = 0.705, 
             colour = "black", linetype = "dashed", size = 0.4) +
  theme_bw(base_size = 8) +
  labs(x = "catch rule multiplier", y = "") #+
p
ggsave(filename = "output/plots/pol_risk_mult_obs-levels_diff.pdf",
       width = 17, height = 7, units = "cm")
p + coord_cartesian(ylim = c(0, 0.1), xlim = c(0, 1))
ggsave(filename = "output/plots/pol_risk_mult_obs-levels_diff_zoom.pdf",
       width = 17, height = 7, units = "cm")



### ------------------------------------------------------------------------ ###
### default rfb-rule risk over time - 100 years ####
### ------------------------------------------------------------------------ ###

input <- readRDS("input/500_100/OM_2_mp_input/random/pol.rds")
res <- readRDS("output/500_100/risk/random/pol/mp_1_2_3_1_1_1_1_2_1_Inf_0_0.2_0.2.rds")

### collapse correction
### stock metrics
SSBs <- FLCore::window(ssb(res$pol@stock), start = 101)
Fs <- FLCore::window(fbar(res$pol@stock), start = 101)
Cs <- FLCore::window(catch(res$pol@stock), start = 101)
yrs <- dim(SSBs)[2]
its <- dim(SSBs)[6]
### collapse correction
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

### Blim risk
df_risk <- data.frame(year = 1:100)
df_risk$first <- sapply(1:100, function(x) {
  mean(c(SSBs[, ac(seq(from = 101, length.out = x))] < Blim), na.rm = TRUE)
})
df_risk$last <- sapply(1:100, function(x) {
  mean(c(SSBs[, ac(seq(from = 100 + x, to = 200))] < Blim), na.rm = TRUE)
})
df_risk <- df_risk %>%
  pivot_longer(c(first, last))

df_risk %>%
  ggplot(aes(x = year, y = value)) +
  geom_line() +
  facet_wrap(~ name) +
  theme_bw(base_size = 8) +
  labs(x = "year", y = expression(italic(B)[lim]~risk)) +
  ylim(c(0, NA))
ggsave(filename = "output/plots/pol_risk_time.pdf",
       width = 17, height = 6, units = "cm")

### ------------------------------------------------------------------------ ###
### risk over time of optimised rfb-rule with 0.75 multiplier ####
### ------------------------------------------------------------------------ ###

input <- readRDS("input/500_100/OM_2_mp_input/random/pol.rds")
res <- readRDS("output/500_100/risk/random/pol/mp_1_2_3_1_1_1_1_2_0.75_Inf_0_0.2_0.2.rds")

### collapse correction
### stock metrics
SSBs <- FLCore::window(ssb(res$pol@stock), start = 101)
Fs <- FLCore::window(fbar(res$pol@stock), start = 101)
Cs <- FLCore::window(catch(res$pol@stock), start = 101)
yrs <- dim(SSBs)[2]
its <- dim(SSBs)[6]
### collapse correction
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

### Blim risk
df_risk <- data.frame(year = 1:100)
Blim <- input$Blim
df_risk$risk <- sapply(1:100, function(x) {
  mean(c(SSBs[, ac(seq(from = 101, length.out = x))] < Blim), na.rm = TRUE)
})
df_risk %>%
  ggplot(aes(x = year, y = risk)) +
  geom_line() +
  theme_bw(base_size = 8) +
  labs(x = "year", y = expression(cumulative~B[lim]~risk)) +
  ylim(c(0, NA)) +
  geom_hline(yintercept = 0.05, colour = "red", size = 0.4)
ggsave(filename = "output/plots/pol_risk_time_0.75.pdf",
       width = 8.5, height = 6, units = "cm")
ggsave(filename = "output/plots/pol_risk_time_0.75.png",
       width = 8.5, height = 6, units = "cm", dpi = 600, type = "cairo")


### ------------------------------------------------------------------------ ###
### observation uncertainty with 0.75 multiplier ####
### ------------------------------------------------------------------------ ###


runs <- readRDS("output/500_50/risk/random/pol/all_runs.rds")


df_unc <- runs2df(runs)
df_unc <- df_unc %>%
  filter(lag_idx == 1 & range_idx_1 == 2 & range_idx_2 == 3 & range_catch == 1 &
         exp_r == 1 & exp_f == 1 & exp_b == 1 & interval == 2 &
           multiplier == 0.75 &
         upper_constraint == Inf & lower_constraint == 0)

df_unc %>%
  ggplot(aes(x = sigmaB, y = risk_Blim)) +
  geom_smooth(method = "loess", span = 0.1, n = 10000, colour = "grey",
              level = 0, size = 0.5) +
  geom_point(size = 0.5, stroke = 0) +
  geom_hline(yintercept = 0.05, colour = "red", size = 0.4) +
  geom_vline(xintercept = 0.2, colour = "black", linetype = "dashed",
             size = 0.4) +
  theme_bw() +
  labs(x = "observation uncertainty (CV)", 
       y = expression(B[lim]~risk)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0), add = c(0, 0)),
                     limits = c(0, 0.4)) +
  scale_x_continuous(expand = expansion(mult = 0, add = 0))
ggsave(filename = "output/plots/pol_risk_obs_unc.pdf",
       width = 8.5, height = 6, units = "cm")
ggsave(filename = "output/plots/pol_risk_obs_unc.png",
       width = 8.5, height = 6, units = "cm", dpi = 600, type = "cairo")


### ------------------------------------------------------------------------ ###
### risk vs. starting status ####
### ------------------------------------------------------------------------ ###
### use 10,000 replicates to split into groups

input <- readRDS("input/10000_50/OM_2_mp_input/random/pol.rds")
res_stats <- readRDS("output/10000_50/risk/random/pol/1_2_3_1_1_1_1_2_0.75_Inf_0_0.2_0.2.rds")
res <- readRDS("output/10000_50/risk/random/pol/mp_1_2_3_1_1_1_1_2_0.75_Inf_0_0.2_0.2.rds")

### collapse correction
### stock metrics
SSBs <- FLCore::window(ssb(res$pol@stock), start = 101)
yrs <- dim(SSBs)[2]
its <- dim(SSBs)[6]
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

### plot
as.data.frame(SSBs) %>%
  ggplot(aes(x = data/c(refpts["msy", "ssb"]))) +
  geom_histogram(aes(y = stat(count)/sum(count)),
                 binwidth = 0.05, colour = "grey", fill = "white",
                 show.legend = TRUE, size = 0.1) +
  stat_ecdf(aes(y = ..y.. * 0.02), 
            pad = FALSE, geom = "step", size = 0.5) +
  scale_y_continuous(labels = scales::percent, name = "frequency",
  sec.axis = sec_axis(trans = ~ . / 0.02,
                      name = expression(cumulative~B[lim]~risk),
                      labels = scales::percent),
                      expand = expansion(mult = c(0, 0.05), add = 0)) +
  scale_x_continuous(expand = c(0, 0), breaks = 0:6) +
  coord_cartesian(xlim = c(NA, 7)) +
  geom_hline(yintercept = 0.05*0.02, colour = "red", size = 0.4) +
  geom_vline(xintercept = Blim/c(refpts["msy", "ssb"]), 
             colour = "black", linetype = "dashed", 
             size = 0.4) + 
  theme_bw(base_size = 8) +
  labs(x = expression(SSB/B[MSY]))
ggsave(filename = "output/plots/pol_risk_SSB_Blim.pdf",
       width = 8.5, height = 6, units = "cm")
ggsave(filename = "output/plots/pol_risk_SSB_Blim.png",
       width = 8.5, height = 6, units = "cm", dpi = 600, type = "cairo")

### ------------------------------------------------------------------------ ###
### risk vs. starting condition ####
### ------------------------------------------------------------------------ ###
### use results from previous section

### starting condition
SSBs0 <- ssb(res$pol@stock)[, ac(100)]
SSBs0 <- SSBs0/c(refpts["msy", "ssb"])
SSBs0 <- c(SSBs0)
SSB_breaks <- seq(from = 0, to = max(SSBs0), by = 0.1)
SSB_groups <- cut(SSBs0, breaks = SSB_breaks)

SSB_levels <- unique(as.character(SSB_groups))

### Blim risk per group
### SSB is on absolut scale 
risk_group <- sapply(SSB_levels, function(x) {
  tmp <- SSBs[,,,,, which(SSB_groups %in% x)]
  mean(tmp < (Blim))
})
### get starting conditions
SSB_levels <- sapply(SSB_levels, function(x) {
  x <- gsub(x = x, pattern = "\\(|\\]", replacement = "")
  x <- unlist(strsplit(x, split = ","))
  mean(as.numeric(x))
})

df_risk <- data.frame(SSB0 = unlist(SSB_levels)[-length(SSB_levels)],
                      risk = unlist(risk_group)[-length(risk_group)])


df_risk %>%
  ggplot(aes(x = SSB0, y = risk)) +
  geom_smooth(method = "loess", span = 0.5, colour = "grey",
              level = 0, size = 0.5) +
  geom_point(size = 0.5, stroke = 0) +
  geom_hline(yintercept = 0.05, colour = "red", size = 0.4) +
  # geom_vline(xintercept = 0.75, colour = "black", linetype = "dashed",
  #            size = 0.4) +
  theme_bw() +
  labs(x = expression(initial~SSB/B[MSY]), 
       y = expression(italic(B)[lim]~risk)) +
  scale_y_continuous(expand = expansion(mult = c(0.00, 0.00), add = 0),
                     limits = c(0, 0.2)) +
  scale_x_continuous(expand = expansion(mult = 0, add = 0),
                     limits = c(0, 1.5))
ggsave(filename = "output/plots/pol_risk_start_level.pdf",
       width = 8.5, height = 6, units = "cm")
ggsave(filename = "output/plots/pol_risk_start_level.png",
       width = 8.5, height = 6, units = "cm", dpi = 600, type = "cairo")


