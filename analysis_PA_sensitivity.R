### ------------------------------------------------------------------------ ###
### analyse sensitivity runs with PA fitness function ####
### ------------------------------------------------------------------------ ###

library(mse)
library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)
library(Cairo)


Blim <- readRDS("input/brps.rds")[["pol"]]@Blim
refpts <- refpts(readRDS("input/brps.rds")[["pol"]])

### ------------------------------------------------------------------------ ###
### collate data ####
### ------------------------------------------------------------------------ ###

runs2df <- function(x) {
  res <- lapply(x, function(y) {
    stats <- as.data.frame(t(y$stats))
    stats$stock <- colnames(y$stats)
    stats <- lapply(stats, unlist)
    pars <- y$pars
    df <- c(pars, stats)
    return(df)
  })
  res <- do.call(bind_rows, res)
  res <- as.data.frame(res)
  res$file <- names(x)
  return(res)
}

### get results
runs <- readRDS("output/500_50/risk/random/pol/all_runs.rds")
runs_df <- runs2df(runs)
### recruitment steepness
runs_steepness <- readRDS("output/500_50/PA_steepness/random/pol/all_runs.rds")
runs_steepness <- runs2df(runs_steepness)
### combine
runs_df <- bind_rows(runs_df, runs_steepness)

### save
saveRDS(runs_df, file = "output/pol_PA_sensitivity.rds")
write.csv(runs_df, file = "output/pol_PA_sensitivity.csv", row.names = FALSE)

### ------------------------------------------------------------------------ ###
### pollack: multiplier vs risk ####
### ------------------------------------------------------------------------ ###

df_mult <- readRDS("output/pol_PA_sensitivity.rds")
# df_mult <- runs_df
df_mult <- df_mult %>%
  filter(lag_idx == 1 & range_idx_1 == 2 & range_idx_2 == 3 & range_catch == 1 &
         exp_r == 1 & exp_f == 1 & exp_b == 1 & interval == 2 &
         upper_constraint == Inf & lower_constraint == 0 &
         sigmaL %in% c(0.2, 0.4, 0.6) & sigmaB %in% c(0.2, 0.4, 0.6) &
         (sigmaR == 0.6 | is.na(sigmaR)) & 
         (sigmaR_rho == 0 | is.na(sigmaR)) &
           is.na(steepness))

### plot Blim risk vs multiplier
# plot(df_mult$risk_Blim ~ df_mult$multiplier)

p_mult <- df_mult %>%
  filter(sigmaL == 0.2) %>%
  ggplot(aes(x = multiplier, y = risk_Blim)) +
  geom_hline(yintercept = 0.055, colour = "red", size = 0.4) +
  geom_smooth(method = "loess", span = 0.1, n = 10000, colour = "black", 
              level = 0, size = 0.5) +
  geom_point(size = 0.0, stroke = 0) +
  geom_vline(xintercept = 0.75, colour = "black", linetype = "dotted",
             size = 0.4) +
  geom_text(x = 0.78, y = 1, hjust = 0, size = 2, check_overlap = TRUE,
            label = expression(italic(x) == 0.75)) +
  theme_bw(base_size = 8) +
  labs(x = expression(multiplier~(italic(x))), y = expression(italic(B)[lim]~risk)) +
  ylim(c(0, 1))
p_mult
# ### for different levels of observation uncertainty -- not used
# df_mult %>%
#   filter(sigmaB %in% c(0.2, 0.4, 0.6)) %>%
#   ggplot(aes(x = multiplier, y = risk_Blim)) +
#   geom_smooth(method = "loess", span = 0.1, n = 10000, colour = "grey", 
#               level = 0, size = 0.5) +
#   geom_point(size = 0.5, stroke = 0) +
#   geom_hline(yintercept = 0.05, colour = "red", size = 0.4) +
#   geom_vline(xintercept = 0.705, colour = "black", linetype = "dotted",
#              size = 0.4) +
#   theme_bw() +
#   facet_wrap(~ paste0("CV = ", sigmaB), ncol = 1) +
#   labs(x = "catch rule multiplier", y = expression(italic(B)[lim]~risk)) +
#   scale_y_continuous(expand = expansion(mult = c(0, 0), add = 0),
#                      limits = c(0, 1)) +
#   scale_x_continuous(expand = expansion(mult = 0, add = 0))

### ------------------------------------------------------------------------ ###
### numerical approximation of risk slope -- not used ####
### ------------------------------------------------------------------------ ###

# NumDiffTable <- function(par, val, method = "central") {
#   ### forward difference approximation
#   if (identical(method, "forward")) {
#     diff <- sapply(seq_along(par)[1:(length(par) - 1)], function(x) {
#       (val[x + 1] - val[x])/(par[x + 1] - par[x])
#     })
#     diff <- c(diff, NA)
#   ### backward difference approximation
#   } else if (identical(method, "backward")) {
#     diff <- sapply(seq_along(par)[2:(length(par))], function(x) {
#       (val[x] - val[x - 1])/(par[x] - par[x - 1])
#     })
#     diff <- c(NA, diff)
#   ### use central difference approximation
#   } else if (identical(method, "central")) {
#     diff <- sapply(seq_along(par)[2:(length(par) - 1)], function(x) {
#       (val[x + 1] - val[x - 1])/(par[x + 1] - par[x - 1])
#     })
#     ### first/last values
#     diff0 <- (val[2] - val[1])/(par[2] - par[1])
#     diff1 <- (val[length(val)] - val[length(val) - 1]) /
#       (par[length(val)] - par[length(val) - 1])
#     diff <- c(diff0, diff, diff1)
#   }
#   return(diff)
# }
# 
# df_diff2 <- df_mult %>%
#   group_by(sigmaB) %>%
#   mutate(risk_Blim_diff = NumDiffTable(par = multiplier, val = risk_Blim,
#                                        method = "central")) %>%
#   mutate(risk_Blim_diff2 = NumDiffTable(par = multiplier, val = risk_Blim_diff,
#                                        method = "central"))
# df_plot2 <- df_diff2 %>%
#   select(multiplier, risk_Blim, risk_Blim_diff, sigmaB) %>%
#   pivot_longer(c(risk_Blim, risk_Blim_diff)) %>%
#   mutate(name = factor(name, levels = c("risk_Blim", "risk_Blim_diff"),
#                        labels = c("italic(B)[lim]~risk",
#                                   "B[lim]~risk~bold('\\'')")))
# 
# p <- df_plot2 %>%
#   ggplot(aes(x = multiplier, y = value)) +
#   geom_smooth(method = "loess", span = 0.2, n = 10000, colour = "grey",
#               level = 0, size = 0.5) +
#   geom_point(size = 0.5, stroke = 0) +
#   facet_grid(name ~ paste0("CV==", sigmaB), scales = "free", labeller = label_parsed) +
#   geom_hline(yintercept = 0.05, colour = "red", size = 0.4) +
#   geom_vline(xintercept = 0.705, 
#              colour = "black", linetype = "dotted", size = 0.4) +
#   theme_bw(base_size = 8) +
#   labs(x = "catch rule multiplier", y = "") #+
# p
# ggsave(filename = "output/plots/pol_risk_mult_obs-levels_diff.pdf",
#        width = 17, height = 7, units = "cm")
# p + coord_cartesian(ylim = c(0, 0.1), xlim = c(0, 1))
# ggsave(filename = "output/plots/pol_risk_mult_obs-levels_diff_zoom.pdf",
#        width = 17, height = 7, units = "cm")

### ------------------------------------------------------------------------ ###
### risk over time of optimised rfb-rule with 0.75 multiplier ####
### ------------------------------------------------------------------------ ###

input <- readRDS("input/sensitivity/500_100/OM_2_mp_input/random/pol.rds")
res <- readRDS(paste0("output/500_100/risk/random/pol/",
                      "mp_1_2_3_1_1_1_1_2_0.75_Inf_0_0.2_0.2.rds"))

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

### Blim risk
df_risk <- data.frame(year = 1:100)
Blim <- input$Blim
df_risk$risk <- sapply(1:100, function(x) {
  mean(c(SSBs[, ac(seq(from = 101, length.out = x))] < Blim), na.rm = TRUE)
})
df_risk$annual <- c(apply(SSBs < Blim, 2, mean, na.rm = TRUE))
df_risk <- df_risk %>%
  pivot_longer(c(risk, annual))


saveRDS(df_risk, file = "output/pol_PA_sensitivity_risk_100yrs.rds")
write.csv(df_risk, file = "output/pol_PA_sensitivity_risk_100yrs.csv", 
          row.names = FALSE)
df_risk <- readRDS("output/pol_PA_sensitivity_risk_100yrs.rds")

df_risk$name <- factor(df_risk$name, 
                       levels = c("annual", "risk"), 
                       labels = c("annual", "total"))

p_years <- df_risk %>%
  ggplot(aes(x = year, y = value, linetype = name)) +
  geom_hline(yintercept = 0.055, colour = "red", size = 0.4) +
  geom_line(size = 0.5) +
  theme_bw(base_size = 8) +
  scale_linetype("") +
  labs(x = "year", y = expression(italic(B)[lim]~risk)) +
  ylim(c(0, 0.5)) +
  geom_vline(xintercept = 50, colour = "black", linetype = "dotted",
             size = 0.4) +
  geom_text(x = 52, y = 0.5, hjust = 0, size = 2, check_overlap = TRUE,
            label = "default: 50 years") +
  theme(legend.position = c(0.74, 0.8),
        legend.background = element_blank(),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(1, "lines"),
        legend.key = element_blank())
p_years

### ------------------------------------------------------------------------ ###
### observation uncertainty with 0.75 multiplier ####
### ------------------------------------------------------------------------ ###

df_unc <- readRDS("output/pol_PA_sensitivity.rds")
# df_unc <- runs_df
df_unc <- df_unc %>%
  filter(lag_idx == 1 & range_idx_1 == 2 & range_idx_2 == 3 & range_catch == 1 &
         exp_r == 1 & exp_f == 1 & exp_b == 1 & interval == 2 &
         multiplier == 0.75 &
         upper_constraint == Inf & lower_constraint == 0 &
         (is.na(sigmaR) | sigmaR == 0.6) &
         (is.na(sigmaR_rho) | sigmaR_rho == 0.6) &
           is.na(steepness))

p_obs <- df_unc %>%
  ggplot(aes(x = sigmaB, y = risk_Blim)) +
  geom_hline(yintercept = 0.055, colour = "red", size = 0.4) +
  geom_smooth(method = "loess", span = 0.1, n = 10000, colour = "black",
              level = 0, size = 0.5) +
  geom_point(size = 0.0, stroke = 0) +
  geom_vline(xintercept = 0.2, colour = "black", linetype = "dotted",
             size = 0.4) +
  geom_text(x = 0.22, y = 0.5, hjust = 0, size = 2, check_overlap = TRUE,
            label = "default: CV=0.2") +
  ylim(c(0, 0.5)) + 
  theme_bw(base_size = 8) +
  labs(x = "observation uncertainty (CV)", 
       y = "")# +
p_obs

### ------------------------------------------------------------------------ ###
### recruitment variability with 0.75 multiplier ####
### ------------------------------------------------------------------------ ###

df_rec <- readRDS("output/pol_PA_sensitivity.rds")
# df_rec <- runs_df
df_rec <- df_rec %>%
  filter(lag_idx == 1 & range_idx_1 == 2 & range_idx_2 == 3 & range_catch == 1 &
         exp_r == 1 & exp_f == 1 & exp_b == 1 & interval == 2 &
           multiplier %in% c(0.75, 1) &
         upper_constraint == Inf & lower_constraint == 0 &
           sigmaL == 0.2 & sigmaB == 0.2 & sigmaR_rho == 0 &
           is.na(steepness))

p_rec <- df_rec %>%
  ggplot(aes(x = sigmaR, y = risk_Blim, linetype = as.factor(multiplier))) +
  geom_hline(yintercept = 0.055, colour = "red", size = 0.4) +
  geom_smooth(method = "loess", span = 0.5, n = 10000, colour = "black",
              level = 0, size = 0.5, se = FALSE) +
  geom_point(size = 0.0, stroke = 0) +
  scale_linetype("multiplier") +
  geom_vline(xintercept = 0.6, colour = "black", linetype = "dotted",
             size = 0.4) +
  geom_text(x = 0.62, y = 0.495, hjust = 0, size = 2, check_overlap = TRUE,
            #label = ("default: 0.6")
            label = expression("default:"~sigma[R]==0.6)) +
  ylim(c(0, 0.5)) + 
  theme_bw(base_size = 8) +
  labs(x = expression(recruitment~variability~"("*italic(sigma)[R]*")"), 
       y = "") +
  theme(legend.position = c(0.2, 0.8),
        legend.background = element_blank(),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(1, "lines"),
        legend.key = element_blank(),
        )
p_rec


### ------------------------------------------------------------------------ ###
### risk vs. starting status ####
### ------------------------------------------------------------------------ ###
### use 10,000 replicates to split into groups

res <- readRDS(paste0("output/10000_50/risk/random/pol/",
                      "mp_1_2_3_1_1_1_1_2_0.75_Inf_0_0.2_0.2.rds"))

### collapse correction
### stock metrics
SSBs <- FLCore::window(ssb(res$pol@stock), start = 100)
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

saveRDS(SSBs, file = "output/pol_PA_sensitivity_SSBs_10000.rds")
SSBs <- readRDS("output/pol_PA_sensitivity_SSBs_10000.rds")

### plot
p_Blim <- as.data.frame(window(SSBs, start = 101)) %>%
  ggplot(aes(x = data/c(refpts["msy", "ssb"]))) +
  geom_histogram(aes(y = (stat(count)/sum(count))*50),
                 binwidth = 0.05, colour = "grey", fill = "white",
                 show.legend = TRUE, size = 0.1) +
  geom_hline(yintercept = 0.055, colour = "red", size = 0.4) +
  stat_ecdf(aes(y = ..y..),
            pad = FALSE, geom = "step", size = 0.5) +
  # scale_y_continuous(labels = scales::percent, name = "frequency",
  # sec.axis = sec_axis(trans = ~ . / 0.02,
  #                     name = expression(cumulative~B[lim]~risk),
  #                     labels = scales::percent),
  #                     expand = expansion(mult = c(0, 0.05), add = 0)) +
  # scale_x_continuous(expand = c(0, 0), breaks = 0:6) +
  # coord_cartesian(xlim = c(NA, 7)) +
  xlim(c(0, 7)) +
  geom_vline(xintercept = Blim/c(refpts["msy", "ssb"]),
             colour = "black", linetype = "dotted",
             size = 0.4) +
  annotate(geom = "text",
           x = Blim/c(refpts["msy", "ssb"]) * 1.15, y = 0.99, 
            hjust = 0, size = 2, check_overlap = TRUE,
            label = expression(default:~italic(B)[lim]==0.57~italic(B)[MSY])) +
  theme_bw(base_size = 8) +
  labs(x = expression(italic(B)[lim]~"["*SSB/italic(B)[MSY]*"]"),
       y = "")
p_Blim

### ------------------------------------------------------------------------ ###
### risk vs. starting condition ####
### ------------------------------------------------------------------------ ###
### use results from previous section

### starting condition
SSBs0 <- SSBs[, ac(100)]
SSBs0 <- SSBs0/c(refpts["msy", "ssb"])
SSBs0 <- c(SSBs0)
SSB_breaks <- seq(from = 0, to = max(SSBs0), by = 0.1)
SSB_groups <- cut(SSBs0, breaks = SSB_breaks)

SSB_levels <- unique(as.character(SSB_groups))

### number of replicates per group:
risk_group_n <- sapply(SSB_levels, function(x) {
  length(which(SSB_groups %in% x))
})

### Blim risk per group
### SSB is on absolute scale 
risk_group <- sapply(SSB_levels, function(x) {
  tmp <- SSBs[, ac(101:150),,,, which(SSB_groups %in% x)]
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


p_initial <- df_risk %>%
  ggplot(aes(x = SSB0, y = risk)) +
  geom_hline(yintercept = 0.05, colour = "red", size = 0.4) +
  geom_smooth(method = "loess", span = 0.3, colour = "black",
              level = 0, size = 0.5) +
  geom_point(size = 0.0, stroke = 0) +
  theme_bw(base_size = 8) +
  labs(x = expression(SSB[italic(y) == 0]/italic(B)[MSY]), 
       y = "") +
  ylim(c(0, 0.2)) + #xlim(c(0, 4)) +
  coord_cartesian(xlim = c(0, 2))
p_initial

### ------------------------------------------------------------------------ ###
### recruitment steepness ####
### ------------------------------------------------------------------------ ###
df_h <- readRDS("output/pol_PA_sensitivity.rds")
# df_h <- runs_df
df_h <- df_h %>%
  filter(lag_idx == 1 & range_idx_1 == 2 & range_idx_2 == 3 & range_catch == 1 &
           exp_r == 1 & exp_f == 1 & exp_b == 1 & interval == 2 &
           multiplier == 0.75 &
           upper_constraint == Inf & lower_constraint == 0 &
           sigmaL == 0.2 & sigmaB == 0.2 & 
           sigmaR == 0.6 & sigmaR_rho == 0 &
           !is.na(steepness))

p_h <- df_h %>%
  ggplot(aes(x = steepness, y = risk_Blim)) +
  geom_hline(yintercept = 0.055, colour = "red", size = 0.4) +
  geom_smooth(method = "loess", span = 0.5, n = 10000, colour = "black",
              level = 0, size = 0.5, se = FALSE) +
  geom_point(size = 0.0, stroke = 0) +
  geom_vline(xintercept = 0.75, colour = "black", linetype = "dotted",
             size = 0.4) +
  geom_text(x = 0.37, y = 1.0, hjust = 0, size = 2, check_overlap = TRUE,
            #label = ("default: 0.6")
            label = expression("default:"~italic(h)==0.75)) +
  ylim(c(0, 1)) +
  theme_bw(base_size = 8) +
  labs(x = expression(recruitment~steepness~"("*italic(h)*")"), 
       y = expression(italic(B)[lim]~risk))
p_h

### ------------------------------------------------------------------------ ###
### combine sensitivity plots ####
### ------------------------------------------------------------------------ ###

p <- plot_grid(p_mult, p_Blim, p_initial, 
               p_years, p_obs, p_rec, 
               p_h,
               nrow = 3, ncol = 3, align = "hv", 
               labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)"), 
               hjust = 0, label_size = 10)
p 
ggsave(filename = "output/plots/PA/pol_PA_sensitivity.png", plot = p,
       width = 17, height = 12, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/PA/pol_PA_sensitivity.pdf", plot = p,
       width = 17, height = 12, units = "cm")
