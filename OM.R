library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
#dir.create("R_lib")
.libPaths(c("R_lib", .libPaths()))
.libPaths()
#remotes::install_github("flr/FLash")
#remotes::install_github("flr/FLife")
#remotes::install_github("flr/FLBRP")
library(FLife)
library(FLash)
library(mseDL)
source("funs.R")
source("GA_funs.R")

### parallel environment
library(doParallel)
cl <- makeCluster(10)
registerDoParallel(cl)
clusterEvalQ(cl, {source("funs.R");source("GA_funs.R")})

### ------------------------------------------------------------------------ ###
### create F-based history ####
### ------------------------------------------------------------------------ ###

n_iter <- 25
yrs_hist <- 100
yrs_proj <- 50

set.seed(5)
f_dist <- rlnorm(n = n_iter, sdlog = 1)
hist(f_dist, breaks = 100)

f_array <- array(dim = c(yrs_hist, 3, n_iter),
                 dimnames = list(seq(yrs_hist), c("min","val","max"),
                                 iter = seq(n_iter)))
f_array[, "val", ] <- rep(f_dist, each = yrs_hist)
f_array[, "val", ] <- f_array[, "val", ] * rlnorm(n = length(f_array[, "val", ]),
                                                  meanlog = 0, sdlog = 0.25)
#names(dim(f_array)) <- 

### plot
# f_plot <- melt(f_array)
# names(f_plot) <- c("year", "target", "iter", "value")
# p1 <- ggplot(data = f_plot %>%
#          filter(target == "val") %>%
#          group_by(year) %>%
#          summarise(q5 = quantile(value, 0.05),
#                    q25 = quantile(value, 0.25),
#                    q50 = quantile(value, 0.5),
#                    q75 = quantile(value, 0.75),
#                    q95 = quantile(value, 0.95)),
#        aes(x = year - 100)) +
#   geom_ribbon(aes(ymin = q5, ymax = q95), colour = "grey", alpha = 0.5,
#               linetype = 0) +
#   geom_ribbon(aes(ymin = q25, ymax = q75), colour = "grey", alpha = 0.5,
#               linetype = 0) +
#   geom_line(aes(y = q50)) +
#   geom_line(data = f_plot %>%
#               filter(target == "val" & iter %in% 11:15),
#             aes(x = year - 100, y = value, colour = as.factor(iter)),
#             show.legend = FALSE, alpha = 0.8) +
#   theme_bw() +
#   labs(x = "year", y = expression(F/F[MSY])) +
#   coord_cartesian(ylim = c(0, 6.5), expand = FALSE) +
#   theme(panel.spacing = unit(0, units = "cm"))
# p2 <- ggplot(data = data.frame(x = f_dist), aes(x = x)) +
#   geom_density() +
#   coord_flip(expand = FALSE, xlim = c(0, 6.5), ylim = c(0, 0.55)) +
#   labs(x = "", y = "density") +
#   theme_bw() +
#   theme(axis.title.y = element_blank(),
#         axis.text.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         plot.margin = margin(c(1, 1, 1, 0)),
#         panel.spacing = unit(0, units = "cm"))
# plot_grid(p1, p2, rel_widths = c(1, 0.3), align = "h", nrow = 1)
# ggsave(filename = paste0("input/", n_iter, "_", yrs_proj, "/F_hist.png"),
#        width = 20, height = 14, units = "cm", dpi = 300, type = "cairo")

### ------------------------------------------------------------------------ ###
### create OMs ####
### ------------------------------------------------------------------------ ###
### use latest FLife

### get lhist for stocks
stocks <- read.csv("input/stocks.csv", stringsAsFactors = FALSE)

### create BRPs
brps <- lapply(stocks$stock, function(stock) {
  lh_i <- stocks[stocks$stock == stock, ]
  lh_i <- lh_i[, c("a", "b", "linf", "l50", "a50", "t0", "k")]
  ### default t0
  lh_i$t0 <- ifelse(is.na(lh_i$t0), -0.1, lh_i$t0)
  ### calc l50 if missing
  # lh_i$l50 <- ifelse(!is.na(lh_i$l50), lh_i$l50, 
  #                    vonB(age = lh_i$a50, 
  #                         params = FLPar(k = lh_i$k, linf = lh_i$linf,
  #                                        t0 = lh_i$t0)))
  # lh_i <- lh_i[, !is.na(lh_i)]
  ### estimate steepness h based on Wiff et al. 2018
  # lh_i$s <- h_Wiff(l50 = lh_i$l50, linf = lh_i$linf)
  lh_i$s <- 0.75
  lh_i <- as(lh_i, "FLPar")
  ### create missing values
  lh_i <- lhPar(lh_i)
  max_age <- ceiling(log(0.05)/(-c(lh_i["k"])) + c(lh_i["t0"]))
  ### create brp
  brp <- lhEql(lh_i, range = c(min = 1, max = max_age, 
                               minfbar = 1, maxfbar = max_age, 
                               plusgroup = max_age)
  )
  ### calculate blim, SSB where recruitment is 30% impaired
  bv <- function(SSB, a, b) a*SSB/(b + SSB)
  solve <- function(SSB) {
    rec = bv(a = c(params(brp)["a"]), 
             b = c(params(brp)["b"]), SSB = SSB)
    abs((c(refpts(brp)["virgin", "rec"]) * 0.7) - rec)
  }
  attr(brp, "Blim") <- optimize(f = solve, lower = 1, upper = 1000)$minimum
  
  return(brp)
})
names(brps) <- stocks$stock
saveRDS(brps, file = "input/OMs/brps.rds")
#brps <- readRDS("input/OMs/brps.rds")

### create FLStocks
stocks_subset <- stocks$stock#rev(stocks$stock[-c(24,29)])#stocks$stock[]
stks_hist <- foreach(stock = stocks_subset, .errorhandling = "pass", 
                     .packages = c("FLCore", "FLash", "FLBRP")) %dopar% {
  stk <- as(brps[[stock]], "FLStock")
  refpts <- refpts(brps[[stock]])
  stk <- qapply(stk, function(x) {#browser()
    dimnames(x)$year <- as.numeric(dimnames(x)$year) - 1; return(x)
  })
  stk <- stf(stk, yrs_hist + yrs_proj - dims(stk)$year + 1)
  stk <- propagate(stk, n_iter)
  ### create stock recruitment model
  stk_sr <- FLSR(params = params(brps[[stock]]), model = model(brps[[stock]]))
  ### create residuals for (historical) projection
  set.seed(0)
  residuals(stk_sr) <- rlnoise(dim(stk)[6], rec(stk) %=% 0, 
                               sd = 0.3, b = 0.2)
  ### fishing history
  ### fwcControl template
  ctrl <- fwdControl(data.frame(year = seq(yrs_hist), 
                                quantity = c("f"), val = NA))
  ### add iterations
  ctrl@trgtArray <- f_array
  ### target * Fmsy
  ctrl@trgtArray[,"val",] <- ctrl@trgtArray[,"val",] * 
    c(refpts["msy", "harvest"])
  ### project fishing history
  stk_stf <- fwd(stk, ctrl, sr = stk_sr,
                 sr.residuals.mult = TRUE, maxF = 5) 
  #plot(stk_stf, iter = 1:50)
  #plot(ssb(stk_stf), iter = 1:50)
  ### run a few times to get closer to target
  # for (i in 1:5) {
  #   stk_stf <- fwd(stk_stf, ctrl, sr = stk_sr,
  #                  sr.residuals.mult = TRUE, maxF = 4)
  # }
  name(stk_stf) <- stock
  path <- paste0("input/", n_iter, "_", yrs_proj, "/OM_1_hist/")
  dir.create(path, recursive = TRUE)
  saveRDS(list(stk = stk_stf, sr = stk_sr),
          file = paste0(path, stock, ".rds"))
  return(NULL)
}

# plot(stk_stf)
# plot(stk_stf, iter = 1:100)
# hist(ssb(stk_stf)[, ac(201)])
# ssb(stk_stf)[, ac(199:201),,,, 1]
# ctrl@trgtArray[ac(199:200),,]
# summary(c(ssb(stk_stf)[, ac(201)]) / ctrl@trgtArray[ac(200),"val",])
# plot(fbar(stk_stf)[, ac(101:200)])

### plot history for all stocks
for (stock in stocks_subset) {
  stk1 <- readRDS(paste0("input/500_50/OM_1_hist/", stock, ".rds"))$stk
  stk1 <- window(stk1, start = -100)
  stk1[, ac(-100:50)] <- stk1[, ac(0:150)]
  plot(window(stk1, end = 0),
       probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
       iter = 11:20) +
    ylim(0, NA) +
    labs(x = "year") +
    geom_hline(data = data.frame(qname = "SSB", 
                                 data = c(refpts(brps[[stock]])["msy", "ssb"])),
               aes(yintercept = data), linetype = "dashed", alpha = 0.5)
  ggsave(filename = paste0("input/", n_iter, "_", yrs_proj, "/SSB_hist_", 
                           stock, ".png"),
         width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
}
rm(stk1); gc()

### ------------------------------------------------------------------------ ###
### prepare OMs for flr/mse MP ####
### ------------------------------------------------------------------------ ###

stks_mp <- foreach(stock = stocks_subset, .errorhandling = "pass", 
                   .packages = c("FLCore", "mseDL")) %do% {
  ### load stock
  tmp <- readRDS(paste0("input/", n_iter, "_", yrs_proj, "/OM_1_hist/", stock, 
                        ".rds"))
  stk_fwd <- tmp$stk
  stk_sr <- tmp$sr
  ### life-history data
  lhist <- stocks[stocks$stock == stock, ]
  #range(stk_stf)
  ### cut of history
  stk_fwd <- window(stk_fwd, start = 75)
  stk_sr@residuals <- window(stk_sr@residuals, start = 75)
  ### length data
  pars_l <- FLPar(a = lhist$a,
                  b = lhist$b,
                  Lc = calc_lc(stk = stk_fwd[, ac(75:100)], 
                               a = lhist$a, b = lhist$b))
  ### indices
  q <- 1/(1 + exp(-1*(an(dimnames(stk_fwd)$age) - dims(stk_fwd)$max/10)))
  idx <- FLQuants(
    sel = stk_fwd@mat %=% q,
    idxB = ssb(stk_fwd),
    idxL = lmean(stk = stk_fwd, params = pars_l))
  ### index deviation
  set.seed(696)
  idx_dev <- FLQuants(sel = stk_fwd@mat %=% 1,
                      idxB = rlnoise(n = dims(idx$idxB)$iter, idx$idxB %=% 0, 
                                    sd = 0.2, b = 0),
                      idxL = rlnoise(n = dims(idx$idxL)$iter, idx$idxL %=% 0, 
                                     sd = 0.2, b = 0))
  ### iem deviation
  set.seed(205)
  iem_dev <- FLQuant(rlnoise(n = dims(stk_fwd)$iter,  catch(stk_fwd) %=% 0,
                            sd = 0.1, b = 0))
  ### lowest observed  index in last 25 years
  I_loss <- list()
  I_loss$SSB_idx <- apply(idx$idxB[, ac(75:100)], 6, min)
  I_loss$SSB_idx_dev <- apply((idx$idxB * idx_dev$idxB)[, ac(75:100)], 6, min)
  ### parameters for components
  pars_est <- list(
    comp_r = TRUE, comp_f = TRUE, comp_b = TRUE,
    idxB_lag = 1, idxB_range_1 = 2, idxB_range_2 = 3, idxB_range_3 = 1,
    catch_lag = 1, catch_range = 1,
    multiplier = 1,
    Lref = rep((lhist$linf + 2*1.5*c(pars_l["Lc"])) / (1 + 2*1.5), n_iter),
    idxL_lag = 1, idxL_range = 1,
    I_trigger = rep(brps[[stock]]@Blim, n_iter),
    exp_r = 1, exp_f = 1, exp_b = 1,
    interval = 2)
  
  ### operating model
  om <- FLom(stock = stk_fwd, ### stock 
             sr = stk_sr, ### stock recruitment and precompiled residuals
             fleetBehaviour = mseCtrl(),
             projection = mseCtrl(method = fwd_attr,
                                  args = list(dupl_trgt = TRUE)))
  tracking = c("C_current", "I_current", "comp_r", "comp_f", "comp_b",
               "multiplier", "exp_r", "exp_f", "exp_b")
  oem <- FLoem(method = wklife_3.2.1_obs,
               observations = list(stk = stk_fwd, idx = idx), 
               deviances = list(stk = FLQuant(), idx = idx_dev),
               args = list(idx_dev = FALSE, ssb = TRUE,
                           lngth = TRUE, lngth_dev = FALSE,
                           lngth_par = pars_l))
  ctrl.mp <- mpCtrl(list(
    ctrl.est = mseCtrl(method = wklife_3.2.1_est,
                       args = pars_est),
    ctrl.phcr = mseCtrl(method = phcr_r,
                        args = pars_est),
    ctrl.hcr = mseCtrl(method = hcr_r,
                       args = pars_est),
    ctrl.is = mseCtrl(method = is_r,
                      args = list(upper_constraint = Inf,
                                  lower_constraint = 0))
  ))
  # iem <- NULL
  iem <- FLiem(method = iem_r,
               args = list(use_dev = FALSE, iem_dev = iem_dev))
  ### genArgs
  genArgs <- list(fy = dims(stk_fwd)$maxyear, ### final simulation year
                  y0 = dims(stk_fwd)$minyear, ### first data year
                  iy = 100, ### first simulation (intermediate) year
                  nsqy = 3, ### not used, but has to provided
                  nblocks = 1, ### block for parallel processing
                  seed = 1, ### random number seed before starting MSE
                  seed_part = FALSE
  )
  ### get reference points
  refpts <- refpts(brps[[stock]])
  Blim <- attr(brps[[stock]], "Blim")
  
  ### list with input to mpDL()
  input <- list(om = om, oem = oem, iem = iem, ctrl.mp = ctrl.mp, 
                genArgs = genArgs,
                scenario = "SSB_idx_comp_r", tracking = tracking, 
                verbose = TRUE,
                refpts = refpts, Blim = Blim, I_loss = I_loss)
  
  ### save OM
  path <- paste0("input/", n_iter, "_", yrs_proj, "/OM_2_mp_input/")
  dir.create(path, recursive = TRUE)
  saveRDS(object = input, file = paste0(path, stock, ".rds"))
  return(NULL)
}

# debugonce(wklife_3.2.1_est)
# debugonce(wklife_3.2.1_obs)
# debugonce(input$ctrl.mp$ctrl.hcr@method)
# debugonce(mpDL)
# input$genArgs$nblocks = 2
# res <- do.call(mpDL, input)
# 
# ### timing
# system.time({res1 <- do.call(mpDL, input)})


