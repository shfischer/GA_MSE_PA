library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(FLife) ### GitHub SHA 25f481f1 2020-03-02
library(FLash)
library(mseDL)
source("funs.R")
source("GA_funs.R")

### parallel environment
library(doParallel)
cl <- makeCluster(10)
registerDoParallel(cl)
invisible(clusterEvalQ(cl, {source("funs.R"); source("GA_funs.R")}))

### ------------------------------------------------------------------------ ###
### create F-based history ####
### ------------------------------------------------------------------------ ###

n_iter <- 500
yrs_hist <- 100
yrs_proj <- 50

set.seed(2)


### ------------------------------------------------------------------------ ###
### create OMs ####
### ------------------------------------------------------------------------ ###
### use latest FLife

### get lhist for stocks
stocks <- read.csv("input/stocks.csv", stringsAsFactors = FALSE)

### use BRPs from paper
brps <- readRDS("input/OMs/brps.rds")

### create FLStocks
stocks_subset <- stocks$stock
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
                               sd = 0.6, b = 0)
  ### replicate residuals from catch rule paper for historical period
  set.seed(0)
  residuals <- rlnoise(dim(stk)[6], (rec(stk) %=% 0)[, ac(1:100)], 
                       sd = 0.6, b = 0)
  residuals(stk_sr)[, ac(1:100)] <- residuals[, ac(1:100)]
  
  ### fishing history

  ### 0.5Fmsy until year 75
  fs <- rep(c(refpts["msy", "harvest"]) * 0.5, 74)
  f0 <- c(refpts["msy", "harvest"]) * 0.5
  fmax <- c(refpts["crash", "harvest"]) * 0.8
  rate <- exp((log(fmax) - log(f0)) / (25))
  fs <- c(fs, rate ^ (1:25) * f0)
  
  ### control object
  ctrl <- fwdControl(data.frame(year = 2:100, quantity = "f", val = fs))
  ### project fishing history
  stk_stf <- fwd(stk, ctrl, sr = stk_sr, sr.residuals = residuals(stk_sr),
                 sr.residuals.mult = TRUE, maxF = 5) 
  name(stk_stf) <- stock
  path <- paste0("input/", n_iter, "_", yrs_proj, "/OM_1_hist/one-way/")
  dir.create(path, recursive = TRUE)
  saveRDS(list(stk = stk_stf, sr = stk_sr),
          file = paste0(path, stock, ".rds"))
  return(NULL)
  #return(list(stk = stk_stf, sr = stk_sr))
}
names(stks_hist) <- stocks_subset

### stock status
res <- lapply(stocks_subset, function(stock) {
  stk <- readRDS(paste0("input/500_50/OM_1_hist/one-way/", stock, ".rds"))$stk
  ssb(stk)[, ac(100)] / refpts(brps[[stock]])["msy", "ssb"]
})

### ------------------------------------------------------------------------ ###
### prepare OMs for flr/mse MP ####
### ------------------------------------------------------------------------ ###

stks_mp <- foreach(stock = stocks_subset, 
                   .packages = c("FLCore", "mseDL")) %do% {
  ### load stock
  tmp <- readRDS(paste0("input/", n_iter, "_", yrs_proj, "/OM_1_hist/one-way/", 
                        stock, ".rds"))
  stk_fwd <- tmp$stk
  stk_sr <- tmp$sr
  ### life-history data
  lhist <- stocks[stocks$stock == stock, ]
  #range(stk_stf)
  ### cut of history
  stk_fwd <- window(stk_fwd, start = 50)
  stk_sr@residuals <- window(stk_sr@residuals, start = 50)
  ### length data
  pars_l <- FLPar(a = lhist$a,
                  b = lhist$b,
                  Lc = calc_lc(stk = stk_fwd[, ac(75:100)], 
                               a = lhist$a, b = lhist$b))
  ### indices
  q <- 1/(1 + exp(-1*(an(dimnames(stk_fwd)$age) - dims(stk_fwd)$max/10)))
  idx <- FLQuants(
    sel = stk_fwd@mat %=% q,
    idxB = quantSums(stk_fwd@stock.n * stk_fwd@stock.wt * (stk_fwd@mat %=% q)),
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
  ### lowest observed index in last 50 years
  I_loss <- list()
  I_loss$SSB_idx <- apply(ssb(stk_fwd)[, ac(50:100)], 6, min)
  I_loss$SSB_idx_dev <- apply((ssb(stk_fwd) * idx_dev$idxB)[, ac(50:100)], 
                              6, min)
  I_loss$idx <- apply(idx$idxB[, ac(50:100)], 6, min)
  I_loss$idx_dev <- apply((idx$idxB * idx_dev$idxB)[, ac(50:100)], 6, min)
  ### parameters for components
  pars_est <- list(
    comp_r = TRUE, comp_f = TRUE, comp_b = TRUE,
    idxB_lag = 1, idxB_range_1 = 2, idxB_range_2 = 3, idxB_range_3 = 1,
    catch_lag = 1, catch_range = 1,
    multiplier = 1,
    Lref = rep((lhist$linf + 2*1.5*c(pars_l["Lc"])) / (1 + 2*1.5), n_iter),
    idxL_lag = 1, idxL_range = 1,
    exp_r = 1, exp_f = 1, exp_b = 1,
    interval = 2,
    B_lim = rep(brps[[stock]]@Blim, n_iter),
    I_trigger = c(I_loss$idx_dev * 1.4) ### default, can be overwritten later
  )
  
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
               args = list(idx_dev = TRUE, ssb = FALSE,
                           lngth = TRUE, lngth_dev = TRUE,
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
  iem <- FLiem(method = iem_r,
               args = list(use_dev = TRUE, iem_dev = iem_dev))
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
  path <- paste0("input/", n_iter, "_", yrs_proj, "/OM_2_mp_input/one-way/")
  dir.create(path, recursive = TRUE)
  saveRDS(object = input, file = paste0(path, stock, ".rds"))
  return(NULL)
}

# debugonce(wklife_3.2.1_est)
# debugonce(wklife_3.2.1_obs)
# debugonce(input$ctrl.mp$ctrl.hcr@method)
# debugonce(mpDL)
# debugonce(goFishDL)
# input$genArgs$nblocks = 250
# res <- do.call(mpDL, input)
# 
# ### timing
# system.time({res1 <- do.call(mpDL, input)})


