### ------------------------------------------------------------------------ ###
### create OMs ####
### ------------------------------------------------------------------------ ###

args <- commandArgs(TRUE)
print("arguments passed on to this script:")
print(args)

### evaluate arguments passed to R
for (i in seq_along(args)) eval(parse(text = args[[i]]))

### ------------------------------------------------------------------------ ###
### set up environment ####
### ------------------------------------------------------------------------ ###

### load packages
req_pckgs <- c("FLCore", "FLash", "FLBRP", "mse", "FLife", 
               "tidyr", "dplyr", "foreach", "doParallel")
for (i in req_pckgs) library(package = i, character.only = TRUE)

### load additional functions
source("funs.R")
source("funs_GA.R")

### ------------------------------------------------------------------------ ###
### setup parallel environment ####
### ------------------------------------------------------------------------ ###

if (isTRUE(n_workers > 1)) {
  ### start doParallel cluster
  cl <- makeCluster(n_workers)
  registerDoParallel(cl)
  cl_length <- length(cl)
  ### load packages and functions into parallel workers
  . <- foreach(i = seq(cl_length)) %dopar% {
    for (i in req_pckgs) library(package = i, character.only = TRUE,
                                 warn.conflicts = FALSE, verbose = FALSE,
                                 quietly = TRUE)
    source("funs.R", echo = FALSE)
    source("funs_GA.R", echo = FALSE)
  }
}

### ------------------------------------------------------------------------ ###
### fishing history dimensions ####
### ------------------------------------------------------------------------ ###

# n_iter <- 10000
# yrs_hist <- 100
# yrs_proj <- 100

set.seed(2)

### ------------------------------------------------------------------------ ###
### with uniform distribution and random F trajectories ####
### ------------------------------------------------------------------------ ###
# fhist <- "random"#"one-way"
if (identical(fhist, "random")) {
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
  
  f_array <- array(dim = c(yrs_hist, 3, n_iter),
                   dimnames = list(seq(yrs_hist), c("min","val","max"),
                                   iter = seq(n_iter)))
  f_array[, "val", ] <- c(t(df))
}

### ------------------------------------------------------------------------ ###
### create OMs ####
### ------------------------------------------------------------------------ ###

### get lhist for stocks
stocks <- read.csv("input/stocks.csv", stringsAsFactors = FALSE)

### BRPs from Fischer et al. (2020)
brps <- readRDS("input/brps.rds")

### create FLStocks
stocks_subset <- stocks$stock[stock_id]#"bll"

if (exists("OM")) { 
  
  if (isTRUE(OM)) {
    
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
      
      ### fishing history from previous paper
      if (isTRUE(fhist == "one-way")) {
        
        ### 0.5Fmsy until year 75, then increase to 0.8Fcrash
        fs <- rep(c(refpts["msy", "harvest"]) * 0.5, 74)
        f0 <- c(refpts["msy", "harvest"]) * 0.5
        fmax <- c(refpts["crash", "harvest"]) * 0.8
        rate <- exp((log(fmax) - log(f0)) / (25))
        fs <- c(fs, rate ^ (1:25) * f0)
        
        ### control object
        ctrl <- fwdControl(data.frame(year = 2:100, quantity = "f", val = fs))
        
      ### random F trajectories
      } else if (isTRUE(fhist == "random")) {
        
        ### control object template
        ctrl <- fwdControl(data.frame(year = seq(yrs_hist), 
                                      quantity = c("f"), val = NA))
        ### add iterations
        ctrl@trgtArray <- f_array
        ### target * Fcrash
        ctrl@trgtArray[,"val",] <- ctrl@trgtArray[,"val",] * 
          c(refpts["crash", "harvest"]) * 1
        
      }
      
      ### project fishing history
      stk_stf <- fwd(stk, ctrl, sr = stk_sr, sr.residuals = residuals(stk_sr),
                     sr.residuals.mult = TRUE, maxF = 5) 
      #plot(stk_stf, iter = 1:50)
      #plot(ssb(stk_stf), iter = 1:50)
      ### run a few times to get closer to target
      # for (i in 1:5) {
      #   stk_stf <- fwd(stk_stf, ctrl, sr = stk_sr,
      #                  sr.residuals.mult = TRUE, maxF = 4)
      # }
      name(stk_stf) <- stock
      path <- paste0("input/", n_iter, "_", yrs_proj, "/OM_1_hist/", fhist, "/")
      dir.create(path, recursive = TRUE)
      saveRDS(list(stk = stk_stf, sr = stk_sr), file = paste0(path, stock, ".rds"))
      
      return(NULL)
      #return(list(stk = stk_stf, sr = stk_sr))
    }
  }
}

### ------------------------------------------------------------------------ ###
### prepare OMs for flr/mse MP ####
### ------------------------------------------------------------------------ ###

if (exists("MP")) { 
  if (isTRUE(MP)) {
    
    stks_mp <- foreach(stock = stocks_subset, .errorhandling = "pass", 
                       .packages = c("FLCore", "mse")) %dopar% {
      ### load stock
      tmp <- readRDS(paste0("input/", n_iter, "_", yrs_proj, "/OM_1_hist/", fhist,
                            "/", stock, ".rds"))
      stk_fwd <- tmp$stk
      stk_sr <- tmp$sr
      rm(tmp); gc()
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
        idxB = tsb(stk_fwd),
        # idxB = quantSums(stk_fwd@stock.n * stk_fwd@stock.wt * (stk_fwd@mat %=% q)),
        idxL = lmean(stk = stk_fwd, params = pars_l),
        PA_status = ssb(stk_fwd) %=% NA_integer_)
      ### index deviation
      PA_status_dev <- FLQuant(NA, dimnames = list(age = c("positive", "negative"), 
                                                   year = dimnames(stk_fwd)$year, 
                                                   iter = dimnames(stk_fwd)$iter))
      set.seed(1)
      PA_status_dev["positive"] <- rbinom(n = PA_status_dev["positive"], 
                                          size = 1, prob = 0.9886215)
      set.seed(2)
      PA_status_dev["negative"] <- rbinom(n = PA_status_dev["negative"], 
                                          size = 1, prob = 1 - 0.4216946)
      set.seed(696)
      idx_dev <- FLQuants(sel = stk_fwd@mat %=% 1,
                          idxB = rlnoise(n = dims(idx$idxB)$iter, idx$idxB %=% 0, 
                                        sd = 0.2, b = 0),
                          idxL = rlnoise(n = dims(idx$idxL)$iter, idx$idxL %=% 0, 
                                         sd = 0.2, b = 0),
                          PA_status = PA_status_dev)
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
      ### harvest rates:
      set.seed(33)
      hr_rates <- runif(n = n_iter, min = 0, max = 1)
      ### parameters for components
      pars_est <- list(
        comp_r = FALSE, comp_f = FALSE, comp_b = FALSE,
        comp_c = FALSE, comp_m = hr_rates,
        idxB_lag = 1, idxB_range_1 = 2, idxB_range_2 = 3, idxB_range_3 = 1,
        catch_lag = 1, catch_range = 1,
        interval = 2,
        idxL_lag = 1, idxL_range = 1,
        exp_r = 1, exp_f = 1, exp_b = 1,
        Lref = rep((lhist$linf + 2*1.5*c(pars_l["Lc"])) / (1 + 2*1.5), n_iter),
        B_lim = rep(brps[[stock]]@Blim, n_iter),
        I_trigger = c(I_loss$idx_dev * 1.4), ### default, can be overwritten later
        pa_buffer = FALSE, pa_size = 0.8, pa_duration = 3,
        upper_constraint = Inf,
        lower_constraint = 0
      )
      
      ### operating model
      om <- FLom(stock = stk_fwd, ### stock 
                 sr = stk_sr, ### stock recruitment and precompiled residuals
                 fleetBehaviour = mseCtrl(),
                 projection = mseCtrl(method = fwd_attr,
                                      args = list(dupl_trgt = TRUE)))
      tracking = c("comp_c", "comp_i", "comp_r", "comp_f", "comp_b",
                   "multiplier", "exp_r", "exp_f", "exp_b")
      oem <- FLoem(method = obs_generic,
                   observations = list(stk = stk_fwd, idx = idx), 
                   deviances = list(stk = FLQuant(), idx = idx_dev),
                   args = list(idx_dev = TRUE, ssb = FALSE, tsb_idx = TRUE,
                               lngth = FALSE, lngth_dev = FALSE,
                               lngth_par = pars_l,
                               PA_status = FALSE, PA_status_dev = FALSE,
                               PA_Bmsy = c(refpts(brps[[stock]])["msy", "ssb"]), 
                               PA_Fmsy = c(refpts(brps[[stock]])["msy", "harvest"])))
      ctrl <- mpCtrl(list(
        est = mseCtrl(method = est_comps,
                           args = pars_est),
        phcr = mseCtrl(method = phcr_comps,
                            args = pars_est),
        hcr = mseCtrl(method = hcr_comps,
                           args = pars_est),
        isys = mseCtrl(method = is_comps,
                          args = pars_est)
      ))
      iem <- FLiem(method = iem_comps,
                   args = list(use_dev = FALSE, iem_dev = iem_dev))
      ### args
      args <- list(fy = dims(stk_fwd)$maxyear, ### final simulation year
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
      
      ### list with input to mp()
      input <- list(om = om, oem = oem, iem = iem, ctrl = ctrl, 
                    args = args,
                    scenario = "harvest_rates", tracking = tracking, 
                    verbose = TRUE,
                    refpts = refpts, Blim = Blim, I_loss = I_loss)
      
      ### save OM
      path <- paste0("input/", n_iter, "_", yrs_proj, "/OM_2_mp_input/", fhist, "/")
      dir.create(path, recursive = TRUE)
      saveRDS(object = input, file = paste0(path, stock, ".rds"))
      return(NULL)
    }
  }
}

# debugonce(wklife_3.2.1_est)
# debugonce(wklife_3.2.1_obs)
# debugonce(input$ctrl$hcr@method)
# debugonce(mp)
# debugonce(goFishDL)
# input$args$nblocks = 250
# res <- do.call(mp, input)
# 
# ### timing
# system.time({res1 <- do.call(mp, input)})


