### ------------------------------------------------------------------------ ###
### objective function ####
### ------------------------------------------------------------------------ ###
mse_r <- function(params, input, path, check_file = FALSE,
                  scenario, 
                  return_res = FALSE,
                  ...) {
  
  ### housekeeping
  invisible(gc())
  if (exists("res_mp")) {
    rm(res_mp)
    invisible(gc())
  }
  
  ### rounding of arguments
  if (isTRUE(scenario %in% c("SSB_idx_r_only", "SSB_idx_r_only_error"))) {
    ### arguments for catch rule component r: use round values (years)
    params <- round(params)
  } else if (isTRUE(scenario %in% c("SSB_idx_rfb_exp", 
                                    "SSB_idx_rfb_exp_error",
                                    "SSB_idx_rfb_exp_error_Itrigger")))  {
    params[1:4] <- round(params[1:4])
    params[5:7] <- round(params[5:7], 1)
  }
  
  ### check for files?
  if (isTRUE(check_file)) {
    ### current run
    run_i <- paste0(params, collapse = "_")
    ### check if path exists
    if (!dir.exists(path)) dir.create(path, recursive = TRUE)
    ### check if run already exists
    if (isTRUE(file.exists(paste0(path, run_i, ".rds")))) {
      obj <- readRDS(paste0(path, run_i, ".rds"))$fitness
      ### if run already exists, return value and quit
      return(obj)
    }
  }
  
  ### insert arguments into input object for mp
  if (isTRUE(isTRUE(scenario %in% c("SSB_idx_r_only", "SSB_idx_r_only_error")))) {
    input$ctrl.mp$ctrl.est@args$idxB_lag     <- params[1]
    input$ctrl.mp$ctrl.est@args$idxB_range_1 <- params[2]
    input$ctrl.mp$ctrl.est@args$idxB_range_2 <- params[3]
    input$ctrl.mp$ctrl.est@args$catch_range  <- params[4]
  } else if (isTRUE(scenario %in% c("SSB_idx_rfb_exp", 
                                    "SSB_idx_rfb_exp_error",
                                    "SSB_idx_rfb_exp_error_Itrigger"))) {
    input$ctrl.mp$ctrl.est@args$idxB_lag     <- params[1]
    input$ctrl.mp$ctrl.est@args$idxB_range_1 <- params[2]
    input$ctrl.mp$ctrl.est@args$idxB_range_2 <- params[3]
    input$ctrl.mp$ctrl.est@args$catch_range  <- params[4]
    input$ctrl.mp$ctrl.phcr@args$exp_r <- params[5]
    input$ctrl.mp$ctrl.phcr@args$exp_f <- params[6]
    input$ctrl.mp$ctrl.phcr@args$exp_b <- params[7]
  }

  ### run MP
  res_mp <- do.call(mpDL, input)
  
  if (isTRUE(return_res)) {
    return(res_mp)
  }
  
  ### get some values
  SSBs <- FLCore::window(ssb(res_mp@stock), start = 101)
  Bmsy <- c(input$refpts["msy", "ssb"])
  Fmsy <- c(input$refpts["msy", "harvest"])
  Cmsy <- c(input$refpts["msy", "yield"])
  Blim <- input$Blim
  
  ### objective function
  if (isTRUE(scenario %in% c("SSB_idx_r_only", "SSB_idx_r_only_error"))) {
    ### optimise component r:
    ### keep SSB time series at start value
    ### SSB at begin of simulation
    SSBs_start <- SSBs
    SSBs_start[] <- ssb(res_mp@stock)[, "101"]
    ### compare risk during simulation with start
    risk_start <- mean(SSBs[, ac(101)] < input$Blim)
    risk_sim <- mean(SSBs < Blim)
    penalty <- ifelse(risk_sim > risk_start, risk_sim/risk_start, 1)
    ### objective: negative squared residuals
    ### GA maximises -> use negative
    obj <- -sum((SSBs - SSBs_start)^2) * penalty
  } else if (isTRUE(scenario %in% c("SSB_idx_rfb_exp",
                                    "SSB_idx_rfb_exp_error",
                                    "SSB_idx_rfb_exp_error_Itrigger"))) {
    ### optimise exponents (weighting) of components:
    ### get stock to Bmsy
    obj <- -sum((SSBs - Bmsy)^2)
  }
  
  ### calculate some stats
  ret <- list()
  ret$fitness <- obj
  ret$Bmsy_dev <- median(c(sqrt((SSBs - Bmsy)^2)), na.rm = TRUE)
  ret$Blim_risk <- mean(c(SSBs < Blim), na.rm = TRUE)
  ret$Bmsy_risk  <- mean(c(SSBs < Bmsy), na.rm = TRUE)
  ret$halfBmsy_risk <- mean(c(SSBs < 0.5*Blim), na.rm = TRUE)
  ret$collapse_risk <- mean(c(SSBs < 1), na.rm = TRUE)
  ret$catch_rel <- median(c(FLCore::window(catch(res_mp@stock), start = 101) / Cmsy))
  ret$SSB_rel <- median(c(FLCore::window(ssb(res_mp@stock), start = 101) / Bmsy))
  ret$F_rel <- median(c(FLCore::window(fbar(res_mp@stock), start = 101) / Fmsy))
  ret$icv <- iav(res_mp@stock@catch, period = 2, summary_all = median)


  ### save result in file
  if (isTRUE(check_file)) {
    saveRDS(ret, paste0(path, run_i, ".rds"))
  }
  
  ### housekeeping
  rm(res_mp)
  invisible(gc())
  
  ### return MSE object or fitness value?
  return(obj)
  
}


### ------------------------------------------------------------------------ ###
### objective function for multi species run ####
### ------------------------------------------------------------------------ ###
mse_ms <- function(params, inp_file, path, check_file = FALSE,
                   scenario, 
                   return_res = FALSE,
                   ...) {
  
  ### housekeeping
  invisible(gc())
  if (exists("res_mp")) {
    rm(res_mp)
    invisible(gc())
  }
  if (getDoParWorkers() > 1)
    . <- foreach(i = 1:getDoParWorkers()) %dopar% {invisible(gc())}
  
  ### rounding of arguments
  params[1:4] <- round(params[1:4])
  params[5:7] <- round(params[5:7], 1)
  params[8] <- round(params[8])
  params[9] <- round(params[9], 2)
  
  ### check for files?
  if (isTRUE(check_file)) {
    ### current run
    run_i <- paste0(params, collapse = "_")
    ### check if path exists
    if (!dir.exists(path)) dir.create(path, recursive = TRUE)
    ### check if run already exists
    if (isTRUE(file.exists(paste0(path, run_i, ".rds")))) {
      stats <- readRDS(paste0(path, run_i, ".rds"))
      ### if run already exists, return value and quit
      return(stats$obj)
    }
  }
  
  ### load input file from disk
  input <- readRDS(inp_file)
  
  ### insert arguments into input object for mp
  input <- lapply(input, function(x) {
    x$ctrl.mp$ctrl.est@args$idxB_lag     <- params[1]
    x$ctrl.mp$ctrl.est@args$idxB_range_1 <- params[2]
    x$ctrl.mp$ctrl.est@args$idxB_range_2 <- params[3]
    x$ctrl.mp$ctrl.est@args$catch_range  <- params[4]
    x$ctrl.mp$ctrl.phcr@args$exp_r <- params[5]
    x$ctrl.mp$ctrl.phcr@args$exp_f <- params[6]
    x$ctrl.mp$ctrl.phcr@args$exp_b <- params[7]
    x$ctrl.mp$ctrl.hcr@args$interval <- params[8]
    x$ctrl.mp$ctrl.phcr@args$multiplier <- params[9]
    return(x)
  })
  
  ### run MP for each list element
  res_mp <- lapply(input, function(x) {
    if (getDoParWorkers() > 1)
      . <- foreach(i = 1:getDoParWorkers()) %dopar% {invisible(gc())}
    do.call(mpDL, x)
  })
  
  if (isTRUE(return_res)) {
    return(res_mp)
  }
  
  ### some quants
  stats <- mapply(function(input_i, res_mp_i) {
    
    SSBs <- FLCore::window(ssb(res_mp_i@stock), start = 101)
    Fs <- FLCore::window(fbar(res_mp_i@stock), start = 101)
    Cs <- FLCore::window(catch(res_mp_i@stock), start = 101)
    Bmsy <- c(input_i$refpts["msy", "ssb"])
    Fmsy <- c(input_i$refpts["msy", "harvest"])
    Cmsy <- c(input_i$refpts["msy", "yield"])
    Blim <- input_i$Blim
    ### TAC interval
    TAC_intvl <- input_i$ctrl.mp$ctrl.hcr@args$interval
    ### SSB/F/Catch objectives
    SSB_obj <- median(abs(c(SSBs/Bmsy - 1)), na.rm = TRUE)
    F_obj   <- median(abs(c(Fs/Fmsy - 1)), na.rm = TRUE)
    C_obj   <- median(abs(c(Cs/Cmsy - 1)), na.rm = TRUE)
    ### risk, use Blim
    # risk_obj <- mean(c(SSBs < Bmsy/2), na.rm = TRUE)
    risk_obj <- mean(c(SSBs < Blim), na.rm = TRUE)
    ### ICV
    icv_obj <- iav(catch(res_mp_i@stock), from = 100, period = TAC_intvl,
                summary_all = median)
    ### objective function
    ### add components
    ### negative because GA maximises function
    # obj <- -(SSB_obj + F_obj + C_obj + risk_obj + icv_obj)
    obj <- -(SSB_obj + risk_obj + icv_obj)
    ### combine
    obj_i <- list(obj = obj,
                  obj_SSB = SSB_obj,
                  obj_F = F_obj,
                  obj_C = C_obj,
                  obj_risk = risk_obj,
                  obj_icv = icv_obj)
    
    ### some stats
    stats_i <- list(
      risk_Blim = mean(c(SSBs < Blim), na.rm = TRUE),
      risk_Bmsy = mean(c(SSBs < Bmsy), na.rm = TRUE),
      risk_halfBmsy = mean(c(SSBs < Bmsy/2), na.rm = TRUE),
      risk_collapse = mean(c(SSBs < 1), na.rm = TRUE),
      SSB = median(c(SSBs), na.rm = TRUE), Fbar = median(c(Fs), na.rm = TRUE),
      Catch = median(c(Cs), na.rm = TRUE),
      SSB_rel = median(c(SSBs/Bmsy), na.rm = TRUE),
      Fbar_rel = median(c(Fs/Fmsy), na.rm = TRUE),
      Catch_rel = median(c(Cs/Cmsy), na.rm = TRUE),
      ICV = icv_obj

    )

    ### combine
    res_i <- c(obj_i, stats_i)

    
  }, input, res_mp)
  
  stats <- list(obj = sum(unlist(stats["obj", ])),
                stats = stats)
  
  
  ### save result in file
  if (isTRUE(check_file)) {
    saveRDS(stats, paste0(path, run_i, ".rds"))
  }
  
  ### housekeeping
  rm(res_mp, input)
  invisible(gc())
  if (getDoParWorkers() > 1)
    . <- foreach(i = 1:getDoParWorkers()) %dopar% {invisible(gc())}
  
  ### return MSE object or fitness value?
  return(stats$obj)
  
}
