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
                                    "SSB_idx_rfb_exp_error")))  {
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
                                    "SSB_idx_rfb_exp_error"))) {
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
                                    "SSB_idx_rfb_exp_error"))) {
    ### optimise exponents (weighting) of components:
    ### get stock to Bmsy
    obj <- -sum((SSBs - Bmsy)^2)
  }
  
  ### calculate some stats
  ret <- list()
  ret$fitness <- obj
  ret$Bmsy_dev <- median(c(sqrt((SSBs - Bmsy)^2), na.rm = TRUE))
  ret$Blim_risk <- mean(c(SSBs < Blim, na.rm = TRUE))
  ret$Bmsy_risk  <- mean(c(SSBs < Bmsy, na.rm = TRUE))
  ret$halfBmsy_risk <- mean(c(SSBs < 0.5*Blim, na.rm = TRUE))
  ret$collapse_risk <- mean(c(SSBs < 1, na.rm = TRUE))
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
