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
  if (isTRUE(all.equal(scenario, "SSB_idx_r_only"))) {
    ### arguments for catch rule component r: use round values (years)
    params <- round(params)
  } else if (isTRUE(all.equal(scenario, "SSB_idx_exp"))) {
    ### exponents for catch rule components: round to 1 decimal digit
    params <- round(params, 1)
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
    if (isTRUE(file.exists(paste0(path, run_i)))) {
      obj <- scan(file = paste0(path, run_i),
                  nlines = 1, nmax = 1, quiet = TRUE)
      ### if run already exists, return value and quit
      return(obj)
    }
  }
  
  ### insert arguments into input object for mp
  if (isTRUE(all.equal(scenario, "SSB_idx_r_only"))) {
    input$ctrl.mp$ctrl.est@args$idxB_lag     <- params[1]
    input$ctrl.mp$ctrl.est@args$idxB_range_1 <- params[2]
    input$ctrl.mp$ctrl.est@args$idxB_range_2 <- params[3]
    input$ctrl.mp$ctrl.est@args$catch_lag    <- params[4]
    input$ctrl.mp$ctrl.est@args$catch_range  <- params[5]
  } else if (isTRUE(all.equal(scenario, "SSB_idx_exp"))) {
    input$ctrl.mp$ctrl.phcr@args$exp_r <- params[1]
    input$ctrl.mp$ctrl.phcr@args$exp_f <- params[2]
    input$ctrl.mp$ctrl.phcr@args$exp_b <- params[3]
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
  
  ### objective function
  if (isTRUE(all.equal(scenario, "SSB_idx_r_only"))) {
    ### optimise component r:
    ### keep SSB time series at start value
    ssbs <- FLCore::window(ssb(res_mp@stock), start = 201)
    ### SSB at begin of simulation
    ssbs_start <- ssbs
    ssbs_start[] <- ssb(res_mp@stock)[, "201"]
    ### objective: negative squared residuals
    ### GA maximises -> use negative
    obj <- -sum((ssbs - ssbs_start)^2)
  } else if (isTRUE(scenario %in% c("SSB_idx_exp", "SSB_idx_rfb_exp",
                                    "SSB_idx_rfb_exp_error"))) {
    ### optimise exponents (weighting) of components:
    ### get stock to Bmsy
    Bmsy <- c(input$refpts["msy", "ssb"])
    ssbs <- FLCore::window(ssb(res_mp@stock), start = 201)
    obj <- -sum((ssbs - Bmsy)^2)
  }
  
  ### save result in file
  if (isTRUE(check_file)) {
    write(x = obj,
          file = paste0(path, run_i))
  }
  
  ### housekeeping
  rm(res_mp)
  invisible(gc())
  
  ### return MSE object or fitness value?
  return(obj)
  
}
