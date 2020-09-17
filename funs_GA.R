### ------------------------------------------------------------------------ ###
### objective function for multi species run ####
### ------------------------------------------------------------------------ ###
mp_fitness <- function(params, inp_file, path, check_file = FALSE,
                   scenario, 
                   return_res = FALSE,
                   collapse_correction = TRUE,
                   obj_SSB = TRUE, obj_C = TRUE, obj_F = FALSE,
                   obj_risk = TRUE, obj_ICV = TRUE, obj_ICES_PA = FALSE,
                   obj_ICES_PA2 = FALSE, obj_ICES_MSYPA = FALSE,
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
  params[10:11] <- round(params[10:11], 2)
  ### fix NaN for upper_constraint
  if (is.nan(params[10])) params[10] <- Inf
  
  ### check for files?
  if (isTRUE(check_file)) {
    ### current run
    run_i <- paste0(params, collapse = "_")
    ### check if path exists
    if (!dir.exists(path)) dir.create(path, recursive = TRUE)
    ### check if run already exists
    if (isTRUE(file.exists(paste0(path, run_i, ".rds")))) {
      ### load stats
      stats <- readRDS(paste0(path, run_i, ".rds"))
      ### set flag for running MP
      run_mp <- FALSE
    } else {
      run_mp <- TRUE
    }
  } else {
    run_mp <- TRUE
  }
  
  if (isTRUE(run_mp)) {
    
    ### load input file from disk
    input <- readRDS(inp_file)
    
    ### insert arguments into input object for mp
    input <- lapply(input, function(x) {
      x$ctrl.mp$ctrl.est@args$idxB_lag     <- params[1]
      x$ctrl.mp$ctrl.est@args$idxB_range_1 <- params[2]
      x$ctrl.mp$ctrl.est@args$idxB_range_2 <- params[3]
      x$ctrl.mp$ctrl.est@args$catch_range  <- params[4]
      x$ctrl.mp$ctrl.est@args$comp_m <- params[9]
      x$ctrl.mp$ctrl.phcr@args$exp_r <- params[5]
      x$ctrl.mp$ctrl.phcr@args$exp_f <- params[6]
      x$ctrl.mp$ctrl.phcr@args$exp_b <- params[7]
      x$ctrl.mp$ctrl.hcr@args$interval <- params[8]
      x$ctrl.mp$ctrl.hcr@args$interval <- params[8]
      x$ctrl.mp$ctrl.is@args$upper_constraint <- params[10]
      x$ctrl.mp$ctrl.is@args$lower_constraint <- params[11]
      
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
    
    ### calculate stats
    stats <- mp_stats(input = input, res_mp = res_mp, 
                      collapse_correction = collapse_correction)
    ### save result in file
    if (isTRUE(check_file)) {
      saveRDS(stats, paste0(path, run_i, ".rds"))
    }
    
  }
  
  ### objective function
  obj <- 0
  ### MSY objectives: target MSY reference values
  if (isTRUE(obj_SSB)) obj <- obj - sum(abs(unlist(stats["SSB_rel", ]) - 1))
  if (isTRUE(obj_C)) obj <- obj - sum(abs(unlist(stats["Catch_rel", ]) - 1))
  if (isTRUE(obj_F)) obj <- obj - sum(abs(unlist(stats["Fbar_rel", ]) - 1))
  ### reduce risk & ICV
  if (isTRUE(obj_risk)) obj <- obj - sum(unlist(stats["risk_Blim", ]))
  if (isTRUE(obj_ICV)) obj <- obj - sum(unlist(stats["ICV", ]))
  ### ICES approach: maximise catch while keeping risk <5%
  if (isTRUE(obj_ICES_PA)) {
    obj <- obj + sum(unlist(stats["Catch_rel", ]))
    ### penalise risk above 5%
    obj <- obj - sum(ifelse(test = unlist(stats["risk_Blim", ]) <= 0.05,
                            yes = 0,
                            no = 10)) 
  }
  if (isTRUE(obj_ICES_PA2)) {
    obj <- obj + sum(unlist(stats["Catch_rel", ]))
    ### penalise risk above 5% - gradual
    obj <- obj - sum(penalty(x = unlist(stats["risk_Blim", ]), 
                             negative = FALSE, max = 1, inflection = 0.06, 
                             steepness = 0.5e+3))
  }
  ### MSY target but replace risk with PA objective
    if (isTRUE(obj_ICES_MSYPA)) {
    obj <- obj - sum(abs(unlist(stats["SSB_rel", ]) - 1)) -
      sum(abs(unlist(stats["Catch_rel", ]) - 1)) -
      sum(unlist(stats["ICV", ])) -
      sum(penalty(x = unlist(stats["risk_Blim", ]), 
                             negative = FALSE, max = 5, inflection = 0.06, 
                             steepness = 0.5e+3))
      ### max penalty: 5
      ### for pollack zero catch has fitness of -4.7
  }
  
  ### housekeeping
  rm(res_mp, input)
  invisible(gc())
  if (getDoParWorkers() > 1)
    . <- foreach(i = 1:getDoParWorkers()) %dopar% {invisible(gc())}
  
  ### return objective function (fitness) value
  return(obj)
  
}

### ------------------------------------------------------------------------ ###
### stats from MSE run(s) ####
### ------------------------------------------------------------------------ ###

### function for calculating stats
mp_stats <- function(input, res_mp, collapse_correction = TRUE) {
  
  mapply(function(input_i, res_mp_i) {
    
    ### stock metrics
    SSBs <- FLCore::window(ssb(res_mp_i@stock), start = 101)
    Fs <- FLCore::window(fbar(res_mp_i@stock), start = 101)
    Cs <- FLCore::window(catch(res_mp_i@stock), start = 101)
    yrs <- dim(SSBs)[2]
    its <- dim(SSBs)[6]
    ### collapse correction
    if (isTRUE(collapse_correction)) {
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
    }
    
    Bmsy <- c(input_i$refpts["msy", "ssb"])
    Fmsy <- c(input_i$refpts["msy", "harvest"])
    Cmsy <- c(input_i$refpts["msy", "yield"])
    Blim <- input_i$Blim
    ### TAC interval
    TAC_intvl <- input_i$ctrl.mp$ctrl.hcr@args$interval
    
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
      ICV = iav(catch(res_mp_i@stock), from = 100, period = TAC_intvl,
                summary_all = median)
    )
    return(stats_i)
  }, input, res_mp)
  
}

### ------------------------------------------------------------------------ ###
### penalty function ####
### ------------------------------------------------------------------------ ###

penalty <- function(x, negative = FALSE, max = 1,
                    inflection = 0.06, steepness = 0.5e+3) {
  y <- max / (1 + exp(-(x - inflection)*steepness))
  if (isTRUE(negative)) y <- -y
  return(y)
}
