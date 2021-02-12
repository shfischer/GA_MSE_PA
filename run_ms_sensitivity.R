### ------------------------------------------------------------------------ ###
### run MSE without optimisation ####
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### arguments ####
### ------------------------------------------------------------------------ ###

args <- commandArgs(TRUE)
print("arguments passed on to this script:")
print(args)

### evaluate arguments, if they are passed to R:
if (length(args) > 0) {
  
  ### extract arguments
  for (i in seq_along(args)) eval(parse(text = args[[i]]))
  ### set default arguments
  ### parallelization
  if (!exists("use_MPI")) use_MPI <- FALSE
  if (!exists("n_blocks")) n_blocks <- 1
  if (!exists("n_workers")) n_workers <- 0
  ### split OM into blocks?
  if (!exists("n_parts")) n_parts <- 1
  if (!exists("part")) part <- 1
  ### projection details
  if (!exists("n_iter")) n_iter <- 500
  if (!exists("n_yrs")) n_yrs <- 100
  if (!exists("fhist")) fhist <- "random"
  if (!exists("catch_rule")) catch_rule <- "catch_rule"
  
  if (!exists("scenario")) scenario <- "risk"
  ### uncertainty
  if (!exists("sigmaL")) sigmaL <- 0.2
  if (!exists("sigmaB")) sigmaB <- 0.2
  ### what to save
  if (!exists("saveMP")) saveMP <- FALSE
  if (!exists("cut_hist")) cut_hist <- TRUE
  if (!exists("collate")) collate <- FALSE

} else {
  
  stop("no argument passed to R")
  
}

### ------------------------------------------------------------------------ ###
### set up environment ####
### ------------------------------------------------------------------------ ###

### load packages
### GA fork from GitHub remotes::install_github("shfischer/GA")
### use mse fork from shfischer/mse, branch mseDL2.0 
### remotes::install_github("shfischer/mse", ref = "mseDL2.0")
req_pckgs <- c("FLCore", "FLash", "mse", "GA", "doParallel", "doRNG", "FLBRP")
for (i in req_pckgs) library(package = i, character.only = TRUE)

### load additional functions
source("funs.R")
source("funs_GA.R")

### ------------------------------------------------------------------------ ###
### setup parallel environment ####
### ------------------------------------------------------------------------ ###

### hybrid MPI
if (isTRUE(use_MPI)) {
  ### 1st: doMPI cluster with 1 worker per node
  message("starting doMPI")
  library(doMPI)
  cl1 <- startMPIcluster()
  message("startMPIcluster() succeeded")
  print(cl1)
  registerDoMPI(cl1)
  cl_length_1 <- cl1$workerCount
  cl_length_1
  
  ### 2nd: doParallel workers inside doMPI workers
  . <- foreach(i = seq(cl_length_1)) %dopar% {
    ### load packages and functions into MPI workers
    for (i in req_pckgs) library(package = i, character.only = TRUE,
                                 warn.conflicts = FALSE, verbose = FALSE,
                                 quietly = TRUE)
  }
  message("MPI package loading succeeded")
  . <- foreach(i = seq(cl_length_1)) %dopar% {
    source("funs.R", echo = FALSE)
    source("funs_GA.R", echo = FALSE)
  }
  message("MPI script loading succeeded")
  ### start doParallel inside MPI processes
  if (isTRUE(n_workers > 1)) {
    . <- foreach(i = seq(cl_length_1)) %dopar% {
      cl2 <- makeCluster(n_workers)
      registerDoParallel(cl2)
      cl_length_2 <- length(cl2)
      ### load packages and functions into parallel workers
      . <- foreach(i = seq(cl_length_2)) %dopar% {
        for (i in req_pckgs) library(package = i, character.only = TRUE,
                                     warn.conflicts = FALSE, verbose = FALSE,
                                     quietly = TRUE)
        source("funs.R", echo = FALSE)
        source("funs_GA.R", echo = FALSE)
      }
    }
  }
  message("setting up doParallel inside MPI succeeded")
} else {
  if (isTRUE(n_workers > 1)) {
    ### start doParallel cluster
    cl1 <- makeCluster(n_workers)
    registerDoParallel(cl1)
    cl_length_1 <- length(cl1)
    ### load packages and functions into parallel workers
    . <- foreach(i = seq(cl_length_1)) %dopar% {
      for (i in req_pckgs) library(package = i, character.only = TRUE,
                                   warn.conflicts = FALSE, verbose = FALSE,
                                   quietly = TRUE)
      source("funs.R", echo = FALSE)
      source("funs_GA.R", echo = FALSE)
    }
  } else {
    cl1 <- FALSE
  }
}

### ------------------------------------------------------------------------ ###
### load data ####
### ------------------------------------------------------------------------ ###

stocks <- read.csv("input/stocks.csv", stringsAsFactors = FALSE)
stock <- stocks$stock[stock_id]
names(stock) <- stock
input <- lapply(stock, function(x) {
  readRDS(paste0("input/", n_iter, "_", n_yrs, "/OM_2_mp_input/", fhist, "/", x,
                        ".rds"))
})

### ------------------------------------------------------------------------ ###
### split simulation into junks ####
### ------------------------------------------------------------------------ ###

if (isTRUE(n_parts > 1)) {
  
  its <- split(seq(n_iter), sort(seq(n_iter) %% n_parts) + 1)[[part]]
  input <- lapply(input, function(x) {#browser()
    
    x_i <- x
    ### OM
    x_i$om@stock <- iter_attr(x_i$om@stock, its)
    x_i$om@sr <- FLCore::iter(x_i$om@sr, its)
    x_i$oem <- iters(x_i$oem, its)
    x_i$iem <- iters(x_i$iem, its)
    x_i$ctrl.mp <- iters(x_i$ctrl.mp, its)
    ### subset ctrl.mp arguments
    # for (z in seq_along(length(x_i$ctrl.mp))) {
    #   x_i$ctrl.mp[[z]]@args <- lapply(x_i$ctrl.mp[[z]]@args, function(y) {
    #     if (isTRUE(length(y) == n_iter)) {
    #       return(y[its])
    #     } else {
    #       return(y)
    #     }
    #   })
    # }
    ### no, these are subset inside the catch rule...
    # x_i$I_loss <- lapply(x_i$I_loss, function(y) {
    #   FLCore::iter(y, its)
    # }) ### not used
    return(x_i)
    
  })
  gc()
}

### ------------------------------------------------------------------------ ###
### scenario set-up ####
### ------------------------------------------------------------------------ ###

### catch rule parameters
rfb_names <- c("lag_idx", "range_idx_1", "range_idx_2", "range_catch",
             "exp_r", "exp_f", "exp_b", "interval", "multiplier",
             "upper_constraint", "lower_constraint")
rfb_pars <- c(1, 2, 3, 1, 1, 1, 1, 2, 1, Inf, 0)
### uncertainty parameters
unc_names <- c("sigmaL", "sigmaB")
unc_pars <- c(0.2, 0.2)
### combine
pars_names <- c(rfb_names, unc_names)
pars_def <- as.list(c(rfb_pars, unc_pars))
names(pars_def) <- pars_names
### load parameter values, if provided
pars <- mget(pars_names, ifnotfound = FALSE)
### load default values
pars_overwrite <- which(sapply(pars, isFALSE))
pars[pars_overwrite] <- pars_def[pars_overwrite]
### coerce into table
pars <- as.data.frame(do.call(cbind, pars))

### ------------------------------------------------------------------------ ###
### define dimension for parallelisation ####
### ------------------------------------------------------------------------ ###

### parallelise over parameter values?
`%do_pars%` <- ifelse(isTRUE(nrow(pars) > 1),
                        foreach::`%dopar%`,
                        foreach::`%do%`)

### parallelise over pars?
`%do_parts%` <- ifelse(isTRUE(n_parts > 1), 
                        foreach::`%dopar%`, 
                        foreach::`%do%`)

### output path
path_out <- paste0("output/", n_iter, "_", n_yrs, "/", scenario, "/",
                   fhist, "/",
                   paste0(stock, collapse = "_"), "/")
dir.create(path_out, recursive = TRUE)


. <- foreach(par_i = split(pars, seq(nrow(pars)))) %do_pars% {
  
  input_i <- input
  
  ### within run parallelisation?
  if (isTRUE(n_workers > 1) & isTRUE(n_blocks > 1)) {
    input_i <- lapply(input_i, function(x) {
      x$args$nblocks <- n_blocks
      return(x)
    })
  }
  
  ## ---------------------------------------------------------------------- ###
  ## catch rule parameters ####
  ## ---------------------------------------------------------------------- ###
  input_i <- lapply(input_i, function(x) {
    x$ctrl$est@args$idxB_lag     <- par_i$lag_idx
    x$ctrl$est@args$idxB_range_1 <- par_i$range_idx_1
    x$ctrl$est@args$idxB_range_2 <- par_i$range_idx_2
    x$ctrl$est@args$catch_range  <- par_i$range_catch
    x$ctrl$est@args$comp_m <- par_i$multiplier
    x$ctrl$phcr@args$exp_r <- par_i$exp_r
    x$ctrl$phcr@args$exp_f <- par_i$exp_f
    x$ctrl$phcr@args$exp_b <- par_i$exp_b
    x$ctrl$hcr@args$interval <- par_i$interval
    x$ctrl$isys@args$upper_constraint <- par_i$upper_constraint
    x$ctrl$isys@args$lower_constraint <- par_i$lower_constraint
    return(x)
  })
  
  ## ---------------------------------------------------------------------- ###
  ## uncertainty ####
  ## ---------------------------------------------------------------------- ###
  ### change uncertainty?
  sigmaB_i <- par_i$sigmaB
  sigmaL_i <- par_i$sigmaL
  if (par_i$sigmaB != pars_def$sigmaB | par_i$sigmaL != pars_def$sigmaL) {
  
    input_i <- lapply(input_i, function(x) {
      
      #browser()
      ### create observation noise
      set.seed(695)
      dev_idxB <- x$oem@deviances$idx$idxB
      dev_idxL <- x$oem@deviances$idx$idxL
      dev_idxB[] <- rlnoise(n = dims(dev_idxB)$iter, dev_idxB %=% 0, 
                                      sd = sigmaB_i, b = 0)
      dev_idxL[] <- rlnoise(n = dims(dev_idxL)$iter, dev_idxL %=% 0, 
                                  sd = sigmaL_i, b = 0)
      set.seed(696)
      dev_idxB[, ac(50:150)] <- rlnoise(n = dims(dev_idxB)$iter,
                                        window(dev_idxB, end = 150) %=% 0,
                                        sd = sigmaB_i, b = 0)
      dev_idxL[, ac(50:150)] <- rlnoise(n = dims(dev_idxB)$iter,
                                        window(dev_idxB, end = 150) %=% 0,
                                        sd = sigmaL_i, b = 0)
      ### insert
      x$oem@deviances$idx$idxB <- dev_idxB
      x$oem@deviances$idx$idxL <- dev_idxL
      
      ### update I_trigger
      I_loss_dev <- apply((x$oem@observations$idx$idxB *
                             dev_idxB)[, ac(50:100)], 6, min)
      I_trigger_dev <- I_loss_dev * 1.4
      x$ctrl$est@args$I_trigger <- c(I_trigger_dev)
      x$I_loss$idx_dev <- I_loss_dev
      
      return(x)
      
    })
    
  }
  
  ## ---------------------------------------------------------------------- ###
  ## run MSE ####
  ## ---------------------------------------------------------------------- ###
  
  ### save fishing history?
  input_i <- lapply(input_i, function(x) {
    x$cut_hist <- cut_hist
    return(x)
  })
  
  ### run MP for each list element
  res <- foreach(input_y = input_i) %do_parts% {
    do.call(mp, input_y)
  }
  names(res) <- names(input_i)
  
  ### ---------------------------------------------------------------------- ###
  ### save ####
  ### ---------------------------------------------------------------------- ###
  
  ### file name
  file_out <- paste0(par_i, collapse = "_")

  
  if (isTRUE(saveMP)) {
    saveRDS(res, file = paste0(path_out, "mp_", file_out, ".rds"))
  }
  
  ## ---------------------------------------------------------------------- ###
  ## stats ####
  ## ---------------------------------------------------------------------- ###
  stats_i <- mp_stats(input = input_i, res_mp = res, 
                      collapse_correction = TRUE)
  saveRDS(stats_i, file = paste0(path_out, file_out, ".rds"))
  
  return(NULL)
  
}

### ------------------------------------------------------------------------ ###
### collate results ####
### ------------------------------------------------------------------------ ###
if (isTRUE(collate)) {
  
  ### collate all available files
  files <- list.files(path = path_out, pattern = "[0-9]*[0-9].rds",
                      full.names = FALSE)
  files <- files[grep(x = files, pattern = "mp", invert = TRUE)]
  names(files) <- sapply(files, function(x) {
    sub(x = x, pattern = ".rds", replacement = "", fixed = TRUE)
  })
  scns <- lapply(files, function(x) {
    params <- sub(x = x, pattern = ".rds", replacement = "", fixed = TRUE)
    params <- an(strsplit(params, split = "_")[[1]])
    names(params) <- names(pars)
    stats <- readRDS(paste0(path_out, x))
    list(pars = params, stats = stats)
  })
  scns[sapply(scns, is.null)] <- NULL
  saveRDS(scns, file = paste0(path_out, "all_runs.rds"))
  
}

### ------------------------------------------------------------------------ ###
### quit ####
### ------------------------------------------------------------------------ ###

quit(save = "no")

