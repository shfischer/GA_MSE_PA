# args = c("use_MPI=FALSE", "n_workers=0", "n_blocks=1", "maxiter=1", "popSize=20", "run=1", "stock_id=c(22,26)", "n_iter=10", "n_yrs=50", "opt_step=4")
# cl1 = NULL
### ------------------------------------------------------------------------ ###
### run MSE ####
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
  
  if (!exists("use_MPI")) use_MPI <- FALSE
  if (!exists("n_blocks")) n_blocks <- 1
  if (!exists("n_workers")) n_workers <- 0
  
  if (!exists("popSize")) stop("popSize missing")
  if (!exists("maxiter")) stop("maxiter missing")
  if (!exists("stock_id")) stop("stock_id missing")
  if (!exists("n_iter")) n_iter <- 500
  if (!exists("n_yrs")) n_yrs <- 50
  if (!exists("run")) run <- maxiter
  if (!exists("opt_step")) opt_step <- 4
  if (!exists("collate")) collate <- TRUE
  
} else {
  
  stop("no argument passed to R")
  
}

### ------------------------------------------------------------------------ ###
### set up environment ####
### ------------------------------------------------------------------------ ###

### load packages
req_pckgs <- c("FLCore", "FLash", "mseDL", "GA", "doParallel", "doRNG", "FLBRP")
for (i in req_pckgs) library(package = i, character.only = TRUE)

### load additional functions
source("funs.R")
source("GA_funs.R")

### ------------------------------------------------------------------------ ###
### setup parallel environment ####
### ------------------------------------------------------------------------ ###

### hybrid MPI
if (isTRUE(use_MPI)) {
  ### 1st: doMPI cluster with 1 worker per node
  library(doMPI)
  cl1 <- startMPIcluster()
  print(cl1)
  registerDoMPI(cl1)
  cl_length_1 <- cl1$workerCount
  cl_length_1
  
  ### 2nd: doParallel workers inside doMPI workers
  foreach(i = seq(cl_length_1)) %dopar% {
    ### load packages and functions into MPI workers
    for (i in req_pckgs) library(package = i, character.only = TRUE,
                                 warn.conflicts = FALSE, verbose = FALSE,
                                 quietly = TRUE)
    source("funs.R", echo = FALSE)
    source("GA_funs.R", echo = FALSE)
    ### start doParallel
    if (isTRUE(n_workers > 1)) {
      cl2 <- makeCluster(n_workers)
      registerDoParallel(cl2)
      cl_length_2 <- length(cl2)
      ### load packages and functions into parallel workers
      . <- foreach(i = seq(cl_length_2)) %dopar% {
        for (i in req_pckgs) library(package = i, character.only = TRUE,
                                     warn.conflicts = FALSE, verbose = FALSE,
                                     quietly = TRUE)
        source("funs.R", echo = FALSE)
        source("GA_funs.R", echo = FALSE)
      }
    }
  }
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
      source("GA_funs.R", echo = FALSE)
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
  readRDS(paste0("input/", n_iter, "_", n_yrs, "/OM_2_mp_input/", x,
                        ".rds"))
})

### ------------------------------------------------------------------------ ###
### specify scenario ####
### ------------------------------------------------------------------------ ###

if (isTRUE(opt_step %in% 1:2)) {
  scenario <- "SSB_idx_r_only"
  input <- lapply(input, function(x) {
    ### OEM
    x$oem@args$idx_dev <- FALSE
    x$oem@args$ssb <- TRUE
    x$oem@args$lngth <- FALSE
    x$oem@args$lngth_dev <- FALSE
    ### catch rule components
    x$ctrl.mp$ctrl.est@args$comp_r <- TRUE
    x$ctrl.mp$ctrl.est@args$comp_f <- FALSE
    x$ctrl.mp$ctrl.est@args$comp_b <- FALSE
    return(x)
  })
  ### GA arguments
  ga_suggestions = rbind(c(0, 1, 1, 1),
                         c(1, 1, 1, 1),
                         c(1, 2, 3, 1))
  ga_lower = c(0, 1, 1, 1)
  ga_upper = c(5, 5, 5, 5)
  ga_names = c("lag_idx", "range_idx_1", "range_idx_2",
               "range_catch")
  if (isTRUE(opt_step == 2)) {
    ### include index error
    scenario <- "SSB_idx_r_only_error"
    input <- lapply(input, function(x) {
      x$oem@args$idx_dev <- TRUE
      return(x)
    })
  }
} else if (isTRUE(opt_step == 3)) {
  scenario <- "SSB_idx_rfb_exp"
  input <- lapply(input, function(x) {
    ### OEM: no uncertainty
    x$oem@args$idx_dev <- FALSE
    x$oem@args$ssb <- TRUE
    x$oem@args$lngth <- TRUE
    x$oem@args$lngth_dev <- FALSE
    ### catch rule components
    x$ctrl.mp$ctrl.est@args$comp_r <- TRUE
    x$ctrl.mp$ctrl.est@args$comp_f <- TRUE
    x$ctrl.mp$ctrl.est@args$comp_b <- TRUE
    ### catch lag fixed
    x$ctrl.mp$ctrl.est@args$catch_lag <- 1
    return(x)
  })
  ### GA arguments
  ga_names = c("lag_idx", "range_idx_1", "range_idx_2", "range_catch",
               "exp_r", "exp_f", "exp_b")
  ga_lower = c(0, 1, 1, 1, 0, 0, 0)
  ga_upper = c(1, 5, 5, 5, 2, 2, 2)
  ga_suggestions = rbind(c(0, 1, 1, 1, 1, 1, 1), ### most current data
                         c(0, 1, 1, 1, 0, 0, 0), ### constant catch
                         c(1, 2, 3, 1, 1, 1, 1), ### default
                         c(1, 2, 3, 1, 0, 1, 1), ### f*b
                         c(1, 2, 3, 1, 0, 1, 0), ### f
                         c(1, 2, 3, 1, 0, 0, 1)) ### b
} else if (isTRUE(opt_step %in% c(4, 5, 6))) {
  scenario <- switch(ac(opt_step),
    "4" = "SSB_idx_rfb_exp_error",
    "5" = "SSB_idx_rfb_exp_error_Itrigger",
    "6" = "full"
  )
  input <- lapply(input, function(x) {
    ### OEM: activate uncertainty
    x$oem@args$idx_dev <- TRUE
    x$oem@args$ssb <- TRUE
    x$oem@args$lngth <- TRUE
    x$oem@args$lngth_dev <- TRUE
    ### IEM: do not activate uncertainty
    x$iem@args$use_dev <- FALSE
    ### catch rule components
    x$ctrl.mp$ctrl.est@args$comp_r <- TRUE
    x$ctrl.mp$ctrl.est@args$comp_f <- TRUE
    x$ctrl.mp$ctrl.est@args$comp_b <- TRUE
    ### catch lag fixed
    x$ctrl.mp$ctrl.est@args$catch_lag <- 1
    return(x)
  })
  ### GA arguments
  if (isTRUE(opt_step %in% c(4, 5))) {
    ga_names = c("lag_idx", "range_idx_1", "range_idx_2", "range_catch",
                 "exp_r", "exp_f", "exp_b")
    ga_lower = c(0, 1, 1, 1, 0, 0, 0)
    ga_upper = c(1, 5, 5, 5, 2, 2, 2)
    ga_suggestions = rbind(c(0, 1, 1, 1, 1, 1, 1), ### most current data
                           c(0, 1, 1, 1, 0, 0, 0), ### constant catch
                           c(1, 2, 3, 1, 1, 1, 1), ### default
                           c(1, 2, 3, 1, 0, 1, 1), ### f*b
                           c(1, 2, 3, 1, 0, 1, 0), ### f
                           c(1, 2, 3, 1, 0, 0, 1)) ### b
  } else if (isTRUE(opt_step %in% c(6))) {
    ga_names = c("lag_idx", "range_idx_1", "range_idx_2", "range_catch",
                 "exp_r", "exp_f", "exp_b", "interval")
    ga_lower = c(0, 1, 1, 1, 0, 0, 0, 1)
    ga_upper = c(1, 5, 5, 5, 2, 2, 2, 5)
    ga_suggestions = rbind(c(0, 1, 1, 1, 1, 1, 1, 1), ### most current data
                           c(0, 1, 1, 1, 1, 1, 1, 2),
                           c(0, 1, 1, 1, 0, 0, 0, 2), ### constant catch
                           ### default, annual/biennial, turning off elements
                           expand.grid(1, 2, 3, 1, 0:1, 0:1, 0:1, 1:2))
  }
  if (isTRUE(opt_step %in% c(5))) {
    ### use Iloss
    input <- lapply(input, function(x) {
      x$ctrl.mp$ctrl.est@args$I_trigger[] <- c(input$I_loss$SSB_idx_dev) * 1.4
      return(x)
    })
  }
}
### output path
path_out <- paste0("output/", n_iter, "_", n_yrs, "/ms/", scenario, "/",
                   paste0(stock, collapse = "_"), "/")
dir.create(path_out, recursive = TRUE)

### ------------------------------------------------------------------------ ###
### run MSE with GA ####
### ------------------------------------------------------------------------ ###

### within scenario parallelisation?
if (isTRUE(n_workers > 1) & isTRUE(n_blocks > 1)) {
  ### use Iloss
  input <- lapply(input, function(x) {
    x$genArgs$nblocks <- n_blocks
    return(x)
  })
}

### set random seed for reproducibility
registerDoRNG(123)
set.seed(1)

### run GA
system.time({
  res <- ga(type = "real-valued", fitness = mse_ms, input = input,
            path = path_out, check_file = TRUE,
            scenario = scenario,
            suggestions = ga_suggestions, lower = ga_lower, upper = ga_upper,
            names = ga_names,
            maxiter = maxiter, popSize = popSize, run = run,
            monitor = TRUE, keepBest = TRUE, parallel = cl1, seed = 1)
})

# debugonce(mse_r)
# mse_r(c(0, 1, 1, 0, 1, 1), input = input, path = path_out, check_file = TRUE,
#       scenario = "SSB_idx_r")


### save result
saveRDS(object = res, file = paste0(path_out, "res.rds"))

### ------------------------------------------------------------------------ ###
### collate runs ####
### ------------------------------------------------------------------------ ###

if (isTRUE(collate)) {
  files <- list.files(path = path_out, pattern = "[0-9]*[0-9].rds",
                      full.names = FALSE)
  scns <- lapply(files, function(x) {
    pars <- an(strsplit(sub(x = x, pattern = ".rds", replacement = "", fixed = TRUE), split = "_")[[1]])
    names(pars) <- ga_names
    stats <- readRDS(paste0(path_out, x))
    list(pars = pars, obj = stats$obj, stats = stats$stats)
  })
  scns <- scns[order(sapply(scns, "[[", "obj"), decreasing = TRUE)]
  saveRDS(scns, file = paste0(path_out, "runs.rds"))
}

### ------------------------------------------------------------------------ ###
### quit ####
### ------------------------------------------------------------------------ ###

quit(save = "no")

