# args = c("use_MPI=FALSE", "n_workers=10", "n_blocks=10", "maxiter=1", "popSize=35", "run=1", "stock_id=c(22,26)", "n_iter=10", "n_yrs=50", "multiplier=TRUE")
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
  ### set default arguments
  ### parallelization
  if (!exists("use_MPI")) use_MPI <- FALSE
  if (!exists("n_blocks")) n_blocks <- 1
  if (!exists("n_workers")) n_workers <- 0
  ### scenario definition
  if (!exists("n_iter")) n_iter <- 500
  if (!exists("n_yrs")) n_yrs <- 50
  if (!exists("fhist")) fhist <- "one-way"
  if (!exists("catch_rule")) catch_rule <- "catch_rule"
  ### GA search
  if (!exists("ga_search")) ga_search <- TRUE
  if (isTRUE(ga_search)) {
    if (!exists("popSize")) stop("popSize missing")
    if (!exists("maxiter")) stop("maxiter missing")
    if (!exists("stock_id")) stop("stock_id missing")
    if (!exists("run")) run <- maxiter
    if (!exists("collate")) collate <- TRUE
    ### objective function elements
    if (!exists("obj_SSB")) obj_SSB <- TRUE
    if (!exists("obj_F")) obj_F <- FALSE
    if (!exists("obj_C")) obj_C <- TRUE
    if (!exists("obj_risk")) obj_risk <- TRUE
    if (!exists("obj_ICV")) obj_ICV <- TRUE
  }

} else {
  
  stop("no argument passed to R")
  
}

### ------------------------------------------------------------------------ ###
### set up environment ####
### ------------------------------------------------------------------------ ###

### load packages
## GA fork from GitHub devtools::install_github("shfischer/GA")
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
  . <- foreach(i = seq(cl_length_1)) %dopar% {
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
  readRDS(paste0("input/", n_iter, "_", n_yrs, "/OM_2_mp_input/", fhist, "/", x,
                        ".rds"))
})


### ------------------------------------------------------------------------ ###
### specify scenario ####
### ------------------------------------------------------------------------ ###

### default catch rule
input <- lapply(input, function(x) {
  ### OEM: activate uncertainty
  x$oem@args$idx_dev <- TRUE
  x$oem@args$ssb <- FALSE
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

### default ICES rule: 2 over 3
if (isTRUE(catch_rule == "2over3")) {
  input <- lapply(input, function(x) {
    ### OEM: turn of length index
    x$oem@args$lngth <- FALSE
    x$oem@args$lngth_dev <- FALSE
    ### add PA buffer stock status and deviation
    x$oem@args$PA_status <- TRUE
    x$oem@args$PA_status_dev <- TRUE
    ### catch rule components: turn of f & b, 2 over 3 rule
    x$ctrl.mp$ctrl.est@args$comp_f <- FALSE
    x$ctrl.mp$ctrl.est@args$comp_b <- FALSE
    x$ctrl.mp$ctrl.est@args$idxB_lag <- 1
    x$ctrl.mp$ctrl.est@args$idxB_range_1 <- 2
    x$ctrl.mp$ctrl.est@args$idxB_range_2 <- 3
    ### PA buffer
    x$ctrl.mp$ctrl.est@args$pa_buffer <- TRUE
    ###
    x$ctrl.mp$ctrl.phcr@args$exp_r <- 1
    x$ctrl.mp$ctrl.phcr@args$exp_f <- 0
    x$ctrl.mp$ctrl.phcr@args$exp_b <- 1 ### PA buffer
    ### biennial
    #x$ctrl.mp$ctrl.hcr@args$interval <- 2
    ### uncertainty cap
    x$ctrl.mp$ctrl.is@args$upper_constraint <- 1.2
    x$ctrl.mp$ctrl.is@args$lower_constraint <- 0.8
    return(x)
  })
  # input$pol$oem@method <- wklife_3.2.1_obs
  # input$pol$ctrl.mp$ctrl.est@method <- wklife_3.2.1_est
  # input$pol$ctrl.mp$ctrl.is@method <- is_r
  #debugonce(goFishDL)
  #res <- do.call(mpDL, c(input$pol, cut_hist = FALSE))
  
}

### within scenario parallelisation?
if (isTRUE(n_workers > 1) & isTRUE(n_blocks > 1)) {
  ### use Iloss
  input <- lapply(input, function(x) {
    x$genArgs$nblocks <- n_blocks
    return(x)
  })
}

### ------------------------------------------------------------------------ ###
### GA set-up ####
### ------------------------------------------------------------------------ ###
if (isTRUE(catch_rule == "catch_rule") & isTRUE(ga_search)) {

  ### GA arguments
  ga_names <- c("lag_idx", "range_idx_1", "range_idx_2", "range_catch",
               "exp_r", "exp_f", "exp_b", "interval", "multiplier")
  ga_suggestions <- rbind(c(0, 1, 1, 1, 1, 1, 1, 1, 1), ### most current data
                          c(0, 1, 1, 1, 1, 1, 1, 2, 1),
                          c(0, 1, 1, 1, 0, 0, 0, 2, 1), ### constant catch
                          ### default, annual/biennial, turning off elements
                          expand.grid(1, 2, 3, 1, 0:1, 0:1, 0:1, 1:2, 0:1))
  ga_default <- c(1, 2, 3, 1, 1, 1, 1, 2, 1)
  ga_lower <- c(0, 1, 1, 1, 0, 0, 0, 1, 0)
  ga_upper <- c(1, 5, 5, 1, 2, 2, 2, 5, 2)
  ### turn of parameters not requested, i.e. limit to default value
  pos_default <- which(!unlist(mget(ga_names, ifnotfound = FALSE)))
  ga_lower[pos_default] <- ga_default[pos_default]
  ga_upper[pos_default] <- ga_default[pos_default]
  ### remove not requested parameters from suggestions
  ga_suggestions[, pos_default] <- rep(ga_default[pos_default], 
                                       each = nrow(ga_suggestions))
  ga_suggestions <- unique(ga_suggestions)
  
  ### ---------------------------------------------------------------------- ###
  ### paths ####
  ### ---------------------------------------------------------------------- ###
  
  ### output path
  # ### set name depending on which GA parameters are used
  scn_pars <- paste0(ga_names[setdiff(seq_along(ga_names), pos_default)],
                     collapse = "-")
  scenario <- "trial"
  path_out <- paste0("output/", n_iter, "_", n_yrs, "/ms/", scenario, "/",
                     fhist, "/",
                     paste0(stock, collapse = "_"), "/")
  dir.create(path_out, recursive = TRUE)
  
  ### objective function elements
  obj_fun_elements <- c("SSB", "F", "C", "risk", "ICV")
  obj_desc <- obj_fun_elements[c(obj_SSB, obj_F, obj_C, obj_risk, obj_ICV)]
  obj_desc <- paste0("obj_", paste0(obj_desc, collapse = "_"), collapse = "")
  
  ### store input data in temp file
  inp_file <- tempfile()
  saveRDS(object = input, file = inp_file, compress = FALSE)
  rm(input)
  gc()
  
  ### ---------------------------------------------------------------------- ###
  ### run MSE with GA ####
  ### ---------------------------------------------------------------------- ###
  
  ### set random seed for reproducibility
  registerDoRNG(123)
  set.seed(1)
  
  ### run GA
  system.time({
    res <- ga(type = "real-valued", fitness = mse_ms, inp_file = inp_file,
              obj_SSB = obj_SSB, obj_F = obj_F, obj_C = obj_C, 
              obj_risk = obj_risk, obj_ICV = obj_ICV,
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
  saveRDS(object = res, file = paste0(path_out, scn_pars, 
                                      "--", obj_desc, "_res.rds"))
  
  ### ---------------------------------------------------------------------- ###
  ### collate runs ####
  ### ---------------------------------------------------------------------- ###
  
  if (isTRUE(collate)) {
    files <- list.files(path = path_out, pattern = "[0-9]*[0-9].rds",
                        full.names = FALSE)
    names(files) <- sapply(files, function(x) {
      sub(x = x, pattern = ".rds", replacement = "", fixed = TRUE)
    })
    scns <- lapply(files, function(x) {
      pars <- an(strsplit(sub(x = x, pattern = ".rds", replacement = "", fixed = TRUE), split = "_")[[1]])
      names(pars) <- ga_names
      ### only keep scenarios where requested parameters are changed
      if (!all(ga_default[pos_default] == pars[pos_default])) return(NULL)
      stats <- readRDS(paste0(path_out, x))
      list(pars = pars, stats = stats)
    })
    scns[sapply(scns, is.null)] <- NULL
    #scns <- scns[order(sapply(scns, "[[", "obj"), decreasing = TRUE)]
    saveRDS(scns, file = paste0(path_out, scn_pars, "--", obj_desc, "_runs.rds"))
  }

### other catch rules
} else if (isTRUE(catch_rule == "2over3")) {
  
  ### output path
  path_out <- paste0("output/", n_iter, "_", n_yrs, "/ms/", catch_rule, "/",
                     fhist, "/")
  dir.create(path_out, recursive = TRUE)
  
  ### run MSE
  ### run MP for each list element
  res_mp <- lapply(input, function(x) {
    do.call(mpDL, x)
  })
  saveRDS(res_mp, paste0(path_out, paste0(stock, collapse = "_"), "_mp.rds"))
  
  ### stats
  stats <- mp_stats(input = input, res_mp = res_mp)
  saveRDS(stats, paste0(path_out, paste0(stock, collapse = "_"), "_stats.rds"))
  
}

### ------------------------------------------------------------------------ ###
### quit ####
### ------------------------------------------------------------------------ ###

quit(save = "no")

