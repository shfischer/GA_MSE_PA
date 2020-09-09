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
  if (!exists("n_yrs")) n_yrs <- 100
  if (!exists("fhist")) fhist <- "random"
  if (!exists("catch_rule")) catch_rule <- "catch_rule"
  if (!exists("scenario")) scenario <- "risk"
  if (!exists("collate")) collate <- TRUE
  if (!exists("saveMP")) saveMP <- TRUE
  
  if (!exists("n_parts")) n_parts <- 1
  if (!exists("part")) part <- 1

} else {
  
  stop("no argument passed to R")
  
}

### ------------------------------------------------------------------------ ###
### set up environment ####
### ------------------------------------------------------------------------ ###

### load packages
### GA fork from GitHub remotes::install_github("shfischer/GA")
### use mse fork from shfischer/mse, branch mseDL2.0 
### remotes::install_github("shfischer/mse", ref = "mseDL2.0)
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
### specify scenario ####
### ------------------------------------------------------------------------ ###

### within scenario parallelisation?
if (isTRUE(n_workers > 1) & isTRUE(n_blocks > 1)) {
  input <- lapply(input, function(x) {
    x$args$nblocks <- n_blocks
    return(x)
  })
}

### ------------------------------------------------------------------------ ###
### GA set-up ####
### ------------------------------------------------------------------------ ###

### GA arguments
ga_names <- c("lag_idx", "range_idx_1", "range_idx_2", "range_catch",
             "exp_r", "exp_f", "exp_b", "interval", "multiplier",
             "upper_constraint", "lower_constraint")
params <- c(1, 2, 3, 1, 1, 1, 1, 2, 1, Inf, 0)
### fix parameters?
pos_fixed <- which(sapply(mget(ga_names, ifnotfound = FALSE), is.numeric))
par_fixed <- names(pos_fixed)
val_fixed <- mget(par_names, ifnotfound = FALSE)[pos_fixed]
params[pos_fixed] <- unlist(val_fixed)

### ------------------------------------------------------------------------ ###
### set catch rule parameters ####
### ------------------------------------------------------------------------ ###

### insert arguments into input object for mp
input <- lapply(input, function(x) {
  x$ctrl$est@args$idxB_lag     <- params[1]
  x$ctrl$est@args$idxB_range_1 <- params[2]
  x$ctrl$est@args$idxB_range_2 <- params[3]
  x$ctrl$est@args$catch_range  <- params[4]
  x$ctrl$est@args$comp_m <- params[9]
  x$ctrl$phcr@args$exp_r <- params[5]
  x$ctrl$phcr@args$exp_f <- params[6]
  x$ctrl$phcr@args$exp_b <- params[7]
  x$ctrl$hcr@args$interval <- params[8]
  x$ctrl$isys@args$upper_constraint <- params[10]
  x$ctrl$isys@args$lower_constraint <- params[11]
  return(x)
})
  
  
### ---------------------------------------------------------------------- ###
### paths ####
### ---------------------------------------------------------------------- ###

### output path
path_out <- paste0("output/", n_iter, "_", n_yrs, "/", scenario, "/",
                   fhist, "/",
                   paste0(stock, collapse = "_"), "/")
dir.create(path_out, recursive = TRUE)


### ---------------------------------------------------------------------- ###
### run MSE ####
### ---------------------------------------------------------------------- ###

### set random seed for reproducibility
registerDoRNG(123)
set.seed(1)

### run MP for each list element
res_mp <- lapply(input, function(x) {
  gc()
  if (getDoParWorkers() > 1)
    . <- foreach(i = 1:getDoParWorkers()) %dopar% {invisible(gc())}
  do.call(mp, x)
})
gc()


### save result
saveRDS(object = res_mp, file = paste0(path_out, paste0(params, collapse = "_"),
                                    "_", part, "-", n_parts, ".rds"))

### ------------------------------------------------------------------------ ###
### quit ####
### ------------------------------------------------------------------------ ###

quit(save = "no")

