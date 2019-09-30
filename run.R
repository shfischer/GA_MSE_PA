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
  
  ### parallelisation environment
  if (!exists("par_env")) par_env <- 1
  if (!exists("n_workers")) n_workers <- 1
  if (!exists("popSize")) stop("popSize missing")
  if (!exists("maxiter")) stop("maxiter missing")
  if (!exists("stock_id")) stop("stock_id missing")
  if (!exists("n_iter")) n_iter <- 500
  if (!exists("n_yrs")) n_yrs <- 50
  if (!exists("run")) run <- maxiter
  if (!exists("opt_step")) opt_step <- 1
  if (!exists("collate")) collate <- 1
  if (!exists("Itrigger")) Itrigger <- 0
  
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
### par_env=1 -> MPI (Rmpi, DoMPI)
### par_env=2 -> DoParallel

if (par_env == 1) {
  
  library(doMPI)
  cl <- startMPIcluster()
  print(cl)
  registerDoMPI(cl)
  cl_length <- cl$workerCount
  cl_length
  
} else if (par_env == 2) {
  
  library(doParallel)
  cl <- makeCluster(n_workers)
  print(cl)
  registerDoParallel(cl)
  cl_length <- length(cl)
  
}

### load packages and functions into workers
. <- foreach(i = seq(cl_length)) %dopar% {
  #devtools::load_all("../mse/")
  for (i in req_pckgs) library(package = i, character.only = TRUE)
  source("funs.R")
  source("GA_funs.R")
}

### ------------------------------------------------------------------------ ###
### load data ####
### ------------------------------------------------------------------------ ###

stocks <- read.csv("input/stocks.csv", stringsAsFactors = FALSE)
stock <- stocks$stock[stock_id]
input <- readRDS(paste0("input/", n_iter, "_", n_yrs, "/OM_2_mp_input/", stock, 
                        ".rds"))

if (isTRUE(opt_step %in% 1:2)) {
  ### OEM
  input$oem@args$idx_dev <- FALSE
  input$oem@args$ssb <- TRUE
  input$oem@args$lngth <- FALSE
  input$oem@args$lngth_dev <- FALSE
  
  scenario <- "SSB_idx_r_only"
  ### catch rule components
  input$ctrl.mp$ctrl.est@args$comp_r <- TRUE
  input$ctrl.mp$ctrl.est@args$comp_f <- FALSE
  input$ctrl.mp$ctrl.est@args$comp_b <- FALSE
  ### GA arguments
  ga_suggestions = rbind(c(0, 1, 1, 1),
                         c(1, 1, 1, 1),
                         c(1, 2, 3, 1))
  ga_lower = c(0, 1, 1, 1)
  ga_upper = c(5, 5, 5, 5)
  ga_names = c("lag_idx", "range_idx_1", "range_idx_2",
               "range_catch")
  if (isTRUE(opt_step == 2)) {
    scenario <- "SSB_idx_r_only_error"
    input$oem@args$idx_dev <- TRUE
  }
} else if (isTRUE(opt_step == 3)) {
  scenario <- "SSB_idx_rfb_exp"
  ### OEM: no uncertainty
  input$oem@args$idx_dev <- FALSE
  input$oem@args$ssb <- TRUE
  input$oem@args$lngth <- TRUE
  input$oem@args$lngth_dev <- FALSE
  ### catch rule components
  input$ctrl.mp$ctrl.est@args$comp_r <- TRUE
  input$ctrl.mp$ctrl.est@args$comp_f <- TRUE
  input$ctrl.mp$ctrl.est@args$comp_b <- TRUE
  ### catch lag fixed
  input$ctrl.mp$ctrl.est@args$catch_lag <- 1
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
} else if (isTRUE(opt_step %in% c(4, 5))) {
  scenario <- "SSB_idx_rfb_exp_error"
  ### OEM: activate uncertainty
  input$oem@args$idx_dev <- TRUE
  input$oem@args$ssb <- TRUE
  input$oem@args$lngth <- TRUE
  input$oem@args$lngth_dev <- TRUE
  ### IEM: do not activate uncertainty
  input$iem@args$use_dev <- FALSE
  ### catch rule components
  input$ctrl.mp$ctrl.est@args$comp_r <- TRUE
  input$ctrl.mp$ctrl.est@args$comp_f <- TRUE
  input$ctrl.mp$ctrl.est@args$comp_b <- TRUE
  ### catch lag fixed
  input$ctrl.mp$ctrl.est@args$catch_lag <- 1
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
  if (isTRUE(opt_step %in% c(5))) {
    scenario <- "SSB_idx_rfb_exp_error_Itrigger"
    ### use Iloss
    input$ctrl.mp$ctrl.est@args$I_trigger[] <- c(input$I_loss$SSB_idx_dev) * 1.4
  }
}
path_out <- paste0("output/", n_iter, "_", n_yrs, "/", scenario, "/", 
                   stock, "/")
dir.create(path_out, recursive = TRUE)

### ------------------------------------------------------------------------ ###
### run MSE with GA ####
### ------------------------------------------------------------------------ ###

### set random seed for reproducibility
registerDoRNG(123)
set.seed(1)

### run GA
system.time({
  res <- ga(type = "real-valued", fitness = mse_r, input = input,
            path = path_out, check_file = TRUE,
            scenario = scenario,
            suggestions = ga_suggestions, lower = ga_lower, upper = ga_upper, 
            names = ga_names,
            maxiter = maxiter, popSize = popSize, run = run,
            monitor = TRUE, keepBest = TRUE, parallel = cl, seed = 1)
})

# debugonce(mse_r)
# mse_r(c(0, 1, 1, 0, 1, 1), input = input, path = path_out, check_file = TRUE,
#       scenario = "SSB_idx_r")


### save result
saveRDS(object = res, file = paste0(path_out, "res.rds"))

### ------------------------------------------------------------------------ ###
### collate runs ####
### ------------------------------------------------------------------------ ###

if (isTRUE(all.equal(collate, 1))) {
  files <- list.files(path = path_out, pattern = "[0-9]*[0-9].rds",
                      full.names = FALSE)
  scns <- lapply(files, function(x) {
    pars <- an(strsplit(sub(x = x, pattern = ".rds", replacement = "", fixed = TRUE), split = "_")[[1]])
    names(pars) <- ga_names
    stats <- readRDS(paste0(path_out, x))
    mtime <- file.mtime(paste0(path_out, x))
    c(stock = stock, file = x, pars, stats, mtime = mtime)
  })
  scns <- do.call(rbind, scns)
  scns <- as.data.frame(scns, stringsAsFactors = FALSE)
  for (i in seq(ncol(scns))) scns[, i] <- unlist(scns[, i])
  scns <- scns[order(scns$fitness, decreasing = TRUE), ]
  saveRDS(scns, file = paste0(path_out, "runs.rds"))
  ### add GA results to summary
  path_smry <- paste0("output/", n_iter, "_", n_yrs, "/", scenario, "/")
  if (isTRUE(scenario %in% c("SSB_idx_r_only", "SSB_idx_r_only_error"))) {
    rnd_dig <- round(res@solution)
  } else if (isTRUE(scenario %in% c("SSB_idx_rfb_exp", 
                                    "SSB_idx_rfb_exp_error"))) {
    rnd_dig <- res@solution
    if (isTRUE(nrow(rnd_dig) > 1)) rnd_dig <- rnd_dig[1,]
    rnd_dig[1:4] <- round(rnd_dig[1:4])
    rnd_dig[5:7] <- round(rnd_dig[5:7], 1)
  }
  smry_add <- scns[1,]
  if (file.exists(paste0(path_smry, "summary.csv"))) {
    smry <- read.csv(file = paste0(path_smry, "summary.csv"), 
                     stringsAsFactors = FALSE)
  } else {
    smry <- smry_add[0, ]
  }
  if (isTRUE(stock %in% smry$stock)) {
    smry[smry$stock == stock, ] <- smry_add
  } else {
    smry <- rbind(smry, smry_add)
  }
  write.csv(smry, file = paste0(path_smry, "summary.csv"), row.names = FALSE)
}

### ------------------------------------------------------------------------ ###
### quit ####
### ------------------------------------------------------------------------ ###

quit(save = "no")

