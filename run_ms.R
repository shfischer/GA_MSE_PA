### ------------------------------------------------------------------------ ###
### run MSE ####
### ------------------------------------------------------------------------ ###

args <- commandArgs(TRUE)
print("arguments passed on to this script:")
print(args)

### evaluate arguments passed to R
for (i in seq_along(args)) eval(parse(text = args[[i]]))

if (!exists("saveMP")) saveMP <- TRUE
if (!exists("stats")) stats <- TRUE
if (!exists("collate")) collate <- FALSE

### ------------------------------------------------------------------------ ###
### set up environment ####
### ------------------------------------------------------------------------ ###

### load packages
req_pckgs <- c("mse", "tidyr", "dplyr", "doParallel")
for (i in req_pckgs) library(package = i, character.only = TRUE)

### load additional functions
source("funs.R")

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
  }
}

### ------------------------------------------------------------------------ ###
### load input  ####
### ------------------------------------------------------------------------ ###

### stock list
stocks <- read.csv("input/stocks.csv", stringsAsFactors = FALSE)
stock <- stocks$stock[stock_id]

### path to input files
path_in <- paste0("input/", n_iter, "_", n_yrs, "/OM_2_mp_input/", fhist, "/")
input <- readRDS(paste0(path_in, stock, ".rds"))

input$args$nblocks <- n_blocks

### ------------------------------------------------------------------------ ###
### MP parameters ####
### ------------------------------------------------------------------------ ###

### default values
if (!exists("hr")) hr <- "uniform"
if (!exists("multiplier")) multiplier <- 1
if (!exists("comp_b")) comp_b <- FALSE
if (!exists("interval")) interval <- 1
if (!exists("idxB_lag")) idxB_lag <- 1
if (!exists("idxB_range_3")) idxB_range_3 <- 1
if (!exists("upper_constraint")) upper_constraint <- Inf
if (!exists("lower_constraint")) lower_constraint <- 0

### load reference values
hr_ref <- readRDS("input/catch_rates.rds")[[stock]]
brp <- readRDS("input/brps.rds")[[stock]]
lhist <- stocks[stocks$stock == stock, ]

### HR rule parameters
hr_params <- data.frame(multiplier = multiplier,
                        comp_b = comp_b,
                        idxB_lag = idxB_lag,
                        idxB_range_3 = idxB_range_3,
                        interval = interval,
                        upper_constraint = upper_constraint,
                        lower_constraint = lower_constraint)


if (isTRUE(n_workers > 1 & n_blocks == 1)) {
  `%do_tmp%` <- `%dopar%`
} else {
  `%do_tmp%` <- `%do%`
}

. <- foreach(hr_i = seq(nrow(hr_params))) %do_tmp% {
  
  par_i <- hr_params[hr_i, ]
  
  input <- hr_par(input = input, brp = brp, lhist = lhist,
                  hr = hr, hr_ref = hr_ref, 
                  multiplier = par_i$multiplier,
                  comp_b = par_i$comp_b, idxB_lag = par_i$idxB_lag, 
                  idxB_range_3 = par_i$idxB_range_3,
                  interval = par_i$interval, 
                  upper_constraint = par_i$upper_constraint,
                  lower_constraint = par_i$lower_constraint)
  
  ### ------------------------------------------------------------------------ ###
  ### run  ####
  ### ------------------------------------------------------------------------ ###
  
  res <- do.call(mp, input)
  
  ### ------------------------------------------------------------------------ ###
  ### save ####
  ### ------------------------------------------------------------------------ ###
  
  ### generate file name
  file_out <- paste0(c(hr, par_i$multiplier, par_i$comp_b, par_i$idxB_lag, 
                       par_i$idxB_range_3, par_i$interval, 
                       par_i$upper_constraint, par_i$lower_constraint), 
                     collapse = "_")
  path_out <- paste0("output/", n_iter, "_", n_yrs, "/", scenario, "/",
                     fhist, "/", paste0(stock, collapse = "_"), "/")
  dir.create(path_out, recursive = TRUE)
  if (isTRUE(saveMP))
    saveRDS(object = res, file = paste0(path_out, "mp_", file_out, ".rds"))
  
  ### ---------------------------------------------------------------------- ###
  ### stats ####
  ### ---------------------------------------------------------------------- ###
  
  if (isTRUE(stats)) {
    res_stats <- mp_stats(input = list(input), res_mp = list(res), 
                          collapse_correction = TRUE)
    res_stats <- cbind(stock = stock, par_i, t(res_stats))
    saveRDS(object = res_stats, 
            file = paste0(path_out, "stats_", file_out, ".rds"))
  }

}

### ------------------------------------------------------------------------ ###
### collate stats ####
### ------------------------------------------------------------------------ ###

if (isTRUE(stats) & isTRUE(collate) & isTRUE(nrow(hr_params) > 1)) {
  files <- paste0("stats_", hr, "_", 
                  sapply(seq(nrow(hr_params)), 
                         function(x) paste0(hr_params[x,], collapse = "_")),
                  ".rds")
  files <- paste0("output/", n_iter, "_", n_yrs, "/", scenario, "/",
                  fhist, "/", paste0(stock, collapse = "_"), "/",
                  files)
  stats_all <- lapply(files, readRDS)
  stats_all <- do.call(rbind, stats_all)
  
  saveRDS(stats_all, file = paste0(
    "output/", n_iter, "_", n_yrs, "/", scenario, "/", fhist, "/", 
    paste0(stock, collapse = "_"), "/",
    "collated_stats_", hr, "_", 
    paste0(apply(hr_params, 2, function(x) {
      ifelse(isTRUE(length(unique(x)) > 1), paste0(range(x), collapse = "-"), x[1])
    }), collapse = "_"), ".rds"))
  
}

### ------------------------------------------------------------------------ ###
### quit ####
### ------------------------------------------------------------------------ ###

quit(save = "no")



# input <- readRDS("input/10000_100/OM_2_mp_input/random/bll.rds")
# debugonce(wklife_3.2.1_est)
# debugonce(wklife_3.2.1_obs)
# debugonce(input$ctrl.mp$ctrl.hcr@method)
# debugonce(mpDL)
# debugonce(goFishDL)
# input$genArgs$nblocks = 10
# input$cut_hist = FALSE
# res <- do.call(mpDL, input)
# 
# ### timing
# system.time({res1 <- do.call(mpDL, input)})
# path_out <- paste0("output/", n_iter, "_", yrs_proj, "/", fhist, "/")
# dir.create(path_out, recursive = TRUE)
# saveRDS(res1, file = paste0(path_out, stock, ".rds"))
