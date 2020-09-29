### ------------------------------------------------------------------------ ###
### run MSE ####
### ------------------------------------------------------------------------ ###

args <- commandArgs(TRUE)
print("arguments passed on to this script:")
print(args)

### evaluate arguments passed to R
for (i in seq_along(args)) eval(parse(text = args[[i]]))

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


input <- hr_par(input = input, brp = brp, lhist = lhist,
                hr = hr, hr_ref = hr_ref, multiplier = multiplier,
                comp_b = comp_b, idxB_lag = idxB_lag, 
                idxB_range_3 = idxB_range_3,
                interval = interval, upper_constraint = upper_constraint,
                lower_constraint = lower_constraint)


### ------------------------------------------------------------------------ ###
### run  ####
### ------------------------------------------------------------------------ ###

res <- do.call(mp, input)

### ------------------------------------------------------------------------ ###
### save ####
### ------------------------------------------------------------------------ ###

### generate file name
file_out <- paste0(c(hr, multiplier, comp_b, idxB_lag, idxB_range_3, interval, 
                     upper_constraint, lower_constraint), collapse = "_")
path_out <- paste0("output/", n_iter, "_", n_yrs, "/", scenario, "/",
                   fhist, "/", paste0(stock, collapse = "_"), "/")
dir.create(path_out, recursive = TRUE)
saveRDS(object = res, file = paste0(path_out, file_out, ".rds"))


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
