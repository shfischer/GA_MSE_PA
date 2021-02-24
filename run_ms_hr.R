### ------------------------------------------------------------------------ ###
### run MSE ####
### ------------------------------------------------------------------------ ###

args <- commandArgs(TRUE)
print("arguments passed on to this script:")
print(args)

### evaluate arguments passed to R
for (i in seq_along(args)) eval(parse(text = args[[i]]))

### evaluate arguments, if they are passed to R:
if (length(args) > 0) {
  
  ### parallelization
  if (!exists("use_MPI")) use_MPI <- FALSE
  if (!exists("n_blocks")) n_blocks <- 1
  if (!exists("n_workers")) n_workers <- 0
  
  ### split OM into blocks?
  if (!exists("n_parts")) n_parts <- 1
  if (!exists("part")) part <- 1
  
  ### projection details
  if (!exists("n_iter")) n_iter <- 500
  if (!exists("n_yrs")) n_yrs <- 50
  if (!exists("fhist")) fhist <- "random"
  
  ### MP parameters
  if (!exists("hr")) hr <- "length"
  if (!exists("multiplier")) multiplier <- 1
  if (!exists("comp_b")) comp_b <- TRUE
  if (!exists("interval")) interval <- 1
  if (!exists("idxB_lag")) idxB_lag <- 1
  if (!exists("idxB_range_3")) idxB_range_3 <- 1
  if (!exists("upper_constraint")) upper_constraint <- Inf
  if (!exists("lower_constraint")) lower_constraint <- 0
  if (!exists("cap_below_b")) cap_below_b <- TRUE
  
  if (!exists("stat_yrs")) stat_yrs <- "all"
  if (!exists("scenario")) scenario <- "sensitivity"
  
  ### observation uncertainty
  if (!exists("sigmaL")) sigmaL <- 0.2
  if (!exists("sigmaB")) sigmaB <- 0.2
  ### recruitment variability
  if (!exists("sigmaR")) sigmaR <- 0.6
  if (!exists("sigmaR_rho")) sigmaR_rho <- 0.0
  
  ### what to save
  if (!exists("saveMP")) saveMP <- TRUE
  if (!exists("stats")) stats <- TRUE
  if (!exists("collate")) collate <- FALSE

}

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
path_in <- paste0("input/hr/", n_iter, "_", n_yrs, "/OM_2_mp_input/", fhist, "/")
input <- readRDS(paste0(path_in, stock, ".rds"))

input$args$nblocks <- n_blocks

### ------------------------------------------------------------------------ ###
### MP parameters ####
### ------------------------------------------------------------------------ ###

### load reference values
hr_ref <- readRDS("input/catch_rates.rds")[[stock]]
brp <- readRDS("input/brps.rds")[[stock]]
lhist <- stocks[stocks$stock == stock, ]

### HR rule parameters & uncertainty
hr_params <- data.frame(multiplier = multiplier,
                        comp_b = comp_b,
                        idxB_lag = idxB_lag,
                        idxB_range_3 = idxB_range_3,
                        interval = interval,
                        upper_constraint = upper_constraint,
                        lower_constraint = lower_constraint,
                        sigmaL = sigmaL,
                        sigmaB = sigmaB,
                        sigmaR = sigmaR,
                        sigmaR_rho = sigmaR_rho)

### ------------------------------------------------------------------------ ###
### go through runs ####
### ------------------------------------------------------------------------ ###

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
                  lower_constraint = par_i$lower_constraint,
                  cap_below_b = cap_below_b)
  
  ## ---------------------------------------------------------------------- ###
  ## observation uncertainty ####
  ## ---------------------------------------------------------------------- ###
  ### change uncertainty?
  sigmaB_i <- par_i$sigmaB
  sigmaL_i <- par_i$sigmaL
  if (par_i$sigmaB != 0.2 | par_i$sigmaL != 0.2) {
    
    input_i <- lapply(list(input), function(x) {
      
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
    input <- input_i[[1]]
    
  }
  
  ## ---------------------------------------------------------------------- ###
  ## recruitment variability ####
  ## ---------------------------------------------------------------------- ###
  ### change variability?
  sigmaR_i <- par_i$sigmaR
  sigmaR_rho_i <- par_i$sigmaR_rho
  if (sigmaR_i != 0.6 | sigmaR_rho_i != 0) {
    
    input_i <- lapply(list(input), function(x) {
      
      #browser()
      ### retrieve original residuals
      dev_R_original <- x$om@sr@residuals
      
      ### create recruitment residuals for projection period
      set.seed(1)
      dev_R_new <- rlnoise(dims(dev_R_original)$iter, dev_R_original %=% 0, 
                           sd = sigmaR_i, b = sigmaR_rho_i)
      ### replicate residuals from GA paper
      qnt_150 <- FLQuant(NA, 
                         dimnames = list(age = "all", year = 0:150,
                                         iter = dimnames(dev_R_original)$iter))
      qnt_100 <- FLQuant(NA, 
                         dimnames = list(age = "all", year = 1:100,
                                         iter = dimnames(dev_R_original)$iter))
      set.seed(0)
      res_150 <- rlnoise(dims(dev_R_original)$iter,
                         qnt_150 %=% 0, 
                         sd = sigmaR_i, b = sigmaR_rho_i)
      set.seed(0)
      res_100 <- rlnoise(dims(dev_R_original)$iter,
                         qnt_100 %=% 0, 
                         sd = sigmaR_i, b = sigmaR_rho_i)
      ### insert into template
      yrs_150 <- seq(from = dims(dev_R_original)$minyear,
                     to = ifelse(dims(dev_R_original)$maxyear >= 150, 
                                 150, dims(dev_R_original)$maxyear))
      dev_R_new[, ac(yrs_150)] <- res_150[, ac(yrs_150)]
      
      yrs_100 <- seq(from = dims(dev_R_original)$minyear,
                     to = 100)
      dev_R_new[, ac(yrs_100)] <- res_100[, ac(yrs_100)]
      
      
      ### insert
      x$om@sr@residuals[] <- dev_R_new
      
      return(x)
      
    })
    input <- input_i[[1]]
    
  }
  
  
  ### ---------------------------------------------------------------------- ###
  ### paths ####
  ### ---------------------------------------------------------------------- ###
  ### generate file name
  file_out <- paste0(c(hr, par_i$multiplier, par_i$comp_b, par_i$idxB_lag, 
                       par_i$idxB_range_3, par_i$interval, 
                       par_i$upper_constraint, par_i$lower_constraint,
                       par_i$sigmaL, par_i$sigmaB, par_i$sigmaR, 
                       par_i$sigmaR_rho), 
                     collapse = "_")
  path_out <- paste0("output/hr/", n_iter, "_", n_yrs, "/", scenario, "/",
                     fhist, "/", paste0(stock, collapse = "_"), "/")
  dir.create(path_out, recursive = TRUE)
  ### skip if run already exists
  if (file.exists(paste0(path_out, "stats_", file_out, ".rds"))) return(NULL)
  
  ### ---------------------------------------------------------------------- ###
  ### run  ####
  ### ---------------------------------------------------------------------- ###
  
  res <- do.call(mp, input)
  
  ### ---------------------------------------------------------------------- ###
  ### save ####
  ### ---------------------------------------------------------------------- ###
  
  if (isTRUE(saveMP))
    saveRDS(object = res, file = paste0(path_out, "mp_", file_out, ".rds"))
  
  ### ---------------------------------------------------------------------- ###
  ### stats ####
  ### ---------------------------------------------------------------------- ###
  
  if (isTRUE(stats)) {
    res_stats <- mp_stats(input = list(input), res_mp = list(res), 
                          collapse_correction = TRUE,
                          stat_yrs = stat_yrs)
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
  files <- paste0("output/hr/", n_iter, "_", n_yrs, "/", scenario, "/",
                  fhist, "/", paste0(stock, collapse = "_"), "/",
                  files)
  stats_all <- lapply(files, readRDS)
  stats_all <- do.call(rbind, stats_all)
  
  saveRDS(stats_all, file = paste0(
    "output/hr/", n_iter, "_", n_yrs, "/", scenario, "/", fhist, "/", 
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
