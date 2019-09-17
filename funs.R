### ------------------------------------------------------------------------ ###
### observations ####
### ------------------------------------------------------------------------ ###
wklife_3.2.1_obs <- function(stk, observations, deviances, genArgs, tracking,
                             ssb = FALSE, ### use SSB as idx
                             idx_dev = FALSE,
                             lngth = FALSE, ### catch length data?
                             lngth_dev = FALSE, 
                             lngth_par,
                             ...) {

  #ay <- genArgs$ay
  ### update observations
  observations$stk <- stk
  ### use SSB as index?
  if (isTRUE(ssb)) {
    observations$idx$idxB <- ssb(observations$stk)
  ### otherwise calculate biomass index
  } else {
    observations$idx$idxB <- quantSums(stk@stock.n * stk@stock.wt * 
                                       observations$idx$sel)
  }
  ### use mean length in catch?
  if (isTRUE(lngth)) {
    observations$idx$idxL <- lmean(stk = stk, params = lngth_par)
  }
  
  ### observation model
  stk0 <- observations$stk
  idx0 <- observations$idx
  ### add deviances to index?
  if (isTRUE(idx_dev)) {
    if (isTRUE(ssb)) {
      idx0$idxB <- observations$idx$idxB * deviances$idx$idxB
    } else {
      idx0$idxB <- quantSums(stk@stock.n * stk@stock.wt * 
                             observations$idx$sel * deviances$idx$sel)
      if (isTRUE("idx" %in% names(deviances$idx)) & 
          all.equal(dim(deviances$idx$idxB), dim(idx0$idxB)))
        idx0$idxB <- idx0$idxB * deviances$idx$idxB
    }
  }
  ### uncertainty for catch length
  if (isTRUE(lngth) & isTRUE(lngth_dev)) {
    idx0$idxL <- observations$idx$idxL * deviances$idx$idxL
  }
  
  return(list(stk = stk0, idx = idx0, observations = observations,
              tracking = tracking))
  
}

### ------------------------------------------------------------------------ ###
### estimator ####
### ------------------------------------------------------------------------ ###

wklife_3.2.1_est <- function(stk, idx, tracking, genArgs,
                             comp_r = FALSE, comp_f = FALSE, comp_b = FALSE, 
                             idxB_lag = 1, idxB_range_1 = 2, idxB_range_2 = 3,
                             idxB_range_3 = 1,
                             catch_lag = 1, catch_range = 1,
                             Lref, I_trigger,
                             idxL_lag = 1, idxL_range = 1,
                             ...) {
  
  ay <- genArgs$ay
  
  ### component r: index trend
  if (isTRUE(comp_r)) {
    r_res <- est_r(idx = idx$idxB, ay = ay,
                   idxB_lag = idxB_lag, idxB_range_1 = idxB_range_1, 
                   idxB_range_2 = idxB_range_2)
  } else {
    r_res <- 1
  }
  tracking["comp_r", ac(ay)] <- r_res
  
  ### component f: length data
  if (isTRUE(comp_f)) {
    f_res <- est_f(idx = idx$idxL, ay = ay,
                   Lref = Lref, idxL_range = idxL_range, idxL_lag = idxL_lag)
  } else {
    f_res <- 1
  }
  tracking["comp_f", ac(ay)] <- f_res
  
  ### component b: biomass safeguard
  if (isTRUE(comp_b)) {
    b_res <- est_b(idx = idx$idxB, ay = ay,
                   I_trigger = I_trigger, idxB_lag = idxB_lag, 
                   idxB_range_3 = idxB_range_3)
  } else {
    b_res <- 1
  }
  tracking["comp_b", ac(ay)] <- b_res
  
  ### current catch
  catch_yrs <- seq(to = ay - catch_lag, 
                 length.out = catch_range)
  catch_current <- yearMeans(catch(stk)[, ac(catch_yrs)])
  tracking["C_current", ac(ay)] <- catch_current
  
  return(list(stk = stk, tracking = tracking))
  
}

### biomass index trend
est_r <- function(idx, ay,
                  idxB_lag, idxB_range_1, idxB_range_2,
                  ...) {
  
  ### index ratio
  yrs_a <- seq(to = c(ay - idxB_lag), length.out = idxB_range_1)
  yrs_b <- seq(to = min(yrs_a) - 1, length.out = idxB_range_2)
  idx_a <- yearMeans(idx[, ac(yrs_a)])
  idx_b <- yearMeans(idx[, ac(yrs_b)])
  idx_ratio <- c(idx_a / idx_b)
  
  return(idx_ratio)
  
}

### length data
est_f <- function(idx, ay, 
                  Lref, idxL_range, idxL_lag,
                  ...) {
  
  ### if fewer iterations provided expand
  if (isTRUE(length(Lref) < dims(idx)$iter)) {
    Lref <- rep(Lref, dims(idx)$iter)
    ### if more iterations provided, subset
  } else if (isTRUE(length(Lref) > dims(idx)$iter)) {
    Lref <- Lref[an(dimnames(idx)$iter)]
  }
  
  ### get mean length in catch
  idx_yrs <- seq(to = ay - idxL_range, length.out = idxL_lag)
  idx_mean <- yearMeans(idx[, ac(idx_yrs)])
  ### length relative to reference
  idx_ratio <- c(idx_mean / Lref)
  ### avoid negative values
  idx_ratio <- ifelse(idx_ratio > 0, idx_ratio, 0)
  return(idx_ratio)
}

### biomass index trend
est_b <- function(idx, ay, 
                  I_trigger, idxB_lag, idxB_range_3,
                  ...) {
  
  ### if fewer iterations provided expand
  if (isTRUE(length(I_trigger) < dims(idx)$iter)) {
    I_trigger <- rep(I_trigger, dims(idx)$iter)
  ### if more iterations provided, subset
  } else if (isTRUE(length(I_trigger) > dims(idx)$iter)) {
    I_trigger <- I_trigger[an(dimnames(idx)$iter)]
  }
  
  ### calculate index mean
  idx_yrs <- seq(to = ay - idxB_lag, length.out = idxB_range_3)
  idx_mean <- yearMeans(idx[, ac(idx_yrs)])
  ### ratio
  idx_ratio <- c(idx_mean / I_trigger)
  ### b is 1 or smaller
  idx_ratio <- ifelse(idx_ratio < 1, idx_ratio, 1)
  
  return(idx_ratio)
  
}

### ------------------------------------------------------------------------ ###
### phcr ####
### ------------------------------------------------------------------------ ###
### parametrization of HCR

phcr_r <- function(tracking, genArgs, multiplier = 1,
                   exp_r = 1, exp_f = 1, exp_b = 1,
                   ...){
  
  ay <- genArgs$ay
  
  hcrpars <- tracking[c("comp_r", "comp_f", "comp_b", "C_current", "multiplier",
                        "exp_r", "exp_f", "exp_b"), ac(ay)]
  hcrpars["multiplier", ] <- multiplier
  hcrpars["exp_r", ] <- exp_r
  hcrpars["exp_f", ] <- exp_f
  hcrpars["exp_b", ] <- exp_b
  
  if (multiplier != 1) tracking["multiplier", ] <- multiplier
  if (exp_r != 1) tracking["exp_r", ] <- exp_r
  if (exp_f != 1) tracking["exp_f", ] <- exp_f
  if (exp_b != 1) tracking["exp_b", ] <- exp_b
  
  
  ### return results
  return(list(tracking = tracking, hcrpars = hcrpars))
  
}


### ------------------------------------------------------------------------ ###
### hcr ####
### ------------------------------------------------------------------------ ###
### apply catch rule

hcr_r <- function(hcrpars, genArgs, tracking, interval = 1, 
                  ...) {
  
  ay <- genArgs$ay ### current year
  iy <- genArgs$iy ### first simulation year
  
  ### check if new advice requested
  if ((ay - iy) %% interval == 0) {
  
    ### calculate advice
    advice <- hcrpars["C_current", ] *
                (hcrpars["comp_r", ]^hcrpars["exp_r", ]) *
                (hcrpars["comp_f", ]^hcrpars["exp_f", ]) *
                (hcrpars["comp_b", ]^hcrpars["exp_b", ]) *
                hcrpars["multiplier", ] 
    #advice <- apply(X = hcrpars, MARGIN = 6, prod, na.rm = TRUE)
    
  } else {
    
    ### use last year's advice
    advice <- tracking["metric.hcr", ac(ay - 1)]
    
  }

  ctrl <- getCtrl(values = c(advice), quantity = "catch", years = ay + 1, 
                  it = dim(advice)[6])
  
  return(list(ctrl = ctrl, tracking = tracking))
  
}

### ------------------------------------------------------------------------ ###
### implementation ####
### ------------------------------------------------------------------------ ###
### no need to convert, already catch in tonnes
### apply TAC constraint, if required

is_r <- function(ctrl, genArgs, tracking, 
                  upper_constraint = Inf, lower_constraint = 0, ...) {
  
  ay <- genArgs$ay
  advice <- ctrl@trgtArray[ac(ay + genArgs$management_lag),"val",]
  
  ### apply TAC constraint, if requested
  if (!is.infinite(upper_constraint) | lower_constraint != 0) {
    
    ### get last advice
    adv_last <- tracking["metric.is", ac(ay - 1)]
    ### ratio of new advice/last advice
    adv_ratio <- advice/adv_last
    
    ### upper constraint
    if (!is.infinite(upper_constraint)) {
      ### find positions
      pos_upper <- which(adv_ratio > upper_constraint)
      ### limit advice
      if (length(pos_upper) > 0) {
        advice[,,,,, pos_upper] <- adv_last[,,,,, pos_upper] * upper_constraint
      }
      ### lower constraint
    }
    if (lower_constraint != 0) {
      ### find positions
      pos_lower <- which(adv_ratio < lower_constraint)
      ### limit advice
      if (length(pos_lower) > 0) {
        advice[,,,,, pos_lower] <- adv_last[,,,,, pos_lower] * lower_constraint
      }
    }
  }
  ctrl@trgtArray[ac(ay + genArgs$management_lag),"val",] <- advice
  
  return(list(ctrl = ctrl, tracking = tracking))
  
}

### ------------------------------------------------------------------------ ###
### implementation error ####
### ------------------------------------------------------------------------ ###

iem_r <- function(ctrl, genArgs, tracking, 
                  iem_dev = FALSE, use_dev, ...) {
  
  ay <- genArgs$ay
  
  ### only do something if requested
  if (isTRUE(use_dev)) {
    
    ### get advice
    advice <- ctrl@trgtArray[ac(ay + genArgs$management_lag), "val", ]
    ### get deviation
    dev <- c(iem_dev[, ac(ay)])
    ### implement deviation
    advice <- advice * dev
    ### insert into ctrl object
    ctrl@trgtArray[ac(ay + genArgs$management_lag),"val",] <- advice
    
  }
  
  return(list(ctrl = ctrl, tracking = tracking))
  
}

### ------------------------------------------------------------------------ ###
### projection ####
### ------------------------------------------------------------------------ ###
fwd_attr <- function(stk, ctrl,
                     sr, ### stock recruitment model
                     sr.residuals, ### recruitment residuals
                     sr.residuals.mult = TRUE, ### are res multiplicative?
                     maxF = 5, ### maximum allowed Fbar
                     dupl_trgt = FALSE,
                     ...) {
  
  ### avoid the issue that the catch is higher than the targeted catch
  ### can happen due to bug in FLash if >1 iteration provided
  ### sometimes, FLash struggles to get estimates and then uses F estimate from
  ### previous iteration
  ### workaround: target same value several times and force FLash to try again
  if (isTRUE(dupl_trgt)) {

    ### duplicate target
    ctrl@target <- rbind(ctrl@target, ctrl@target, ctrl@target)
    ### replace catch in second row with landings
    ctrl@target$quantity[1] <- "landings"
    ctrl@target$quantity[3] <- "catch"

    ### extract target values
    val_temp <- ctrl@trgtArray[, "val", ]

    ### extend trgtArray
    ### extract dim and dimnames
    dim_temp <- dim(ctrl@trgtArray)
    dimnames_temp <- dimnames(ctrl@trgtArray)
    ### duplicate years
    dim_temp[["year"]] <- dim_temp[["year"]] * 3
    dimnames_temp$year <- rep(dimnames_temp$year, 3)

    ### create new empty array
    trgtArray <- array(data = NA, dim = dim_temp, dimnames = dimnames_temp)

    ### fill with values
    ### first as target
    trgtArray[1, "val", ] <- val_temp
    ### then again, but as max
    trgtArray[2, "max", ] <- val_temp
    ### min F
    trgtArray[3, "max", ] <- val_temp

    ### insert into ctrl object
    ctrl@trgtArray <- trgtArray
  }
  
  ### project forward with FLash::fwd
  stk[] <- fwd(object = stk, control = ctrl, sr = sr, 
               sr.residuals = sr.residuals, 
               sr.residuals.mult = sr.residuals.mult,
               maxF = maxF)
  
  ### return stock
  return(list(object = stk))
  
}


### ------------------------------------------------------------------------ ###
### iter subset  ####
### ------------------------------------------------------------------------ ###

iter_attr <- function(object, iters, subset_attributes = TRUE) {
  
  ### subset object to iter
  res <- FLCore::iter(object, iters)
  
  if (isTRUE(subset_attributes)) {
    
    ### get default attributes of object class
    attr_def <- names(attributes(new(Class = class(object))))
    
    ### get additional attributes
    attr_new <- setdiff(names(attributes(object)), attr_def)
    
    ### subset attributes
    for (attr_i in attr_new) {
      attr(res, attr_i) <- FLCore::iter(attr(res, attr_i), iters)
    }
    
  }
  
  return(res)
  
}

### ------------------------------------------------------------------------ ###
### estimtate steepness based on l50/linf ratio ####
### according to Wiff et al. 2018
### ------------------------------------------------------------------------ ###
h_Wiff <- function(l50, linf) {
  l50linf <- l50/linf
  ### linear model
  lin <- 2.706 - 3.698*l50linf
  ### logit
  h <- (0.2 + exp(lin)) / (1 + exp(lin))
  return(h)
}

### ------------------------------------------------------------------------ ###
### mean length in catch ####
### ------------------------------------------------------------------------ ###
lmean <- function(stk, params) {
  
  ### calculate length from age with a & b
  weights <- c(catch.wt(stk)[, 1,,,, 1])
  lengths <- (weights / c(params["a"]))^(1 / c(params["b"]))
  catch.n <- catch.n(stk)
  dimnames(catch.n)$age <- lengths
  ### subset to lengths > Lc
  catch.n <- catch.n[lengths > c(params["Lc"]),]
  
  ### calculate mean length
  lmean <- apply(X = catch.n, MARGIN = c(2, 6), FUN = function(x) {
    ### calculate
    res <- weighted.mean(x = an(dimnames(x)$age), 
                         w = ifelse(is.na(x), 0, x), na.rm = TRUE)
    ### check if result obtained
    ### if all catch at all lengths = 0, return 0 as mean length
    # if (is.nan(res)) {
    #   if (all(ifelse(is.na(x), 0, x) == 0)) {
    #     res[] <- 0
    #   }
    # }
    return(res)
  })
  return(lmean)
}

### ------------------------------------------------------------------------ ###
### length at first capture ####
### ------------------------------------------------------------------------ ###
calc_lc <- function(stk, a, b) {
  ### find position in age vector
  Ac <- apply(catch.n(stk), MARGIN = c(2, 6), function(x) {
    head(which(x >= (max(x, na.rm = TRUE)/2)), 1)
  })
  Ac <- an(median(Ac))
  ### calculate lengths
  weights <- c(catch.wt(stk)[, 1,,,, 1])
  lengths <- (weights / a)^(1 / b)
  ### length at Ac
  Lc <- floor(lengths[Ac]*10)/10
  return(Lc)
}

