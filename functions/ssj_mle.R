ssj_mle <- function(y, x = seq(1,length(y)), maxsum = sd(y, na.rm = TRUE)/mean(diff(x)),
                    #y, x, maxsum = sd(y),
                    edf = TRUE, parinit = NULL, last_delta = 1) {
  n <- length(y)
  n1 <- n+1
  nobs <- sum(!is.na(y))
  vy <- var(y, na.rm = TRUE)
  vdy <- var(diff(y), na.rm = TRUE)
  sdy <- sqrt(vy)
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  delta <- c(diff(x), last_delta)
  
  var_eps <- numeric(n)
  var_eta <- numeric(n)
  var_zeta <- numeric(n)
  cov_eta_zeta <- numeric(n)
  a1  <- rep(y[!is.na(y)][1], n1)
  a2  <- numeric(n1)
  p11 <- rep(vy*1.0e5, n1)
  p12 <- numeric(n1)
  p22 <- rep(vdy*1.0e5, n1)
  k1  <- numeric(n)
  k2  <- numeric(n)
  i   <- numeric(n)
  f   <- numeric(n)
  r1  <- numeric(n1)
  r2  <- numeric(n1)
  n11 <- numeric(n1)
  n12 <- numeric(n1)
  n22 <- numeric(n1)
  e   <- numeric(n1)
  d   <- numeric(n1)
  w   <- numeric(n)
  cnst <- log(2*pi)/2
  
  D <- array(0, c(2, 2, n))
  D[1, 1, ] <- delta^3/3
  D[1, 2, ] <- D[2, 1, ] <- delta^2/2
  D[2, 2, ] <- delta
  
  # pars[1] sigma
  # pars[2] sigma_eps
  # pars[3:(n+2)] sigma_t
  vt_eta_ndx <- 3:(n+2)
  obj <- function(pars, wgt = FALSE) {
    sig2  <- pars[1]*pars[1] # common variance
    sig2t <- pars[vt_eta_ndx]*pars[vt_eta_ndx] # additional variances
    sum_sig2t <- sig2 + sig2t # common variance + additional variances
    var_eps[]  <- pars[2]*pars[2] # variance of the noise
    var_eta[]  <- D[1, 1, ]*sum_sig2t # variance of eta
    var_zeta[] <- D[2, 2, ]*sum_sig2t # variance of zeta
    cov_eta_zeta[] <- D[1, 2, ]*sum_sig2t # covariance of eta and zeta
    
    if (wgt) { # compute weights for edf
      mloglik <- -llt_delta(y, delta, var_eps, var_eta, var_zeta, cov_eta_zeta,
                            a1, a2, p11, p12, p22,
                            k1, k2, i, f, r1, r2,
                            n11, n12, n22, e, d, w)
    } else { # no weights computed
      mloglik <- -llt_delta(y, delta, var_eps, var_eta, var_zeta, cov_eta_zeta,
                            a1, a2, p11, p12, p22,
                            k1, k2, i, f, r1, r2,
                            n11, n12, n22, e, d)
    }
    rn1  <- (r1[-1]*r1[-1] - n11[-1]) # per calcolare il gradiente
    rn2  <- (r2[-1]*r2[-1] - n22[-1]) # =
    rn12 <- (r1[-1]*r2[-1] - n12[-1]) # =
    rn_sum <- rn1*D[1, 1, ] + 2*rn12*D[1, 2, ] + rn2*D[2, 2, ]    # =
    list(
      # Minus mean log-likelihood for computational stability
      objective = mloglik/nobs,
      # Average gradient
      gradient  = c(
        # w.r.t sigma_zeta
        -pars[1]*sum(rn_sum),
        # w.r.t sigma_epsilon
        -pars[2]*sum(e*e-d),
        # w.r.t. sigma_eta_t
        -pars[vt_eta_ndx]*rn_sum
      )/nobs
    )
    
  }
  
  ##### Constraints to the object function
  g <- function(pars, wgt) {
    list(
      constraints = sum(pars[vt_eta_ndx]) - maxsum, # constraint <= 0
      jacobian = c(0, 0, rep(1, n)) # derivatives of constr. w.t. to σ_t
    )
  }
  
  ##### Optimization step
  ## Starting values
  if (is.null(parinit)) {
    inits <- c(sd_zeta = sdy/10, sd_eps = sdy, rep(1, n-1), 0) #
  } else {
    inits <- parinit
  }
  ## Check on starting values
  # lb <- c(0, 0, rep(0, n))
  quasizero <- sdy*1.0e-9
  lb <- c(quasizero, quasizero, rep(0, n)) # lower bound
  inits[inits < lb] <- lb[inits < lb]      # fix to lower bound those that are zerO
  
  ## Optimization with CCSA ("conservative convex separable approximation")
  # see https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/
  opt <- nloptr::nloptr(x0 = inits,
                        eval_f = obj,
                        lb = lb,
                        eval_g_ineq = g, # constraint
                        opts = list(algorithm = "NLOPT_LD_CCSAQ",
                                    xtol_rel = 1.0e-5,
                                    check_derivatives = FALSE,
                                    maxeval = 2000),
                        wgt = FALSE #... -> passa ad obj e g
  )
  if (edf == TRUE) {
    obj(opt$solution, TRUE) # update pointed variables using estimated values
    df <- sum(w)
  } else {
    df <- 2 + sum(opt$solution[vt_eta_ndx] > quasizero) # come LASSO per lin reg
  }
  loglik <- -nobs*(opt$objective + cnst)
  
  ##### Output list
  list(opt = opt,
       nobs = n,
       df = df,
       maxsum = maxsum,
       loglik = loglik,
       pars = c(sigma_slope = opt$solution[1],
                sigma_noise = opt$solution[2],
                lambda = opt$solution[2]*opt$solution[2]/opt$solution[1]/opt$solution[1]),
       sigmas = c(NA, opt$solution[-c(1, 2, n+2)]), # il primo è vuoto perché c'è un lag negli shock
       weights = if (edf) w else NULL,
       ic = c(aic  = 2*(df - loglik),
              aicc = 2*(df*n/(n-df-1) - loglik),
              bic  = df*log(n) - 2*loglik,
              hq   = 2*(log(log(n))*df - loglik)),
       smoothed_level = (a1 + p11*r1 + p12*r2)[-(n+1)],
       var_smoothed_level = (p11 - p11*p11*n11 - 2*p11*p12*n12 - p12*p12*n22)[-(n+1)],
       x = x
  )
}

auto_ssj_mle <- function(y, x = seq(1,length(y)),
                         grid = seq(0, sd(y, na.rm = TRUE)/mean(diff(x))*10, sd(y, na.rm = TRUE)/mean(diff(x))/10),
                         #grid = seq(0, sd(y)*10, sd(y)/10), #101 punti
                         ic = c("bic", "hq", "aic", "aicc"),
                         edf = TRUE, last_delta = 1) {
  ic <- match.arg(ic)
  k <- length(grid)
  last_ic <- Inf
  for (M in grid) {
    out <- ssj_mle(y = y, x = x, maxsum = M,
                   edf = edf, last_delta = last_delta)
    current_ic <- switch (ic,
                          bic = out$ic["bic"],
                          hq = out$ic["hq"],
                          aic = out$ic["aic"],
                          aicc = out$ic["aicc"]
    )
    if (current_ic < last_ic) {
      best <- out
      last_ic <- current_ic
    }
  }
  best
}
