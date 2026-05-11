## tuning di lambda:
auto_lambda_ssj <- function(y, x = seq(1,length(y)), #lambda_min=1e-6, lambda_max=1e12, #3000
                            maxsum = sd(y, na.rm = TRUE)/mean(diff(x)),
                            edf = TRUE, last_delta = 1,
                            alpha = NULL, beta = 10) {
  
  bic_objective <- function(loglambda, y, x, maxsum, edf, last_delta) {
    ssj_fix(y = y, x = x, lambda = exp(loglambda), maxsum = maxsum,
            edf = edf, last_delta = last_delta)$ic['bic']
  }
  
  lambda_min <- smooth.spline(x=x,y=y)$lambda *diff(range(x))^3
  lambda_max <- lambda_min*beta
  
  opt_lambda <- optimize(
    bic_objective,
    interval = log(c(lambda_min, lambda_max)),
    y = y,
    x = x,
    maxsum = maxsum,
    edf = edf,
    last_delta = last_delta
  )
  
  lambda_fin <- exp(opt_lambda$minimum)
  ## somma pesata attraverso peso alfa:
  if (is.null(alpha)) {
    temp <- ssj_fix(y = y, x = x, lambda = lambda_fin, maxsum = maxsum, #opt_lambda$minimum
                    edf = edf, last_delta = last_delta)
    diffs <- abs(diff(temp$smoothed_level)) #proxy dei salti
    alpha <- max(diffs)/sum(diffs)
  }

  lambda_fin <- (1-alpha)*lambda_min +alpha*lambda_fin
  
  ssj_fix(y = y, x = x, lambda = lambda_fin, maxsum = maxsum, #opt_lambda$minimum
          edf = edf, last_delta = last_delta)
}

## tuning congiunto di M e lambda:
auto_ssj <- function(y, x = seq(1,length(y)), #lambda_min=1, lambda_max=1e4, #lambda_min=1e-6, lambda_max=1e12
                     grid = seq(0, sd(y, na.rm = TRUE)/mean(diff(x))*10, sd(y, na.rm = TRUE)/mean(diff(x))/10),
                     ic = c("bic", "hq", "aic", "aicc"),
                     edf = TRUE, last_delta = 1,
                     alpha = NULL, beta = 10) {
  ic <- match.arg(ic)
  last_ic <- Inf
  for (M in grid) {
    out <- auto_lambda_ssj(y = y, x = x,
                           #lambda_min = lambda_min, lambda_max = lambda_max,
                           maxsum = M, edf = edf, last_delta = last_delta,
                           alpha = alpha, beta = beta)
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


