# FUNZIONI DI BASSO LIVELLO ####

## Stagionalità sinusoidale di frequenza 3 ####

seas <- function(x,range,min,freq=3) {
  sin((x-min)/range * 2*pi*freq) #-sin, cos
}


## Variabile indipendente ####

x_gen <- function(n,nj,bound,jumps_t, s=292) {
  set.seed(s)
  x <- sort(c(runif(n-2*nj, bound[1],bound[2]), jumps_t, jumps_t*.999))
}


## Simulazione trend LLT continuo ####

trend_LLT <- function(x,σ_ε,σ,bound,range,n,j,nj,jumps_ndx, s=18) {
  δ <- diff(x)
  
  set.seed(s)
  ζ <- rnorm(n, 0,σ)
  temp <- cumsum(ζ)
  for (i in 1:nj) {
    ζ[jumps_ndx[i]] <- ζ[jumps_ndx[i]] +j[i]*σ *ifelse(rbinom(1,1, .5), 1,-1)
  }
  β <- cumsum(ζ)
  μ <- c(0, cumsum(δ*β[-1]))
  
  μ_j <- numeric(n)
  for (i in 1:nj) {
    μ_j[jumps_ndx[i]:n] <- μ_j[jumps_ndx[i]:n] +range/20*j[i]*σ_ε
  }
  trend <- -3*x*σ
  stag <- range/1.5*seas(x,range,bound[1])*σ

  list(μ +trend +stag +μ_j, #trend con salti
       c(0,cumsum(δ*temp[-1])) +trend +stag) #trend senza salti
}


## Simulazione di N serie LLT continue ####

series_LLT <- function(N,tr, range,σ_ε,n, s=18) {
  Y <- list()
  Y_detr <- list()
  set.seed(s)
  for (i in 1:N) {
    ε <- range/20*rnorm(n, 0,σ_ε)
    Y[[i]] <- tr[[1]] +ε
    Y_detr[[i]] <- tr[[2]] +ε
  }
  list(Y, Y_detr)
}


## Applicazione iterativa del metodo sulle serie simulate ####

sim_LLT <- function(N, x,Y,Y_detr, range,n,nj,jumps_ndx,
                    func=c('ssj','alt','exp'),
                    method=c('fix','mle','auto')) {
  func <- match.arg(func)
  method <- match.arg(method)
  
  sample_laSS <- numeric(1)
  sample_la <- numeric(1)
  sample_hla <- numeric(1)
  sample_tr <- list()
  sample_j <- list()
  sample_t <- numeric(1)
  
  f <- switch (func,
    'ssj' = switch (method,
      'fix'  = function(λ) auto_ssj_fix(x=x,y=y,lambda=λ), #λ*(range/n)^3 #smooth.spline(x=x,y=y)$lambda*range^3
      'mle'  = function(λ) auto_ssj_mle(x=x,y=y),
      'auto' = function(λ) auto_ssj(x=x,y=y)
    ),
    'alt' = switch (method,
      'fix'  = function(λ) auto_ssjAlt_fix(x=x,y=y,lambda=λ),
      'mle'  = function(λ) auto_ssjAlt_mle(x=x,y=y),
      'auto' = function(λ) auto_ssjAlt(x=x,y=y)
    ),
    'exp' = switch (method,
      'fix'  = function(λ) auto_ssjExp_fix(x=x,y=y,lambda=λ),
      'mle'  = function(λ) auto_ssjExp_mle(x=x,y=y),
      'auto' = function(λ) auto_ssjExp(x=x,y=y)
    )
  )
  
  g <- switch (func,
    'exp' = function() c(0, as.numeric(filter$gammas > 0)),
    function() c(0, as.numeric(filter$sigmas > 0)[-1])
  )
  
  for (i in 1:N) {print(paste0(i,"/",N))
    y <- Y[[i]]
    lambda <- smooth.spline(x=x,y=Y_detr[[i]])$lambda*range^3
    #print(c(smooth.spline(x=x,y=y)$lambda*range^3, lambda))
    
    start <- Sys.time()
    #
    filter <- f(lambda) #λ=40
    #
    end <- Sys.time()
    
    #print((filter$pars[['sigma_noise']]/filter$pars[['sigma_slope']])^2)
    #filter <- auto_ssj_fix(x=x,y=y,lambda=1.299201e-05)λ/n^3
    sample_laSS[[i]] <- smooth.spline(x=x,y=y)$lambda*range^3
    sample_la[[i]] <- lambda
    sample_hla[[i]] <- filter$pars[['lambda']]
    sample_tr[[i]] <- filter$smoothed_level
    sample_j[[i]] <- g() #c(0,as.numeric(filter$sigmas > 0)[-1])
    sample_t[[i]] <- as.numeric(difftime(end, start, units = "secs"))
  }
  #smooth_dist <- Reduce("+", sample_tr)/N #previsione media
  mat_tr <- do.call(rbind, sample_tr)
  lower <- apply(mat_tr, 2, quantile, probs = .025)
  med <- apply(mat_tr, 2, mean) #previsione media
  #med <- apply(mat_tr, 2, quantile, probs = .5) #previsione mediana
  upper <- apply(mat_tr, 2, quantile, probs = .975)
  
  jump_dist <- Reduce("+", sample_j)/N
  
  TPR <- sum(jump_dist[jumps_ndx])/nj
  TNR <- sum(-(jump_dist-1)[-jumps_ndx])/(n-nj)
  BA <- mean(c(TPR,TNR)) # balanced accuracy
  
  list(
    smooth = med,
    lower = lower,
    upper = upper,
    jump_dist = jump_dist,
    metr = c(
      lambdaSS = mean(sample_laSS),
      lambda = mean(sample_la),
      hlambda = mean(sample_hla),
      BA  = BA,
      TPR = TPR,
      TNR = TNR,
      time = mean(sample_t)
    )
  )
}


## Metodo grafico ####

sim_graph <- function(x,μ, sim_ris,jumps_t, bound,range,n, BW=T) {
  metrics <- round(sim_ris$metr, 2)
  # tit <- paste0("BA: ", metrics['BA'],
  #              ",  TPR: ", metrics['TPR'],
  #              ",  TNR: ", metrics['TNR'],
  #              ",  time: ", metrics['time'],
  #              ",  ", expression(lambda),": ", metrics['lambda'],
  #              ",  ", expression(hat(lambda)),": ", metrics['hlambda'])
  
  tit <- bquote(lambda[SS] ~ ":" ~ .(metrics['lambdaSS']) ~ 
              ",  " ~ lambda ~ ":" ~ .(metrics['lambda']) ~ 
                ",  " ~ hat(lambda) ~ ":" ~ .(metrics['hlambda']) ~ 
                ",  BA:" ~ .(metrics['BA']) ~ 
                ",  TPR:" ~ .(metrics['TPR']) ~ 
                ",  TNR:" ~ .(metrics['TNR']) ~ 
                ",  time:" ~ .(metrics['time']))
  
  #par(mfrow=c(2,1), mar = c(2, 3, 2, 2))#, cex.main=3, cex.axis=2
  #par(cex.main=2)
  
  if (BW) {
    clrs = c('black','black', rgb(0, 0, 0, 0.2), 'darkgray')#darkgray
    ltps = c(2, 1)
  } else {
    clrs = c('red','blue', rgb(0, 0, 1, 0.2), 'cornflowerblue')
    ltps = c(1, 1)
  }
  
  mrgn <- .025*diff(range(c(μ, sim_ris$lower, sim_ris$upper)))
  lims <- c(min(c(μ, sim_ris$lower, sim_ris$upper)) -mrgn,
            max(c(μ, sim_ris$lower, sim_ris$upper)) +mrgn)
  
  layout(mat=matrix(c(1,2), nrow=2,ncol=1), heights = c(3,2))
  #par(mar = c(1.3, 3, 2, 2))
  par(mar = c(1, 1, 2, 1))
  
  plot(x=x, y=μ, type='n', main=tit, cex.main=.85, xaxt = "n",yaxt = "n", ylim=lims)
  lines(x=x, y=μ, lty=ltps[1], col=clrs[1], lwd=2)
  
  polygon(
    c(x, rev(x)),
    c(sim_ris$lower, rev(sim_ris$upper)),
    col = clrs[3],
    border = NA
  )
  lines(x=x, y=sim_ris$smooth, col=clrs[2], lwd=2.5, lty=as.numeric(ltps[2]))
  lines(x=x, y=sim_ris$lower, col=clrs[2], lwd=1.5)
  lines(x=x, y=sim_ris$upper, col=clrs[2], lwd=1.5)
  
  abline(v=jumps_t, lty=3, lwd=2)
  
  par(mar = c(2, 1, 0, 1))
  
  plot(
    x=x, y=sim_ris$jump_dist,
    type = "h",
    lwd = 4,
    col = clrs[4],#steelblue,deepskyblue4,
    xaxt = "n",
    yaxt = "n",
    xlab = "x",
    #ylab = "Probabilità",
    ylim = c(0, 1)
  )
  #axis(1, at = seq(0, 1, by = .1), labels = round(seq(0, 1, by = .1), 2))
  #axis(1, at = seq(bound[1], bound[2], by = range/10), labels = seq(bound[1], bound[2], by = range/10))
  axis(1, at=c(bound[1], jumps_t, bound[2]), labels=c(bound[1], jumps_t, bound[2]))
  #axis(2, at = seq(0, 1, by = 0.5), labels = seq(0, 1, by = 0.5))
  
  abline(v=jumps_t, lty=3, lwd=2)
  abline(h=c(0,1), col='white', lwd=4)
  abline(h=c(0,1), col='black', lwd=1.5)
  
  box()
  
  #par(mfrow=c(1,1), cex.main=1, cex.axis=1)
  layout(mat=matrix(1, nrow=1,ncol=1), heights = 1)
  #par(cex.main=1)
}



# FUNZIONI AD ALTO LIVELLO ####

## Esempio minimo di generazione dati ####

esemp <- function(
    n = 100, #numero di osservazioni
    bound = c(0,n), #intervallo di campionamento
    range = diff(bound),
    
    σ = 1, #scarto quadratico medio di ζ
    λ = 40, #σ_ε2/σ2
    
    j = c(15, -5, -10), #salti
    jt = c(.3, .4, .6), #posizioni relative dei salti
    nj = length(j),

    save = NULL #nome del grafico da salvare
    ) {
  σ_ε <- sqrt(λ)*σ
  
  jumps_t <- bound[1] +jt*range
  
  x <- x_gen(n,nj,bound,jumps_t)
  
  jumps_ndx <- which(x %in% jumps_t)
  
  tr <- trend_LLT(x=x,σ_ε=σ_ε,σ=σ,
                 bound=bound,range=range,n=n,
                 j=j,nj=nj,jumps_ndx=jumps_ndx)
  
  Y <- series_LLT(N=1,tr=tr,range=range,σ_ε=σ_ε,n=n)
  y <- Y[[1]][[1]]
  
  smooth <- smooth.spline(x,y)
  start <- Sys.time()
  smooth_nojumps <- smooth.spline(x,Y[[2]][[1]])
  end <- Sys.time()
  
  print(paste0("λ (SS): ", round(smooth$lambda*diff(range(x))^3,2),
               ", λ: ", round(smooth_nojumps$lambda*diff(range(x))^3,2),
               ", SS time (s): ", as.numeric(difftime(end, start, units = "secs"))))
  
  if (!is.null(save)) pdf(save,width=6, height=3)
  par(mar = c(4, 4, 1, 1))
  plot(x=x,y=y, type='p',pch=1,cex=.7, xlab='x',ylab='y')
  lines(x=x,y=tr[[1]], col='black', lwd=2)
  lines(x=smooth$x,y=smooth$y, col='black', lty=6,lwd=2)
  #lines(x=smooth$x,y=smooth$y, col='gray', lty=1,lwd=1.5)
  abline(v=jumps_t, lty=3, lwd=2)
  if (!is.null(save)) dev.off()
}

## Bayesian changepoint detection ####

esemp_bayes <- function(
    n = 100, #numero di osservazioni
    bound = c(0,n), #intervallo di campionamento
    range = diff(bound),
    
    σ = 1, #scarto quadratico medio di ζ
    λ = 40, #σ_ε/σ
    
    j = c(15, -5, -10), #salti
    jt = c(.3, .4, .6), #posizioni relative dei salti
    nj = length(j),
    
    R = 8e2,#1.5e3
    burn = floor(R*.4),
    
    s = 1,
    save = NULL #nome del grafico da salvare
) {
  σ_ε <- sqrt(λ)*σ
  
  jumps_t <- bound[1] +jt*range
  
  x <- x_gen(n,nj,bound,jumps_t)
  
  jumps_ndx <- which(x %in% jumps_t)
  
  tr <- trend_LLT(x=x,σ_ε=σ_ε,σ=σ,
                 bound=bound,range=range,n=n,
                 j=j,nj=nj,jumps_ndx=jumps_ndx)
  
  y <- series_LLT(N=1,tr=tr,range=range,σ_ε=σ_ε,n=n)[[1]][[1]]
  
  set.seed(s)
  sim <- changepoint_detection(y, R=R)
  
  if (!is.null(save)) pdf(save,width=6, height=3)
  conv(sim,x, burn=burn,cut=T)
  if (!is.null(save)) dev.off()
  
  ris <- get_opt_part(sim,burn)
  print(paste0("φ estim: ", round(ris$phi,2)))
}

bayes <- function(
    N = 1e3, #numero di iterazioni
    
    n = 100, #numero di osservazioni
    bound = c(0,n), #intervallo di campionamento
    range = diff(bound),
    
    σ = 1, #scarto quadratico medio di ζ
    λ = 40, #σ_ε/σ
    
    j = c(15, -5, -10), #salti
    jt = c(.3, .4, .6), #posizioni relative dei salti
    nj = length(j),
    
    save = NULL, #nome del grafico da salvare
    
    s = 1, #seed
    
    φ = NULL,
    R = 400,
    burn = floor(R*.4)
) {
  σ_ε <- sqrt(λ)*σ
  
  jumps_t <- bound[1] +jt*range
  
  x <- x_gen(n,nj,bound,jumps_t)
  
  jumps_ndx <- which(x %in% jumps_t)
  
  tr <- trend_LLT(x=x,σ_ε=σ_ε,σ=σ,
                 bound=bound,range=range,n=n,
                 j=j,nj=nj,jumps_ndx=jumps_ndx)
  μ <- tr[[1]]
  
  Y <- series_LLT(N=N,tr=tr,range=range,σ_ε=σ_ε,n=n)[[1]]
  
  sample_j <- list()
  sample_t <- numeric(1)
  
  set.seed(s)
  for (i in 1:N) {print(paste0(i,"/",N))
    y <- Y[[i]]
    
    start <- Sys.time()
    #
    sim <- changepoint_detection(y, R=R,iter=F, φ=φ)#2000, x
    conv(sim,x, burn=burn,cut=T)
    ris <- get_opt_part(sim,burn, iter=F)
    #
    end <- Sys.time()
    
    graph(ris,y,x=x)
    
    print(paste0("prob change φ: ", round(sim$pchange,2)))
    print(paste0("φ estim: ", round(ris$phi,2)))
    
    jumps <- numeric(n)
    jumps[ris$jumps] <- 1
    
    sample_j[[i]] <- jumps
    sample_t[[i]] <- as.numeric(difftime(end, start, units = "secs"))
  }
  jump_dist <- Reduce("+", sample_j)/N
  
  TPR <- round(sum(jump_dist[jumps_ndx])/nj,2)
  TNR <- round(sum(-(jump_dist-1)[-jumps_ndx])/(n-nj),2)
  BA <- round(mean(c(TPR,TNR)),2)
  time = round(mean(sample_t),2)
  
  tit <- bquote("BA: " ~ .(BA) ~ 
                ",    TPR: " ~ .(TPR) ~ 
                ",    TNR: " ~ .(TNR) ~ 
                ",    time: " ~ .(time))
  
  # tit <- paste0("BA: ", round(BA,2),
  #               ",    TPR: ", round(TPR,2),
  #               ",    TNR: ", round(TNR,2),
  #               ",    time: ", round(time,2))
  
  if (!is.null(save)) pdf(save,width=6, height=3)
  
  layout(mat=matrix(c(1,2), nrow=2,ncol=1), heights = c(3,2))
  par(mar = c(1, 2, 2, 1))
  
  plot(x=x, y=μ, type='n', main=tit, cex.main=.85, xaxt = "n")
  lines(x=x, y=μ, lty=1, col='black', lwd=2)
  
  abline(v=jumps_t, lty=3, lwd=2)
  
  par(mar = c(2, 2, 0, 1))
  
  plot(
    x=x, y=jump_dist,
    type = "h",
    lwd = 4,
    col = 'darkgray',#steelblue,deepskyblue4,
    xaxt = "n",
    yaxt = "n",
    xlab = "x",
    ylab = "Probabilità",
    ylim = c(0, 1)
  )
  #axis(1, at = seq(0, 1, by = .1), labels = round(seq(0, 1, by = .1), 2))
  #axis(1, at = seq(bound[1], bound[2], by = range/10), labels = seq(bound[1], bound[2], by = range/10))
  axis(1, at=c(bound[1], jumps_t, bound[2]), labels=c(bound[1], jumps_t, bound[2]))
  axis(2, at = seq(0, 1, by = 0.5), labels = seq(0, 1, by = 0.5))
  
  abline(v=jumps_t, lty=3, lwd=2)
  abline(h=c(0,1), col='white', lwd=4)
  abline(h=c(0,1), col='black', lwd=1.5)
  
  box()
  
  layout(mat=matrix(1, nrow=1,ncol=1), heights = 1)
  
  if (!is.null(save)) dev.off()
}

## Simulazione ####

simulation <- function(
    func=c('ssj','alt','exp'), #funzione da testare
    method=c('fix','mle','auto'), #lambda
    N = 1e3, #numero di iterazioni
    
    n = 100, #numero di osservazioni
    bound = c(0,n), #intervallo di campionamento
    range = diff(bound),
    
    σ = 1, #scarto quadratico medio di ζ
    λ = 40, #σ_ε/σ
    
    j = c(15, -5, -10), #salti
    jt = c(.3, .4, .6), #posizioni relative dei salti
    nj = length(j),
    
    #s = 1, #seed
    
    save = NULL #nome del grafico da salvare
    ) {
  σ_ε <- sqrt(λ)*σ
  jumps_t <- bound[1] +jt*range
  
  x <- x_gen(n,nj,bound,jumps_t)
  jumps_ndx <- which(x %in% jumps_t)
  
  tr <- trend_LLT(x=x,σ_ε=σ_ε,σ=σ,
                  bound=bound,range=range,n=n,
                  j=j,nj=nj,jumps_ndx=jumps_ndx)
  
  series <- series_LLT(N=N,tr=tr, range=range,σ_ε=σ_ε,n=n)
  
  sim <- sim_LLT(N=N,func=func,method=method,
                 x=x,Y=series[[1]],Y_detr=series[[2]],
                 range=range,n=n,nj=nj,jumps_ndx=jumps_ndx)
  
  if (!is.null(save)) pdf(save,width=6, height=3)
  sim_graph(x=x,μ=tr[[1]], sim_ris=sim,jumps_t=jumps_t,
            bound=bound,range=range,n=n)
  if (!is.null(save)) dev.off()
}

