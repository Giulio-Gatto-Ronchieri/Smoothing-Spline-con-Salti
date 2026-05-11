# FUNZIONI AUSILIARIE ####

## Logit ####
logit <- function(x) log(x/(1-x))
logit_inv <- function(x) exp(x)/(1+exp(x))

## log-verosimiglianza ####

loglik_seg <- function(yj, δj, φ=.3, c=0.001, a=2, b=5){
  nj <- length(yj)
  
  # costruzione della matrice Sj
  d <- rep(1,nj)
  phi <- φ^δj
  if (nj > 2) {
    d[2:(nj-1)] <- 1 +(phi[1:(nj-2)]^2 +phi[2:(nj-1)]^2) / 2
  }
  Sj <- diag(d)
  if (nj>1) {
    for (i in 1:(nj-1)) {
      Sj[i,i+1] <- Sj[i+1,i] <- -phi[i]
    }
  }
  
  # (log-)verosimiglianza del cluster, formula 9 pagina 11
  if (nj==1) {
    somma <- 0
  } else {
    somma <- sum(yj)-yj[1]-yj[nj]
  }
  
  φj <- φ^mean(δj)
  unomph <- 1-φj
  unopph <- 1+φj
  
  Q <- as.numeric( sum(yj * (Sj %*% yj)) #t(yj)%*%Sj%*%yj
                   -((unomph)*(sum(yj)-φj*somma)*(sum(yj)-φj*somma)) / (c+nj-φj*(nj-c-2))
                   +2*b*unopph*unomph )
  if (Q<=0) return(NA)
  
  lver = as.numeric( a*(log(2)+log(b)+log(unopph)+log(unomph)) +lgamma(nj/2+a) -(nj/2)*log(pi) -lgamma(a)
              +(1/2) * ( log(c)+2*log(unopph)+log(unomph) -log(c+nj-φj*(nj-c-2)) )
              -(nj/2+a) * log(Q) )
  #print(paste("lver =",lver))
  lver
}

loglik <- function(y, δ, ρ, φ){
  start <- 1
  loglik <- 0
  for (nj in ρ) {
    yj <- y[start:(start+nj-1)]
    if (nj==1) {
      δj <- 1
    } else {
      δj <- δ[start:(start+nj-2)]
    }
    loglik <- loglik + loglik_seg(yj=yj,δj=δj,φ=φ)
    start <- start + nj
  }
  loglik#; print(loglik)
}

## Probabilità di accettazione ####

α_split <- function(y,δ,ρ,q,k,ρ_new,nj,φ,monit) {
  ODDS <- 1-q
  if (k>1) {
    ODDS <- ODDS/q
  }
  Ratio <- sum(ρ>1)*(nj-1)/k #se k=1, sum(ρ>1) vale 1 e nj vale n
  loglikRatio <- loglik(y,δ,ρ_new,φ=φ) - loglik(y,δ,ρ,φ=φ)
  
  if (is.na(loglikRatio)) return(0)
  α <- exp(log(ODDS) +log(Ratio) + loglikRatio)
  if (monit) print(paste("α_split =", α))
  return(α)
}

α_merge <- function(y,δ,n,ρ,q,k,ρ_new,nj,φ,monit) {
  ODDS <- q
  if (k<n) {
    ODDS <- ODDS/(1-q)
  }
  Ratio <- (k-1) / ( max(sum(ρ>1),1)*(nj-1) ) #se k=n, nj vale 2 e sum(ρ>1) vale 0
  loglikRatio <- loglik(y,δ,ρ_new,φ=φ) -loglik(y,δ,ρ,φ=φ)
  
  if (is.na(loglikRatio)) return(0)
  α <- exp(log(ODDS) +log(Ratio) +loglikRatio)
  if (monit) print(paste("α_merge =", α))
  return(α)
}

α_shuffle <- function(y,δ,ρ,ρ_new,φ,monit) {
  loglikRatio <- loglik(y,δ,ρ_new,φ=φ) -loglik(y,δ,ρ,φ=φ)
  
  if (is.na(loglikRatio)) return(0)
  α <- exp(loglikRatio)
  if (monit) print(paste("α_shuffle =", α))
  return(α)
}

α_φ <- function(y,δ,ρ,φ,φ_new,monit) {
  loglikRatio <- loglik(y,δ,ρ,φ=φ_new) -loglik(y,δ,ρ,φ=φ)
  
  if (is.na(loglikRatio)) return(0)
  α <-exp(loglikRatio
      +log(φ_new*(1-φ_new)) -log(φ*(1-φ)))
  if (monit) print(paste("α_φ =", α))
  return(α)
}

## Entropia ####
entr <- function(x,k=length(x)) {
  if (k==1){
    0
  } else {
    p <- x/sum(x)
    -sum(p*log(p)) / log(k)
  }
}



# FUNZIONI AD ALTO LIVELLO ####

## Campionamento MC ####

changepoint_detection <- function(y, x=NULL, φ=NULL, R=500, q=.5, σ=1, # se φ viene dato in input, rimane fisso per tutte le iterazioni
                                  iter=TRUE, monit=FALSE) {
  if (monit) iter <- TRUE
  
  n <- length(y)
  ρ <- n
  k <- length(ρ)
  if (is.null(x)) {
    δ <- rep(1,n-1)
  } else {
    δ <- diff(x)/mean(diff(x))
  }
  
  # campioni
  sample_ρ <- list()
  sample_j <- list()
  sample_k <- numeric(1)
  sample_h <- numeric(1)
  
  # φ automatico
  fix <- TRUE
  if (is.null(φ)) {
    φ <- .5
    fix <- FALSE
  }
  cont <- 0
  sample_φ <- φ
  
  for (r in 1:R) {if (iter) print(paste0(r,"/",R))
    if (runif(1) <= q*(k<n)+(k==1)) { # SPLIT
      valid <- which(ρ>1)
      if (length(valid)>1) j <- sample(which(ρ>1), 1)
      else j <- valid
      nj <- ρ[j]
      l <- sample(1:(nj-1), 1)
      
      ρ_new <- ρ
      ρ_new[j] <- l
      ρ_new <- append(ρ_new, nj-l, j)
      #temp=α_split(y,ρ,q,k,ρ_new,nj,φ=φ)
      if (runif(1) <= α_split(y,δ,ρ,q,k,ρ_new,nj,φ=φ,monit=monit)) {
        ρ <- ρ_new
      }
    } else { # MERGE
      j <- sample(1:(k-1), 1)
      nj <- ρ[j] +ρ[j+1]
      
      ρ_new <- ρ
      ρ_new[j] <- nj
      ρ_new <- ρ_new[-(j+1)]
      #temp=α_merge(y,n,ρ,q,k,ρ_new,nj,φ=φ)
      if (runif(1) <= α_merge(y,δ,n,ρ,q,k,ρ_new,nj,φ=φ,monit=monit)) {
        ρ <- ρ_new
      }
    }
    k <- length(ρ)
    sample_k[r] <- k
    
    if (k>1) { # SHUFFLE
      i <- sample(1:(k-1), 1)
      j <- sample(1:(ρ[i]+ρ[i+1]-1), 1)
      
      ρ_new <- ρ
      ρ_new[i+1] <- ρ[i] +ρ[i+1] -j
      ρ_new[i] <- j
      #temp=α_shuffle(y,ρ,ρ_new,φ=φ)
      if (runif(1) <= α_shuffle(y,δ,ρ,ρ_new,φ=φ,monit=monit)) {
        ρ <- ρ_new
      }
    }
    sample_ρ[[r]] <- ρ
    sample_j[[r]] <- head(cumsum(ρ)+1, -1)
    sample_h[[r]] <- entr(ρ,k)
    
    # MH per φ
    if (!fix) {
      φ_new <- logit_inv(rnorm(1, logit(φ),σ))
      #temp=α_φ(y,ρ,φ=φ,φ_new=φ_new)
      if (runif(1) <= α_φ(y,δ,ρ,φ=φ,φ_new=φ_new,monit=monit)) {
        φ <- φ_new
      }
    }
    sample_φ[r+1] <- φ
    if (φ!=sample_φ[r]) cont <- cont+1
    if (cont/r < .25) σ <- σ*.99
    else if (cont/r > .5) σ <- σ*1.01
    #print(σ)
  }
  
  # output
  list(
    partitions = sample_ρ,
    jumps = sample_j,
    npart = sample_k,
    entr  = sample_h,
    phi_chain = sample_φ[-1],
    pchange = cont/R,
    n = n,
    R = R
  )
}

## Monitoraggio della convergenza ####

conv <- function(change_det,x=NULL, burn=floor(change_det$R*.4),
                 band=T,cut=F,BW=T) {
  if (is.null(x)) x <- 1:change_det$n
  if (BW) clrs <- c('black','black','gray','black')#gray30
  else clrs <- c('black','blue','red','green')
  par(mfrow=c(2,1),mar = c(1.5, 4, 1, 1))
  
  R <- change_det$R#length(change_det$npart)
  asc <- seq(1,R) #ascissa
  plot(x=rep(asc,times=sapply(change_det$jumps,length)),y=x[unlist(change_det$jumps)],
       type='p',xlab="",ylab="Punti di cambio",xaxt="n", #,main="Convergence monitoring"
       pch=15,cex=40/length(x), xlim=c(0,R), ylim=c(min(x),max(x)), col=clrs[1]) #cex=.5
  abline(h=c(min(x),max(x)),col=clrs[2],lwd=1.5)
  #if (cut) abline(v=burn,col=clrs[3],lwd=2,lty=1)
  if (cut) {
    rect(0,min(x),burn,max(x), col = rgb(0, 0, 0, 0.1), border = NA)
    abline(v=burn,lwd=1.5,lty=3)
  }
  
  par(mar = c(4, 4, 0, 1))
  
  plot(change_det$phi_chain, type='l',xlab="Iterazioni",ylab=expression(phi), lwd=1.5, col=clrs[1])
  abline(h=c(0,1),col=clrs[2])
  if (band) abline(h=quantile(change_det$phi_chain[burn:R],probs=c(.1,.9)),col=clrs[4],lty=2,lwd=1.5)
  #if (cut) abline(v=burn,col=clrs[3],lwd=2,lty=1)
  if (cut) {
    rect(0,0,burn,1, col = rgb(0, 0, 0, 0.1), border = NA)
    abline(v=burn,lwd=1.5,lty=3)
  }
  
  par(mfrow=c(1,1))
}

## Estrazione della partizione ottimale ####

get_opt_part <- function(change_det,burn=floor(change_det$R*.4),iter=TRUE) {
  n <- change_det$n
  
  if (burn==0) {
    parts <- change_det$partitions
    chain <- change_det$phi_chain
  } else {
    parts <- change_det$partitions[-(1:burn)]
    chain <- change_det$phi_chain[-(1:burn)]
  }
  
  R = length(parts)
  
  # calcolo delle matrici di similarità
  sample_W <- list()
  
  for (r in 1:R) {if (iter) print(paste0(r,"/",R))
    ρ <- parts[[r]]
    W <- diag(rep(1,n))
    for (i in 2:n) {
      for (j in 1:(i-1)) {
        starts <- cumsum(ρ)+1
        W[i,j] <- which(starts>i)[1] == which(starts>j)[1]
      }
    }
    sample_W[[r]] <- W
  }
  
  # matrice di similarità globale a posteriori
  M <- Reduce("+",sample_W)/R
  
  # scarti quadratici
  scarti <- numeric(1)
  for (r in 1:R) {
    scarti[r] <- sum((sample_W[[r]]-M)^2)
  }
  
  # partizione ottimale
  part <- parts[[which.min(scarti)]]
  
  # lista dei salti
  starts <- head(cumsum(part)+1, -1)
  
  # output
  list(
    part = part,
    jumps = starts,
    phi = mean(chain),
    n_kept = R
  )
}

## Metodo grafico ####

graph <- function(opt_part,y,x=1:length(y),BW=T) { #jumps=T
  if (BW) clrs <- c('black','black')
  else clrs <- c('blue','red')
    
  plot(y=y,x=x,type='b',cex=1,pch=1)#pch=19
  
  jumps <- opt_part$jumps
  if (length(jumps) > 0) {
    points(x=c(x[jumps-1],x[jumps]), y=c(y[jumps-1],y[jumps]),lwd=1.5,pch=3, col=clrs[2])
    rect(x[jumps-1], par("usr")[3], x[jumps], par("usr")[4],
         col = rgb(0, 0, 0, 0.1), border = NA)
  }
}
