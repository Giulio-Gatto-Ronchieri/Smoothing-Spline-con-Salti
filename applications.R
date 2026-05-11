# CARICAMENTO FUNZIONI ####
setwd(".../ssj functions")

## Funzione C++ ####
Rcpp::sourceCpp("llt.cpp")


## Funzioni R ####
source('changepoint detection.R')

### ssj ####
source('ssj_fix.R')
source('ssj_mle.R')
source('auto_ssj.R')

### ssjAlt ####
source('ssjAlt_fix.R')
source('ssjAlt_mle.R')
source('auto_ssjAlt.R')

### ssjExp ####
source('ssjExp_fix.R')
source('ssjExp_mle.R')
source('auto_ssjExp.R')



# CARICAMENTO DATI ####
library(openxlsx)
setwd(".../data")

## Dati con variabile dipendente discreta ####

### senza salti ####
#### debito pubblico italiano ####
#https://db.nomics.world/Eurostat/gov_10q_ggdebt/Q.F31.S13.MIO_EUR.IT?tab=chart
debito <- read.csv("debt.csv")
colnames(debito) <- c('periodo', 'indice')
debito <- debito[1:100,]

### con salti ####
#### produzione industriale italiana ####
#https://db.nomics.world/Eurostat/sts_inpr_q/Q.PRD.B-D.CA.I21.IT?tab=chart
prod <- read.csv("production.csv")
colnames(prod) <- c('periodo', 'indice')
prod <- prod[(10*4+1):(35*4),]

#### indicatore dei prezzi dei consumatori italiani ####
#https://db.nomics.world/Eurostat/prc_hicp_midx/M.I05.CP00.IT?tab=chart
prezzi <- read.csv("price.csv")
colnames(prezzi) <- c('periodo', 'indice')
prezzi <- prezzi[which(prezzi$periodo=='2019-01'):which(prezzi$periodo=='2024-12'),]

#### produzione di bestiame italiana ####
#(unici dati disponibili dal 1971 al 2025: bull, bullock e cow)
#https://agridata.ec.europa.eu/extensions/DashboardBeef/BeefProduction.html
anni <- 1971:2025#as.character(1971:2025)

bull <- read.xlsx("bull.xlsx")[,c(1,3)]
colnames(bull) <- c('anno', 'peso') #in migliaia di tonnellate
bull <- bull[bull$anno%in%anni,]

bullock <- read.xlsx("bullock.xlsx")[,c(1,3)]
colnames(bullock) <- c('anno', 'peso') #in migliaia di tonnellate
bullock <- bullock[bullock$anno%in%anni,]

cow <- read.xlsx("cow.xlsx")[,c(1,3)]
colnames(cow) <- c('anno', 'peso') #in migliaia di tonnellate
cow <- cow[cow$anno%in%anni,]


## Dati con variabile dipendente continua ####

### senza salti ####
#### iris ####
#libreria R base

### con salti ####
#### mcycle ####
#libreria 'MASS'
library(MASS)

#### old faithful ####
#libreria R base

#### prezzi nel mercato dei cereali ####
#(EU weekly cereals prices – detailed > Feed maize, Bologna)
#https://agriculture.ec.europa.eu/data-and-analysis/markets/overviews/market-observatories/crops/cereals-statistics_en
cereals = read.xlsx("cereals-eu-prices.xlsx", sheet = "Data", startRow = 2)



# ANALISI ####
save <- TRUE #salvataggio grafici
setwd("...")

## Dati con variabile dipendente discreta ####

### senza salti ####
#### debito pubblico italiano ####
plot(debito$indice,type='l')

x <- seq(1,nrow(debito))
y <- debito$indice

ndx_2008 <- 8*4+1
ndx_2010 <- 10*4+1
ndx_2012 <- 12*4+1
ndx_2014 <- 14*4+1
ndx_2020 <- 20*4+1

# serie storica
if (save) pdf("debito.pdf",width=6, height=3)
par(mar = c(4, 4, 1.2, 1.2))
plot(x,y,type='l',xaxt = "n",xlab='Anno',ylab='Denaro')
rect(ndx_2008, par("usr")[3], ndx_2010-2, par("usr")[4],
     col = rgb(0, 0, 0, 0.1), border = NA)
rect(ndx_2012, par("usr")[3], ndx_2014-2, par("usr")[4],
     col = rgb(0, 0, 0, 0.1), border = NA)
rect(ndx_2020+1, par("usr")[3], ndx_2020+2, par("usr")[4],
     col = rgb(0, 0, 0, 0.1), border = NA)
axis(1, seq(1,length(x)+1, by=20), seq(2000,2025, by=5))
if (save) dev.off()

# salti
set.seed(1)
sim <- changepoint_detection(y, R=500)
if (save) pdf("debito_bcd.pdf",width=6, height=3)
conv(sim,x=x) # -> nessun salto
if (save) dev.off()

ris <- get_opt_part(sim)
round(ris$phi,2) # -> phi non alto (superiore a .9)

# smoothing
lambda <- 5
smooth <- smooth.spline(x,y,lambda = lambda/diff(range(x))^3)
ssj <- auto_ssj_fix(y,x,lambda=lambda)

if (save) pdf("debito_smooth1.pdf",width=6, height=3)
par(mar = c(2, 2, 1.2, 1.2))
plot(x,y,type='l',xaxt = "n",xlab='',ylab='',lty=1,lwd=1)
#plot(x,y,type='l',xaxt = "n",xlab='',ylab='',lty=3,lwd=2)#gray30
axis(1, seq(1,length(x)+1, by=16), seq(2000,2024, by=4))

lines(smooth$x,smooth$y,col='black',lwd=1.5)

lines(x=x,y=ssj$smoothed_level,col='black',lty=2,lwd=3)
#abline(v=x[which(ssj$sigmas>0)],lty=3)
if (save) dev.off()
ssj$sigmas # -> salti erroneamente identificati

# uguale al filtro HP:
#install.packages("mFilter")
library(mFilter)
hp <- hpfilter(x=y,freq=5)

par(mar = c(2, 2, 1.2, 1.2))
plot(x,y,type='l',xaxt = "n",xlab='',ylab='',lty=1,lwd=1)
axis(1, seq(1,length(x)+1, by=16), seq(2000,2024, by=4))
lines(smooth$x,smooth$y,col='black',lwd=1.5)
lines(x=x,y=hp$trend,col='black',lty=2,lwd=3)

# con salti disattivati
ssj <- ssj_fix(y,x,lambda=lambda,maxsum=0)

plot(x,y,type='l',xaxt = "n",xlab='Anno',ylab='Indice',col='black',lty=3,lwd=2)
axis(1, seq(1,length(x)+1, by=16), seq(2000,2024, by=4))
lines(smooth$x,smooth$y,col='black',lwd=1.5)
lines(x=x,y=ssj$smoothed_level,col='black',lty=2,lwd=3)

# lambda stimato
smooth <- smooth.spline(x,y)
ssj <- auto_ssj_mle(y,x)

if (save) pdf("debito_smooth2.pdf",width=6, height=3)
par(mar = c(2, 2, 1.2, 1.2))
plot(x,y,type='l',xaxt = "n",xlab='',ylab='',col='black',lty=3,lwd=2)
axis(1, seq(1,length(x)+1, by=16), seq(2000,2024, by=4))
lines(smooth$x,smooth$y,col='black',lwd=1.5)
lines(x=x,y=ssj$smoothed_level,col='black',lty=2,lwd=2)
if (save) dev.off()

# metodi alternativi
#ssjAlt <- auto_ssjAlt_fix(y,x,lambda=lambda)
#ssjExp <- auto_ssjExp_fix(y,x,lambda=lambda)
ssjAlt <- auto_ssjAlt_mle(y,x)
ssjExp <- auto_ssjExp_mle(y,x)

plot(x,y,type='l',xaxt = "n",xlab='Anno',ylab='Indice',col='black',lty=3,lwd=2)
axis(1, seq(1,length(x)+1, by=16), seq(2000,2024, by=4))
lines(x=x,y=ssjAlt$smoothed_level,col='black',lwd=1.5)
lines(x=x,y=ssjExp$smoothed_level,col='black',lty=2,lwd=3)

ssjExp$gammas

## BIC
lam <- seq(lambda/10,lambda*100, length.out=40)
out <- numeric(0)
for (i in 1:40) {
  print(paste0(i,"/",40))
  ssj <- auto_ssj_fix(y,x, lambda=lam[i])
  out[i] <- ssj$ic["bic"]
}
lam[which.min(out)]
plot(x=lam,y=out, type='b')
abline(h=min(out),v=lam[which.min(out)],col='red',lwd=2)
abline(v=lambda,col='blue',lwd=2)
##


### con salti ####
#### produzione industriale italiana ####
plot(prod$indice,type='l',xaxt = "n",xlab='Anno',ylab='Indice')
axis(1, seq(1,nrow(prod)+1, by=20), seq(2000,2025, by=5))

x <- seq(1,nrow(prod))
y <- prod$indice

#
ndx_2008 <- 8*4+1
ndx_2010 <- 10*4+1
ndx_2012 <- 12*4+1
ndx_2014 <- 14*4+1
ndx_2020 <- 20*4+1
#
ndx_2011 <- 11*4+1
ndx_2019 <- 19*4+1

if (save) pdf("industria.pdf",width=6, height=3)
par(mar = c(4, 4, 2, 2))
plot(x,y,type='l',xaxt = "n",xlab='Anno',ylab='Indice')
rect(ndx_2008+1, par("usr")[3], ndx_2010-2, par("usr")[4],
     col = rgb(0, 0, 0, 0.1), border = NA)
rect(ndx_2011+1, par("usr")[3], ndx_2011+3, par("usr")[4], #ndx_2012+2
     col = rgb(0, 0, 0, 0.1), border = NA)
rect(ndx_2019+3, par("usr")[3], ndx_2020+3, par("usr")[4],
     col = rgb(0, 0, 0, 0.1), border = NA)
axis(1, seq(1,length(x)+1, by=20), seq(2000,2025, by=5))
if (save) dev.off()

# salti
set.seed(1)
sim <- changepoint_detection(y, R=500)
conv(sim,x=x,cut=T)

ris <- get_opt_part(sim)
j <- ris$jumps

if (save) pdf("industria_bcd.pdf",width=6, height=3)
par(mar = c(2, 2, 1, 1))
plot(x,y,type='b',xaxt = "n",xlab='',ylab='',lty=1,lwd=1,cex=.5)
axis(1, seq(1,length(x)+1, by=20), seq(2000,2025, by=5))
abline(v=(x[j-1]+x[j])/2,lty=1,lwd=1.5)
if (save) dev.off()

# smoothing
smooth <- smooth.spline(x,y)
ssj <- auto_ssj(y,x)

if (save) pdf("industria_smooth.pdf",width=6, height=3)
par(mar = c(2, 2, 1, 1))
plot(x,y,type='l',xaxt = "n",xlab='',ylab='',lty=1,lwd=1)
#plot(x,y,type='l',xaxt = "n",xlab='',ylab='',lty=3,lwd=1.5)#gray30
axis(1, seq(1,length(x)+1, by=20), seq(2000,2025, by=5))

lines(smooth$x,smooth$y,col='black',lty=6,lwd=2)

lines(x=x,y=ssj$smoothed_level,col='black',lty=1,lwd=2)
rect(x[which(ssj$sigmas>0)-1], par("usr")[3], x[which(ssj$sigmas>0)], par("usr")[4],
     col = rgb(0, 0, 0, 0.1), border = NA)
if (save) dev.off()

ssj2 <- auto_ssj(y,x,alpha=1)

par(mar = c(2, 2, 1, 1))
plot(x,y,type='l',xaxt = "n",xlab='',ylab='',lty=1,lwd=1)
#plot(x,y,type='l',xaxt = "n",xlab='',ylab='',lty=3,lwd=1.5)#gray30
axis(1, seq(1,length(x)+1, by=20), seq(2000,2025, by=5))

lines(smooth$x,smooth$y,col='black',lty=6,lwd=2)

lines(x=x,y=ssj2$smoothed_level,col='black',lty=1,lwd=2)
rect(x[which(ssj2$sigmas>0)-1], par("usr")[3], x[which(ssj2$sigmas>0)], par("usr")[4],
     col = rgb(0, 0, 0, 0.1), border = NA)

smooth$lambda*diff(range(x))^3
ssj$pars['lambda']
ssj2$pars['lambda']

# detrendizzazione smoothing spline
detr_ss <- y-smooth$y
detr_ssj <- y-ssj$smoothed_level

if (save) pdf("industria_detr.pdf",width=6, height=3)
par(mfrow=c(2,1),mar = c(1, 2, 1, 1))
plot(x,detr_ss,type='l',xaxt = "n",xlab='',ylab='', ylim=c(min(detr_ss)-1,max(detr_ss)+1))
abline(h=0,lty=2,lwd=1.5)
rect(ndx_2008+1, par("usr")[3], ndx_2010-2, par("usr")[4],
     col = rgb(0, 0, 0, 0.1), border = NA)
rect(ndx_2011+1, par("usr")[3], ndx_2011+3, par("usr")[4], #ndx_2012+2
     col = rgb(0, 0, 0, 0.1), border = NA)
rect(ndx_2019+3, par("usr")[3], ndx_2020+3, par("usr")[4],
     col = rgb(0, 0, 0, 0.1), border = NA)

par(mar = c(2, 2, 0, 1))
plot(x,detr_ssj,type='l',xaxt = "n",xlab='',ylab='', ylim=c(min(detr_ss)-1,max(detr_ss)+1))
abline(h=0,lty=2,lwd=1.5)
rect(ndx_2008+1, par("usr")[3], ndx_2010-2, par("usr")[4],
     col = rgb(0, 0, 0, 0.1), border = NA)
rect(ndx_2011+1, par("usr")[3], ndx_2011+3, par("usr")[4], #ndx_2012+2
     col = rgb(0, 0, 0, 0.1), border = NA)
rect(ndx_2019+3, par("usr")[3], ndx_2020+3, par("usr")[4],
     col = rgb(0, 0, 0, 0.1), border = NA)
axis(1, seq(1,length(x)+1, by=20), seq(2000,2025, by=5))
par(mfrow=c(1,1))
if (save) dev.off()

# metodi alternativi
ssj <- auto_ssj_mle(y,x)
ssjAlt <- auto_ssjAlt(y,x)

if (save) pdf("industria_smoothAlt.pdf",width=6, height=3)
par(mar = c(2, 2, 1, 1))
plot(x,y,type='l',xaxt = "n",xlab='',ylab='',lty=1,lwd=1)
#plot(x,y,type='l',xaxt = "n",xlab='',ylab='',lty=3,lwd=1.5)#gray30
axis(1, seq(1,length(x)+1, by=20), seq(2000,2025, by=5))

lines(x=x,y=ssj$smoothed_level,col='black',lty=2,lwd=2)
#rect(x[which(ssj$sigmas>0)-1], par("usr")[3], x[which(ssj$sigmas>0)], par("usr")[4],
#     col = rgb(0, 0, 0, 0.1), border = NA)
lines(x=x,y=ssjAlt$smoothed_level,col='black',lty=1,lwd=2)
rect(x[which(ssjAlt$sigmas>0)-1], par("usr")[3], x[which(ssjAlt$sigmas>0)], par("usr")[4],
     col = rgb(0, 0, 0, 0.1), border = NA)
if (save) dev.off()

#ssj <- auto_ssjExp(y,x)
ssj <- auto_ssjExp_mle(y,x)

plot(x,y,type='l',xaxt = "n",xlab='',ylab='',lty=3,lwd=1.5)#gray30
axis(1, seq(1,length(x)+1, by=20), seq(2000,2025, by=5))

lines(x=x,y=ssj$smoothed_level,col='black',lty=1,lwd=2)

## BIC
lambda <- 30
lam <- seq(lambda/10,lambda*100, length.out=40)
out <- numeric(0)
for (i in 1:40) {
  print(paste0(i,"/",40))
  ssj <- auto_ssj_fix(y,x, lambda=lam[i])
  out[i] <- ssj$ic["bic"]
}
lam[which.min(out)]
plot(x=lam,y=out, type='b')
abline(h=min(out),v=lam[which.min(out)],col='red',lwd=2)
abline(v=lambda,col='blue',lwd=2)
##


#### indicatore dei prezzi dei consumatori italiani ####
x <- seq(1,nrow(prezzi))
y <- prezzi$indice

ndx_2021 <- 2*12+1
ndx_2022 <- 3*12+1
ndx_2023 <- 4*12+1

# grafico
if (save) pdf("prezzi.pdf",width=6, height=3)
par(mar = c(4, 4, 1, 1))
plot(x,y,type='l',xaxt = "n",xlab='Anno',ylab='Indice')
axis(1, seq(1,length(x)+1, by=12), 2019:2025)

rect(c(ndx_2021,ndx_2022), par("usr")[3], c(ndx_2022,ndx_2023), par("usr")[4],
     col = rgb(0, 0, 0, 0.1), border = NA)
abline(v=c(ndx_2021,ndx_2022,ndx_2023),lty=3)
if (save) dev.off()

smooth <- smooth.spline(x,y)
plot(x,y,type='l',xaxt = "n",xlab='Anno',ylab='Indice')
axis(1, seq(1,length(x)+1, by=12), 2019:2025)
lines(smooth$x,smooth$y,col='black',lwd=2)

# lambda fissato
lambda <- 500
smooth <- smooth.spline(x,y,lambda = lambda/diff(range(x))^3)
ssj <- auto_ssj_fix(y,x,lambda=lambda)

if (save) pdf("prezzi_smooth1.pdf",width=6, height=3)
par(mar = c(2, 2, 1, 1))
plot(x,y,type='l',xaxt = "n",xlab='Anno',ylab='Indice')
axis(1, seq(1,length(x)+1, by=12), 2019:2025)

lines(smooth$x,smooth$y,col='black',lty=6,lwd=2)
lines(x=x,y=ssj$smoothed_level,col='black',lwd=2)

rect(x[which(ssj$sigmas>0)-1], par("usr")[3], x[which(ssj$sigmas>0)], par("usr")[4],
     col = rgb(0, 0, 0, 0.1), border = NA)
if (save) dev.off()

ssj <- auto_ssjAlt_fix(y,x,lambda=lambda)
plot(x,y,type='l',xaxt = "n",xlab='Anno',ylab='Indice')
axis(1, seq(1,length(x)+1, by=12), 2019:2025)
lines(x=x,y=ssj$smoothed_level,col='black',lwd=2)

ssj <- auto_ssjExp_fix(y,x,lambda=lambda)
plot(x,y,type='l',xaxt = "n",xlab='Anno',ylab='Indice')
axis(1, seq(1,length(x)+1, by=12), 2019:2025)
lines(x=x,y=ssj$smoothed_level,col='black',lwd=2)

# lambda automatico
ssj <- auto_ssj_mle(y,x)
ssjExp <- auto_ssjExp_mle(y,x)

if (save) pdf("prezzi_smooth2.pdf",width=6, height=3)
par(mar = c(2, 2, 1, 1))
plot(x,y,type='l',xaxt = "n",xlab='Anno',ylab='Indice')
axis(1, seq(1,length(x)+1, by=12), 2019:2025)

lines(x=x,y=ssj$smoothed_level,col='black',lty=2,lwd=2)
lines(x=x,y=ssjExp$smoothed_level,col='black',lwd=2)
if (save) dev.off()

## BIC
lam <- seq(lambda/10,lambda*100, length.out=40)
out <- numeric(0)
for (i in 1:40) {
  print(paste0(i,"/",40))
  ssj <- auto_ssj_fix(y,x, lambda=lam[i])
  out[i] <- ssj$ic["bic"]
}
lam[which.min(out)]
plot(x=lam,y=out, type='b')
abline(h=min(out),v=lam[which.min(out)],col='red',lwd=2)
abline(v=lambda,col='blue',lwd=2)
##


#### allevamento di bestiame ####
##
plot(cow$anno,cow$peso, type='l',xlab='Anno',ylab='Peso')

x <- seq(1,nrow(cow))
y <- cow$peso

lambda <- 10
smooth <- smooth.spline(x,y, lambda=lambda/diff(range(x))^3)
ssj <- ssjAlt_fix(y,x, lambda=lambda,maxsum=1000)
ssj$sigmas

plot(x,y, type='l',xlab='Anno',ylab='Peso')
lines(smooth$x,smooth$y,col='black',lty=6,lwd=1.5)
lines(x=x,y=ssj$smoothed_level,col='black',lwd=2)
# -> impressionante capacità di ssjAlt di seguire i dati, tramite salti nel livello

ssj <- ssj_fix(y,x, lambda=lambda,maxsum=1000)
plot(x,y, type='l',xlab='Anno',ylab='Peso')
lines(smooth$x,smooth$y,col='black',lty=6,lwd=1.5)
lines(x=x,y=ssj$smoothed_level,col='black',lwd=2)
# -> ssj invece usa i cambiamenti di pendenza
##

plot(bull$anno,bull$peso, type='l',xlab='Anno',ylab='Peso')
plot(bullock$anno,bullock$peso, type='l',xlab='Anno',ylab='Peso')

# plot(cow$anno,bull$peso+bullock$peso, type='l',xlab='Anno',ylab='Peso')
# plot(cow$anno,bull$peso+cow$peso, type='l',xlab='Anno',ylab='Peso')
# plot(cow$anno,bullock$peso+cow$peso, type='l',xlab='Anno',ylab='Peso')
# plot(cow$anno,bull$peso+bullock$peso+cow$peso, type='l',xlab='Anno',ylab='Peso')

x <- seq(1,nrow(bull))
y <- bull$peso
y2 <- bullock$peso

plot(x,y+y2,type='l')

# grafico
if (save) pdf("best.pdf",width=6, height=3)
par(mar = c(4, 4, 1, 1))
plot(x,y,type='l',xaxt = "n",xlab='Anno',ylab='Peso', ylim=c(0,750))
lines(x,y2,lty=5)
axis(1, seq(1,length(x)+1, by=9), seq(1971,2025, by=9))
rect(1975-1970, par("usr")[3], 1976-1970, par("usr")[4],
     col = rgb(0, 0, 0, 0.1), border = NA)
abline(v=c(2007,2014)-1970,lty=3,lwd=1.5)
if (save) dev.off()

smooth <- smooth.spline(x,y)
plot(x,y,type='l',xaxt = "n",xlab='Anno',ylab='Indice')
axis(1, seq(1,length(x)+1, by=9), seq(1971,2025, by=9))
lines(smooth$x,smooth$y,col='black',lwd=2)

# lambda fissato
lambda <- 10
smooth <- smooth.spline(x,y, lambda=lambda/diff(range(x))^3)
ssj <- auto_ssj_fix(y,x, lambda=lambda)

if (save) pdf("bull_ssj_fix.pdf",width=6, height=3)
par(mar = c(2, 2, 1, 1))
plot(x,y, type='l',xlab='',ylab='', xaxt='n')
axis(1, seq(1,length(x)+1, by=9), seq(1971,2025, by=9))
lines(smooth$x,smooth$y,col='black',lty=6,lwd=2)
lines(x=x,y=ssj$smoothed_level,col='black',lwd=2)
rect(x[which(ssj$sigmas>0)-1], par("usr")[3], x[which(ssj$sigmas>0)], par("usr")[4],
     col = rgb(0, 0, 0, 0.1), border = NA)
if (save) dev.off()

ssjAlt <- auto_ssjAlt_fix(y,x, lambda=lambda)

plot(x,y, type='l',xlab='',ylab='', xaxt='n')
axis(1, seq(1,length(x)+1, by=9), seq(1971,2025, by=9))
lines(smooth$x,smooth$y,col='black',lty=6,lwd=1.5)
lines(x=x,y=ssjAlt$smoothed_level,col='black',lwd=2)
rect(x[which(ssjAlt$sigmas>0)-1], par("usr")[3], x[which(ssjAlt$sigmas>0)], par("usr")[4],
     col = rgb(0, 0, 0, 0.1), border = NA)

# lambda stimato
#smooth <- smooth.spline(x,y)
ssj <- auto_ssj_mle(y,x)

plot(x,y, type='l',xlab='',ylab='', xaxt='n')
axis(1, seq(1,length(x)+1, by=9), seq(1971,2025, by=9))
lines(smooth$x,smooth$y,col='black',lty=6,lwd=1.5)
lines(x=x,y=ssj$smoothed_level,col='black',lwd=2)
rect(x[which(ssj$sigmas>0)-1], par("usr")[3], x[which(ssj$sigmas>0)], par("usr")[4],
     col = rgb(0, 0, 0, 0.1), border = NA)

ssj <- auto_ssjAlt_mle(y,x)

plot(x,y, type='l',xlab='',ylab='', xaxt='n')
axis(1, seq(1,length(x)+1, by=9), seq(1971,2025, by=9))
lines(smooth$x,smooth$y,col='black',lty=6,lwd=1.5)
lines(x=x,y=ssj$smoothed_level,col='black',lwd=2)
rect(x[which(ssj$sigmas>0)-1], par("usr")[3], x[which(ssj$sigmas>0)], par("usr")[4],
     col = rgb(0, 0, 0, 0.1), border = NA)

ssj <- auto_ssjExp_mle(y,x)

plot(x,y, type='l',xlab='',ylab='', xaxt='n')
axis(1, seq(1,length(x)+1, by=9), seq(1971,2025, by=9))
lines(smooth$x,smooth$y,col='black',lty=6,lwd=1.5)
lines(x=x,y=ssj$smoothed_level,col='black',lwd=2)
rect(x[which(ssj$gammas>0)-1], par("usr")[3], x[which(ssj$gammas>0)], par("usr")[4],
     col = rgb(0, 0, 0, 0.1), border = NA)

## BIC
lam <- seq(lambda/10,lambda*100, length.out=40)
out <- numeric(0)
for (i in 1:40) {
  print(paste0(i,"/",40))
  ssj <- auto_ssj_fix(y,x, lambda=lam[i])
  out[i] <- ssj$ic["bic"]
}
lam[which.min(out)]
plot(x=lam,y=out, type='b')
abline(h=min(out),v=lam[which.min(out)],col='red',lwd=2)
abline(v=lambda,col='blue',lwd=2)
##


## Dati con variabile dipendente continua ####

### senza salti ####
#### iris ####
if (save) pdf("iris1.pdf",width=6, height=3)
par(mfrow=c(2,3),mar=c(1,1,1,1))
for (i in 1:3) {
  for (j in (i+1):4) {
    x <- iris[,i]
    y <- iris[,j]
    
    plot(x,y, xlab='',ylab='', xaxt='n',yaxt='n', cex=.5)
    
    smooth <- smooth.spline(x,y)
    lines(smooth$x,smooth$y,col='black',lwd=3.5)#,lwd=1.5
    
    ssj <- auto_ssj_fix(y,x,lambda=smooth$lambda*diff(range(x))^3)
    #plot(x,y)
    #plot(x=sort(x),y=ssj$smoothed_level,type='l')
    lines(x=ssj$x,y=ssj$smoothed_level,col='gray',lwd=2)#,lty=2,lwd=2.5
    # abline(v=x[which(ssj$sigmas>0)],lty=3)
    print(round(sum(ssj$sigmas[-1]>0)/150,2))
  }
}
if (save) dev.off()

if (save) pdf("iris2.pdf",width=6, height=3)
par(mfrow=c(2,3),mar=c(1,1,1,1))
for (i in 1:3) {
  for (j in (i+1):4) {
    x <- iris[,i]*100
    y <- iris[,j]
    
    plot(x,y, xlab='',ylab='', xaxt='n',yaxt='n', cex=.5)
    
    smooth <- smooth.spline(x,y)
    lines(smooth$x,smooth$y,col='black',lwd=3.5)
    
    ssj <- auto_ssj_fix(y,x,lambda=smooth$lambda*diff(range(x))^3)
    #plot(x,y)
    #plot(x=sort(x),y=ssj$smoothed_level,type='l')
    lines(x=ssj$x,y=ssj$smoothed_level,col='gray',lwd=2)
    # abline(v=x[which(ssj$sigmas>0)],lty=3)
    print(round(sum(ssj$sigmas[-1]>0)/150,2))
  }
}
if (save) dev.off()

par(mfrow=c(1,1))
# sensibile alla scala: la tendenza all'identificare dei salti aumenta
# all'avvicinarsi delle osservazioni

## BIC
x <- iris[,1]
y <- iris[,2]

lambda <- smooth.spline(x,y)$lambda*diff(range(x))^3
lam <- seq(lambda/10,lambda*100, length.out=40)
out <- numeric(0)
for (i in 1:40) {
  print(paste0(i,"/",40))
  ssj <- auto_ssj_fix(y,x, lambda=lam[i])
  out[i] <- ssj$ic["bic"]
}
lam[which.min(out)]
plot(x=lam,y=out, type='b')
abline(h=min(out),v=lam[which.min(out)],col='red',lwd=2)
abline(v=lambda,col='blue',lwd=2)
#
x <- iris[,1]*100 #20

lambda <- smooth.spline(x,y)$lambda*diff(range(x))^3
lam <- seq(lambda/10,lambda*100, length.out=40)
out <- numeric(0)
for (i in 1:40) {
  print(paste0(i,"/",40))
  ssj <- auto_ssj_fix(y,x, lambda=lam[i])
  out[i] <- ssj$ic["bic"]
}
lam[which.min(out)]
plot(x=lam,y=out, type='b')
abline(h=min(out),v=lam[which.min(out)],col='red',lwd=2)
abline(v=lambda,col='blue',lwd=2)
##


### con salti ####
#### mcycle ####
x <- mcycle$times
y <- mcycle$accel

smooth <- smooth.spline(x,y)

if (save) pdf("mcycle.pdf",width=6, height=3)
par(mar=c(4,4,1,1))
plot(x,y, xlab='Tempo',ylab='Accelerazione', cex=.7)
#rect(c(14.2,32), par("usr")[3], c(22.7,par("usr")[2]), par("usr")[4],
rect(c(0,22.7), par("usr")[3], c(14.2,32), par("usr")[4],
     col = rgb(0, 0, 0, 0.1), border = NA)
abline(v=c(14.2,22.7,32,par("usr")[2]),lty=3,lwd=2)
lines(smooth$x,smooth$y,col='black',lwd=1.5)
if (save) dev.off()

lambda <- smooth$lambda*diff(range(x))^3
ssj <- auto_ssj_fix(y,x,lambda=lambda)

if (save) pdf("mcycle_smooth.pdf",width=6, height=3)
par(mar=c(2,2,1,1))
plot(x,y, xlab='',ylab='', cex=.7)
lines(x=ssj$x,y=ssj$smoothed_level,col='black',lwd=1.5)
rect(x[which(ssj$sigmas>0)-1], par("usr")[3], x[which(ssj$sigmas>0)], par("usr")[4],
     col = rgb(0, 0, 0, 0.1), border = NA)
abline(v=c(14.2,22.7,32,par("usr")[2]),lty=3,lwd=2)
if (save) dev.off()

# altri metodi
ssj <- auto_ssj_mle(y,x)
plot(x,y, xlab='',ylab='', cex=.7)
lines(x=ssj$x,y=ssj$smoothed_level,col='black',lwd=1.5)

ssj <- auto_ssj(y,x)
plot(x,y, xlab='',ylab='', cex=.7)
lines(x=ssj$x,y=ssj$smoothed_level,col='black',lwd=1.5)

#
ssj <- auto_ssjAlt_fix(y,x,lambda=lambda)
plot(x,y, xlab='',ylab='', cex=.7)
lines(x=ssj$x,y=ssj$smoothed_level,col='black',lwd=1.5)

ssj <- auto_ssjAlt_mle(y,x)
plot(x,y, xlab='',ylab='', cex=.7)
lines(x=ssj$x,y=ssj$smoothed_level,col='black',lwd=1.5)

#
ssj <- auto_ssjExp_fix(y,x,lambda=lambda)
plot(x,y, xlab='',ylab='', cex=.7)
lines(x=ssj$x,y=ssj$smoothed_level,col='black',lwd=1.5)

ssj <- auto_ssjExp_mle(y,x)
if (save) pdf("mcycle_smoothExp.pdf",width=6, height=3)
par(mar=c(2,2,1,1))
plot(x,y, xlab='',ylab='', cex=.7)
lines(x=ssj$x,y=ssj$smoothed_level,col='black',lwd=1.5)
if (save) dev.off()

## BIC
lam <- seq(lambda/10,lambda*100, length.out=40)
out <- numeric(0)
for (i in 1:40) {
  print(paste0(i,"/",40))
  ssj <- auto_ssj_fix(y,x, lambda=lam[i])
  out[i] <- ssj$ic["bic"]
}
lam[which.min(out)]
plot(x=lam,y=out, type='b')
abline(h=min(out),v=lam[which.min(out)],col='red',lwd=2)
abline(v=lambda,col='blue',lwd=2)
##


#### old faithful ####
x <- faithful$eruptions
y <- faithful$waiting

if (save) pdf("faithful_hist.pdf",width=6, height=2.5)
par(mar=c(4,4,1,1))
hist(x,freq=F, main='',xlab='Durata',ylab='Densità')#,breaks=7
if (save) dev.off()

lim1 <- 2.95; lim2 <- 3.27

ndx1 <- which(x<lim1); ndx2 <- which(x>lim2)
x1 <- x[ndx1]; y1 <- y[ndx1]
x2 <- x[ndx2]; y2 <- y[ndx2]

reg1 <- lm(y1 ~ x1); reg2 <- lm(y2 ~ x2)

smooth <- smooth.spline(x,y)
lambda <- smooth$lambda*diff(range(x))^3

if (save) pdf("faithful.pdf",width=6, height=3)
par(mar=c(4,4,1,1))
plot(x,y, xlab='Durata',ylab='Attesa', cex=.7)
rect(lim1, par("usr")[3], lim2, par("usr")[4],
     col = rgb(0, 0, 0, 0.1), border = NA)
abline(v=c(lim1,lim2),lty=3,lwd=1.5)
curve(coef(reg1)[1] +coef(reg1)[2]*x, from=min(x1), to=lim1,
      add=TRUE, lty=2,lwd=2)
curve(coef(reg2)[1] +coef(reg2)[2]*x, from=lim2, to=max(x2),
      add=TRUE, lty=2,lwd=2)
segments(lim1,coef(reg1)[1] +coef(reg1)[2]*lim1,
         lim2,coef(reg2)[1] +coef(reg2)[2]*lim2,
         lty=2,lwd=2)
lines(smooth$x,smooth$y,col='black',lwd=2)
if (save) dev.off()

# ssj
ssj <- auto_ssj_fix(y,x,lambda=lambda)
plot(x,y, xlab='Durata',ylab='Attesa', cex=.7)
lines(ssj$x,y=ssj$smoothed_level,col='black',lty=2,lwd=2.5)

ssj1 <- auto_ssj_mle(y,x)
plot(x,y, xlab='Durata',ylab='Attesa', cex=.7)
lines(ssj1$x,y=ssj1$smoothed_level,col='black',lty=2,lwd=2.5)
# -> uguale

ssj <- auto_ssj(y,x)
plot(x,y, xlab='Durata',ylab='Attesa', cex=.7)
lines(ssj$x,y=ssj$smoothed_level,col='black',lty=2,lwd=2.5)
# -> uguale

# ssjAlt
ssj2 <- auto_ssjAlt_fix(y,x,lambda=lambda)
plot(x,y, xlab='Durata',ylab='Attesa', cex=.7)
lines(ssj2$x,y=ssj2$smoothed_level,col='black',lty=2,lwd=2.5)

ssj <- auto_ssjAlt_mle(y,x)
plot(x,y, xlab='Durata',ylab='Attesa', cex=.7)
lines(ssj$x,y=ssj$smoothed_level,col='black',lty=2,lwd=2.5)

ssj <- auto_ssjAlt(y,x)
plot(x,y, xlab='Durata',ylab='Attesa', cex=.7)
lines(ssj$x,y=ssj$smoothed_level,col='black',lty=2,lwd=2.5)

# ssjExp
ssj <- auto_ssjExp_fix(y,x,lambda=lambda)
plot(x,y, xlab='Durata',ylab='Attesa', cex=.7)
lines(ssj$x,y=ssj$smoothed_level,col='black',lty=2,lwd=2.5)

ssj <- auto_ssjExp_mle(y,x)
plot(x,y, xlab='Durata',ylab='Attesa', cex=.7)
lines(ssj$x,y=ssj$smoothed_level,col='black',lty=2,lwd=2.5)
rect(ssj$x[which(ssj$gammas>0)-1], par("usr")[3], ssj$x[which(ssj$gammas>0)], par("usr")[4],
     col = rgb(0, 0, 0, 0.1), border = NA)

ssj <- auto_ssjExp(y,x)
plot(x,y, xlab='Durata',ylab='Attesa', cex=.7)
lines(ssj$x,y=ssj$smoothed_level,col='black',lty=2,lwd=2.5)
# -> Exp: risultati insoddisfacenti

if (save) pdf("faithful_smooth.pdf",width=6, height=3)
par(mar=c(2,2,1,1))
plot(x,y, xlab='',ylab='', cex=.7)
lines(ssj1$x,y=ssj1$smoothed_level,col='black',lty=1,lwd=2)
lines(ssj2$x,y=ssj2$smoothed_level,col='black',lty=2,lwd=2)
rect(ssj1$x[which(ssj1$sigmas>0)-1], par("usr")[3], ssj1$x[which(ssj1$sigmas>0)], par("usr")[4],
     col = rgb(0, 0, 0, 0.1), border = NA)
abline(v=c(ssj2$x[which(ssj2$sigmas>0)-1],ssj2$x[which(ssj2$sigmas>0)]), lty=3,lwd=1.5)
if (save) dev.off()

# parametri regolati manualmente:
ssj1 <- ssj_fix(y,x, lambda=10,maxsum=147)
ssj2 <- auto_ssjAlt_mle(y,x, grid = seq(0, sd(y)*10, sd(y)/10))

if (save) pdf("faithful_smoothMan.pdf",width=6, height=3)
par(mar=c(2,2,1,1))
plot(x,y, xlab='',ylab='', cex=.7)
lines(ssj1$x,y=ssj1$smoothed_level,col='black',lty=1,lwd=2)
lines(ssj2$x,y=ssj2$smoothed_level,col='black',lty=2,lwd=2)
rect(ssj1$x[which(ssj1$sigmas>0)-1], par("usr")[3], ssj1$x[which(ssj1$sigmas>0)], par("usr")[4],
     col = rgb(0, 0, 0, 0.1), border = NA)
abline(v=c(ssj2$x[which(ssj2$sigmas>0)-1],ssj2$x[which(ssj2$sigmas>0)]), lty=3,lwd=1.5)
if (save) dev.off()

## BIC
lam <- seq(lambda/10,lambda*100, length.out=40)
out <- numeric(0)
for (i in 1:40) {
  print(paste0(i,"/",40))
  ssj <- auto_ssj_fix(y,x, lambda=lam[i])
  out[i] <- ssj$ic["bic"]
}
lam[which.min(out)]
plot(x=lam,y=out, type='b')
abline(h=min(out),v=lam[which.min(out)],col='red',lwd=2)
abline(v=lambda,col='blue',lwd=2)
##


#### prezzi nel mercato dei cereali ####
week <- rev(as.Date(cereals[['Date.to']], origin = "1899-12-30")) #domeniche, fine delle settimane
year <- sapply(week, function(x) substr(x, 1, 4))
years <- (min(year):max(year))[-1]
years_ndx <- sapply(years, function(x) min(which(x==year)))

bol <- apply(cereals[,which(sapply(names(cereals),function(x) grepl('Bologna',x)))],2, rev)
par(mfrow=c(2,2),mar=c(2,2,1,1))
for (i in 1:4) {
  plot(bol[,i], xlab='Anno',ylab='Prezzo',xaxt='n',
       type='b',cex=.5,pch=19)
  axis(1, years_ndx, years)
}
par(mfrow=c(1,1))
# si sceglie la serie con la discontinuità più eclatante

maize <- as.numeric(bol[,4])
plot(maize, xlab='Anno',ylab='Prezzo',xaxt='n',
     type='l')
axis(1, years_ndx, years) # -> valori mancanti

# 2021-2023
x_full <- years_ndx[which(years==2021)]:(years_ndx[which(years==2024)]-1)
y_full <- maize[x_full]

not_missing <- which(!is.na(y_full))
x <- x_full[not_missing]; y <- y_full[not_missing]

# grafico
if (save) pdf("mais.pdf",width=6, height=3)
par(mar=c(4,4,1,1))
plot(x_full,y_full, xlab='Anno',ylab='Prezzo',xaxt='n',
     type='l')
axis(1, years_ndx, years)
gap <- F
for (i in 1:length(x_full)) {
  if (is.na(y_full[i]) & !gap) {
    start <- max(x_full[1],x_full[i]-1)
    gap <- T
  } else if (!is.na(y_full[i]) & gap) {
    rect(start, par("usr")[3], x_full[i], par("usr")[4],
         col = rgb(0, 0, 0, 0.1), border = NA)
    gap <- F
  }
}
if (gap) rect(start, par("usr")[3], x_full[i], par("usr")[4],
              col = rgb(0, 0, 0, 0.1), border = NA)
if (save) dev.off()

plot(x,y, type='b',cex=.5,pch=19)
plot(x,y, type='l')

# tempo continuo
smooth <- smooth.spline(x,y)
lines(smooth$x,smooth$y,col='red',lwd=1.5)
# -> estremo undersmoothing

# smoothing con lambda scelto a mano
lambda <- 1000 #666
smooth <- smooth.spline(x,y,lambda=lambda/diff(range(x))^3)
# -> bene nelle zone stabili, ma taglia gli angoli nelle discontinuità

ssj <- auto_ssj_fix(y,x,lambda=lambda)
# -> risolve

if (save) pdf("mais_smooth.pdf",width=6, height=3)
par(mar=c(2,2,1,1))
plot(x,y, xlab='',ylab='',xaxt='n',
     type='l')
axis(1, years_ndx, years)
lines(smooth$x,smooth$y,lty=6,lwd=2)
lines(x=x,y=ssj$smoothed_level,lwd=2)
if (save) dev.off()

if (save) pdf("mais_smoothComp.pdf",width=6, height=3)
par(mar=c(1,2,1,1))
par(mfrow=c(1,3))
plot(x,y, type='l', xlab='',ylab='',xaxt='n', ylim=c(195,405))
par(mar=c(1,1,1,1))
plot(smooth$x,smooth$y,lwd=2,type='l', xlab='',ylab='',xaxt='n',yaxt='n', ylim=c(190,420))
plot(x=x,y=ssj$smoothed_level,lwd=2,type='l', xlab='',ylab='',xaxt='n',yaxt='n', ylim=c(190,420))
par(mfrow=c(1,1))
if (save) dev.off()

# altri metodi
ssj <- auto_ssj_mle(y,x)
plot(x,y, xlab='',ylab='',xaxt='n',
     type='l')
axis(1, years_ndx, years)
lines(x=x,y=ssj$smoothed_level,lwd=2)
abline(v=x[which(ssj$sigmas>0)],lty=3)

ssj <- auto_ssjAlt_fix(y,x,lambda=lambda)
plot(x,y, xlab='',ylab='',xaxt='n',
     type='l')
axis(1, years_ndx, years)
lines(x=x,y=ssj$smoothed_level,lwd=2)
abline(v=x[which(ssj$sigmas>0)],lty=3)

ssj <- auto_ssjAlt_mle(y,x)
plot(x,y, xlab='',ylab='',xaxt='n',
     type='l')
axis(1, years_ndx, years)
lines(x=x,y=ssj$smoothed_level,lwd=2)
abline(v=x[which(ssj$sigmas>0)],lty=3)

ssj <- auto_ssjExp_fix(y,x,lambda=lambda)
plot(x,y, xlab='',ylab='',xaxt='n',
     type='l')
axis(1, years_ndx, years)
lines(x=x,y=ssj$smoothed_level,lwd=2)
abline(v=x[which(ssj$gammas>0)],lty=3)

ssj <- auto_ssjExp_mle(y,x)
plot(x,y, xlab='',ylab='',xaxt='n',
     type='l')
axis(1, years_ndx, years)
lines(x=x,y=ssj$smoothed_level,lwd=2)
abline(v=x[which(ssj$gammas>0)],lty=3)

# -> gli altri metodi tendono all'undersmoothing, come SS

# tempo discreto
ssj <- auto_ssj_fix(y_full,x_full,lambda=lambda)
#ssj <- auto_ssj_fix(y_full[-1],x_full[-1],lambda=lambda)
# lly old: la prima osservazione non può essere mancante (stime nel futuro più verosimili)

if (save) pdf("mais_smoothDisc.pdf",width=6, height=3)
par(mar=c(2,2,1,1))
plot(x_full,y_full, xlab='',ylab='',xaxt='n',
     type='l')
axis(1, years_ndx, years)
lines(x=ssj$x,y=ssj$smoothed_level,lwd=2)
gap <- F
for (i in 1:length(x_full)) {
  if (is.na(y_full[i]) & !gap) {
    start <- max(x_full[1],x_full[i]-1)
    gap <- T
  } else if (!is.na(y_full[i]) & gap) {
    rect(start, par("usr")[3], x_full[i], par("usr")[4],
         col = rgb(0, 0, 0, 0.1), border = NA)
    gap <- F
  }
}
if (gap) rect(start, par("usr")[3], x_full[i], par("usr")[4],
              col = rgb(0, 0, 0, 0.1), border = NA)
if (save) dev.off()

# altri metodi
ssj <- auto_ssj_mle(y_full,x_full)
plot(x_full,y_full, xlab='',ylab='',xaxt='n',
     type='l')
axis(1, years_ndx, years)
lines(x=ssj$x,y=ssj$smoothed_level,lwd=2)

#
ssj <- auto_ssjAlt_fix(y_full,x_full,lambda=lambda)
plot(x_full,y_full, xlab='',ylab='',xaxt='n',
     type='l')
axis(1, years_ndx, years)
lines(x=ssj$x,y=ssj$smoothed_level,lwd=2)

ssj <- auto_ssjAlt_mle(y_full,x_full)
plot(x_full,y_full, xlab='',ylab='',xaxt='n',
     type='l')
axis(1, years_ndx, years)
lines(x=ssj$x,y=ssj$smoothed_level,lwd=2)

#
ssj <- auto_ssjExp_fix(y_full,x_full,lambda=lambda)
plot(x_full,y_full, xlab='',ylab='',xaxt='n',
     type='l')
axis(1, years_ndx, years)
lines(x=ssj$x,y=ssj$smoothed_level,lwd=2)

ssj <- auto_ssjExp_mle(y_full,x_full)
plot(x_full,y_full, xlab='',ylab='',xaxt='n',
     type='l')
axis(1, years_ndx, years)
lines(x=ssj$x,y=ssj$smoothed_level,lwd=2)

## BIC
lam <- seq(lambda/10,lambda*100, length.out=40)
out <- numeric(0)
for (i in 1:40) {
  print(paste0(i,"/",40))
  ssj <- auto_ssj_fix(y,x, lambda=lam[i])
  out[i] <- ssj$ic["bic"]
}
lam[which.min(out)]
plot(x=lam,y=out, type='b')
abline(h=min(out),v=lam[which.min(out)],col='red',lwd=2)
abline(v=lambda,col='blue',lwd=2)
##
