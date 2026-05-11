# CARICAMENTO FUNZIONI ####
setwd(".../functions")

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

### funzioni di simulazione ####
setwd("...")
source('sim func.R')



# ESECUZIONE ####

## Esempio minimo di dati simulati ####
esemp(save="esempio.pdf")


## Benchmark: bayesian changepoint detection ####
# esempio minimo
esemp_bayes(save="esemp_bayes.pdf")
# simulazione
bayes(φ=.8,save="bayes.pdf")


## ssj ####
# ssj_fix
simulation(save="ssj_fix.pdf")
# ssj_mle
simulation(method='mle',save="ssj_mle.pdf")
# auto_ssj
simulation(N=100,method='auto',save="auto_ssj.pdf")


## ssjAlt ####
# ssjAlt_fix
simulation(func='alt',save="ssjAlt_fix.pdf")
# ssjAlt_mle
simulation(func='alt',method='mle',save="ssjAlt_mle.pdf")
# auto_ssjAlt
simulation(N=100,func='alt',method='auto',save="auto_ssjAlt.pdf")


## ssjExp ####
# ssjExp_fix
simulation(func='exp',save="ssjExp_fix.pdf")
# ssjExp_mle
simulation(func='exp',method='mle',save="ssjExp_mle.pdf")
# auto_ssjExp
simulation(N=100,func='exp',method='auto',save="auto_ssjExp.pdf")

