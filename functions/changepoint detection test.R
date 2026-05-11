setwd(".../functions")
source("changepoint detection.R")

# SIMULAZIONE ####

## AR(0) ####
set.seed(1)

y_sim <- c(
  rnorm(40,0,1),
  rnorm(40,4,1),
  rnorm(40,1,1)
)

## AR(1) ####
set.seed(1)

y_sim <- rnorm(1, 0, 1)
for(i in 2:40){
  y_sim[i] <- rnorm(1, 0 * 0.7 + 0.3 * y_sim[i-1], 1)
}
for(i in 41:80){
  y_sim[i] <- rnorm(1, 4 * 0.7 + 0.3 * y_sim[i-1], 1)
}
for(i in 81:120){
  y_sim[i] <- rnorm(1, 1 * 0.7 + 0.3 * y_sim[i-1], 1)
}

# ESITO ####
set.seed(1)

sim <- changepoint_detection(y_sim, R=500)
conv(sim,cut=T)
sim$pchange
ris <- get_opt_part(sim)
graph(ris,y_sim)
ris$phi

