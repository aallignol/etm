### Some problems from Tobi

require(etm)
require(data.table)

dt <- data.table(read.csv("/data/Ulm/Teaching/WS_2013_14/Competing_MSM/prg/my_icu_pneu.csv"))


tail(sort(dt[, exit]))

test <- etm(dt, c("0", "1", "2"), tra_ill(), "cens", s = 0)

test <- etm(dt, c("0", "1", "2"), tra_ill(), "cens", s = 87, t = 88.1)

## one n.risk negative
dt[exit > 87]


### Get there because we have covariates
test <- etm(dt, c("0", "1", "2"), tra_ill(), "cens", s = 0, strat_variable = "sex")


######################################################################
library(etm)
data(sir.cont)
# Modification for patients entering and leaving a state
# at the same date
# Change on ventilation status is considered
# to happen before end of hospital stay
sir.cont <- sir.cont[order(sir.cont$id, sir.cont$time), ]
for (i in 2:nrow(sir.cont)) {
  if (sir.cont$id[i]==sir.cont$id[i-1]) {
    if (sir.cont$time[i]==sir.cont$time[i-1]) {
      sir.cont$time[i-1] <- sir.cont$time[i-1] - 0.5
    }
  }
}
### Computation of the transition probabilities
# Possible transitions.
tra <- matrix(ncol=3,nrow=3,FALSE)
tra[1, 2:3] <- TRUE
tra[2, c(1, 3)] <- TRUE
# etm
fit.etm <- etm(sir.cont, c("0", "1", "2"), tra, "cens", s=0)

fit.etm$time

etm(sir.cont, c("0", "1", "2"), tra, "cens", s=fit.etm$time[92]) #--> unity matrix even though we have an event (0->2 transition) at time 183

subset(sir.cont,time==1)

aa <- etm(sir.cont, c("0", "1", "2"), tra, "cens", s=0, t=1) #--> unity matrix even though something happens at t=1?

aa <- etm(sir.cont, c("0", "1", "2"), tra, "cens", s=165)

subset(sir.cont,time > 165)


trprob(fit.etm, "0 0", 0) # should return 1


