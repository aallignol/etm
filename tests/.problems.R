### Some problems from Tobi

require(etm)
require(data.table)

dt <- data.table(read.csv("/data/Ulm/Teaching/WS_2013_14/Competing_MSM/prg/my_icu_pneu.csv"))


tail(sort(dt[, exit]))

test <- etm(dt, c("0", "1", "2"), tra_ill(), "cens", s = 0)
test <- etm(dt, c("0", "1", "2"), tra_ill(), "cens", s = 88, t = 88.1)
