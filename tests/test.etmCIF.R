### test file for etmCIF.
### Really simple tests and comparison with etm

require(etm)

data(abortion)

from <- rep(0, nrow(abortion))
to <- abortion$cause
entry <- abortion$entry
exit <- abortion$exit
id <- 1:nrow(abortion)
data <- data.frame(id, from, to, entry, exit, group = abortion$group)

## Computation of the CIFs with etm
tra <- matrix(FALSE, 4, 4)
tra[1, 2:4] <- TRUE

cif.control <- etm(data[data$group == 0, ], c("0", "1", "2", "3"),
                        tra, NULL, 0)
cif.exposed <- etm(data[data$group == 1, ], c("0", "1", "2", "3"),
                        tra, NULL, 0)


## Computation of the CIFs with etmCIF

netm <- etmCIF(Surv(entry, exit, cause != 0) ~ group, abortion,
               etype = cause, failcode = 3)

### let's do some comparisons :-)

all.equal(trprob(cif.control, "0 3"), netm[[1]]$est["0", "3", ])
all.equal(trprob(cif.control, "0 2"), netm[[1]]$est["0", "2", ])
all.equal(trprob(cif.control, "0 1"), netm[[1]]$est["0", "1", ])

all.equal(trprob(cif.exposed, "0 3"), netm[[2]]$est["0", "3", ])
all.equal(trprob(cif.exposed, "0 2"), netm[[2]]$est["0", "2", ])
all.equal(trprob(cif.exposed, "0 1"), netm[[2]]$est["0", "1", ])


all.equal(trcov(cif.control, "0 3"), netm[[1]]$cov["0 3", "0 3", ])
all.equal(trcov(cif.control, "0 2"), netm[[1]]$cov["0 2", "0 2", ])
all.equal(trcov(cif.control, "0 1"), netm[[1]]$cov["0 1", "0 1", ])

all.equal(trcov(cif.exposed, "0 3"), netm[[2]]$cov["0 3", "0 3", ])
all.equal(trcov(cif.exposed, "0 2"), netm[[2]]$cov["0 2", "0 2", ])
all.equal(trcov(cif.exposed, "0 1"), netm[[2]]$cov["0 1", "0 1", ])


netm

## test on the summary
snetm <- summary(netm)

snetm

all.equal(unname(trprob(cif.control, "0 3")), snetm[[1]][[3]]$P)
all.equal(unname(trprob(cif.control, "0 2")), snetm[[1]][[2]]$P)
all.equal(unname(trprob(cif.control, "0 1")), snetm[[1]][[1]]$P)

