#####################################################################
### The etm machinery
### Arthur Allignol <arthur.allignol@uni-ulm.de
#####################################################################


.etm <- function(entry,
                 exit,
                 from,
                 to,
                 nstate) {

    times <- sort(unique(exit[to != 0]))

    zzz <- .Call("gen_msm", times, entry, exit, from, to, nstate)
    zzz
}
    
