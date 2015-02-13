#####################################################################
### The etm machinery
### Arthur Allignol <arthur.allignol@uni-ulm.de
#####################################################################


.etm <- function(times, 
                 entry,
                 exit,
                 from,
                 to,
                 nstate,
                 s,
                 t) {

    times <- sort(unique(exit[to != 0]))
    times <- times[times > s & times <= t]

    zzz <- .Call("gen_msm", times, entry, exit, from, to, nstate)
    zzz
}
    
