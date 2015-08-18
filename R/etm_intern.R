#####################################################################
### The etm machinery
### Arthur Allignol <arthur.allignol@uni-ulm.de
#####################################################################


.etm <- function(entry,
                 exit,
                 from,
                 to,
                 nstate,
                 s,
                 t,
                 covariance) {

    times <- unique(exit[to != 0])
    times <- times[times > s & times <= t]

    zzz <- .Call("gen_msm", times, entry, exit, from, to, nstate)

    if (covariance) {

        stop("doing")
        var <- .Call("cov_aj",
                     zzz$time,
                     zzz$est,
                     zzz$n.risk,
                     zzz$n.event,
                     zzz$dna)
        
    }
        
    zzz
}
    
