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
                 covariance,
                 const_modif) {

    times <- unique(exit[to != 0])
    times <- times[times > s & times <= t]

    zzz <- .Call("gen_msm", times, entry, exit, from, to, nstate, const_modif)

    if (covariance) {

        cov_etm <- .Call("cov_aj",
                         zzz$time,
                         zzz$est,
                         zzz$n.risk,
                         zzz$n.event,
                         zzz$dna)

        zzz$cov <- cov_etm
        
    }
    
    zzz
}
    
