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

    c_modif <- matrix(const_modif, length(times), nstate, byrow = TRUE)

    zzz <- .Call("gen_msm", times, entry, exit, from, to, nstate, c_modif)

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
    
