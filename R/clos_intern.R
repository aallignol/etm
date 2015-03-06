### To be used for single endpoint
clos.nocp <- function(x, aw, ratio) {

    dims <- dim(x$est)
    los <- matrix(rep(x$time, 3), ncol = 3, byrow = FALSE)
    tau <- max(x$data$exit)
    times <- x$time
    surv <- x$est[1, 1, ]

    ## Call to C++ function
    zzz <- .Call("los_nocp",
                 times,
                 x$delta.na,
                 tau)
    zzz
}
