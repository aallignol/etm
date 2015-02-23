### To be used for single endpoint
clos.nocp <- function(x, tr.mat, aw, ratio) {

    dims <- dim(x$est)
    los <- matrix(rep(x$time, 3), ncol = 3, byrow = FALSE)
    tau <- max(x$data$exit)

    
