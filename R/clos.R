### Expected change of LoS
### Arthur Allignol <arthur.allignol@uni-ulm.de>
#####################################################################


clos <- function(x, aw, ratio, ...) {
    UseMethod("clos")
}

clos.etm <- function(x, aw = FALSE, ratio = FALSE, ...) {
    if (!inherits(x, "etm")) {
        stop("'x' must be an 'etm' object")
    }
    if (is.null(x$delta.na)) {
        stop("Needs the increment of the Nelson-Aalen estimator")
    }

    ## test if we have an illness-death model
    ## if (!(dim(x$est)[1] %in% c(3, 4))) {
    ##     stop("The multistate model should be an illness-death model without recovery (possibly with competing endpoints")
    ## } else {
    ##     if (!(all(x$tra == tra_ill()) | all(x$tra == tra_ill_comp()))) {
    ##         stop("The multistate model should be an illness-death model without recovery (possibly with competing endpoints")
    ##     }
    ## }
    
    dims <- dim(x$est)
    comp.risk <- FALSE
    if (dims[1] == 4) comp.risk <- TRUE
    ## I <- diag(1, dims[1])
    ## tr.mat <- array(apply(x$delta.na, 3, "+", I), dim = dims)
    if (comp.risk) {
        ## res <- clos.cp(x, tr.mat, aw, ratio)
        stop("not yet")
    }
    else res <- clos.nocp(x, aw, ratio)
    
    class(res) <- "clos.etm"
    res
}
