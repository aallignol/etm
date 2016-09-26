pseudo_values <- function(x,
                         statistic = c("elos",
                                       "los",
                                       "state_occupation"),
                         ...) {
    UseMethod("pseudoValues")
}

pseudo_values.data.frame <- function(x,
                                     statistic = c("elos",
                                                   "los",
                                                   "state_occupation"),
                                     state.names,
                                     tra,
                                     cens.name,
                                     s = 0,
                                     t, ...)
{
                                     
    ## All the data verifications
    ##
    ##

    x <- data.table(x)

    reg <- names(x)
    names_msm <- intersect(c("id", "entry", "exit", "time", "from", "to"), reg)
    x <- x[, names_msm, with = FALSE]

    colnames(tra) <- rownames(tra) <- state.names
    t.from <- lapply(1:dim(tra)[2], function(i) {
        rep(rownames(tra)[i], sum(tra[i, ]))
    })
    t.from <- unlist(t.from)
    t.to <- lapply(1:dim(tra)[2], function(i) {
        colnames(tra)[tra[i, ]==TRUE]
    })
    t.to <- unlist(t.to)
    trans <- data.frame(from=t.from, to=t.to)
    namen <- paste(trans[, 1], trans[, 2])
    
    ## test on transitions
    test <- x[, unique(paste(from, to))]
    if (!(is.null(cens.name))) {
        ref <- c(paste(trans$from, trans$to), paste(unique(trans$from), cens.name))
    } else {
        ref <- paste(trans$from, trans$to)
    }
    ref.wo.cens <- paste(trans$from, trans$to)
    if (!(all(test %in% ref)==TRUE))
        stop("There is undefined transitions in the data set")
    if (x[, sum(as.character(from) == as.character(to))] > 0)
        stop("Transitions into the same state are not allowed")
    if (!(all(ref.wo.cens %in% test) == TRUE))
        warning("You may have specified more possible transitions than actually present in the data")

    n <- length(unique(x$id))

    ## data.table transformation
    x[, id := if (is.character(id)) as.factor(id) else id]
    if (!(is.null(cens.name))) {
        x[, from := factor(from, levels = c(cens.name, state.names), ordered = TRUE)]
        levels(x$from) <- 0:length(state.names)
        x[, to := factor(to, levels = c(cens.name, state.names), ordered = TRUE)]
        levels(x$to) <- 0:length(state.names)
    } else{
        x[, from := factor(from, levels = state.names, ordered = TRUE)]
        levels(x$from) <- 1:length(state.names)
        x[, to := factor(x$to, levels = state.names, ordered = TRUE)]
        levels(x$to) <- 1:length(state.names)
    }

### if not, put like counting process data
    if ("time" %in% names(x)) {
        x <- transfo_to_counting(x)
        if (sum(x$entry < x$exit) != nrow(x))
            stop("Exit time from a state must be > entry time")
    } else {
        if (sum(x$entry < x$exit) != nrow(x))
            stop("Exit time from a state must be > entry time")
    }

    x[, from := as.integer(as.character(from))]
    x[, to := as.integer(as.character(to))]
    
    if (t=="last") t <- max(x$exit)
    if (!(0 <= s & s < t))
        stop("'s' and 't' must be positive, and s < t")
    if (t <= x[, min(exit)] | s >= x[, max(exit)])
        stop("'s' or 't' is an invalid time")

    ## remove the lines in which transition before s
    x <- x[exit > s]

    ## call to etm, get the risks set and n.events plus the statistics of interest
    zzz <- .etm(entry = x[, entry],
                exit = x[, exit],
                from = x[, from],
                to = x[, to],
                nstate = dim(tra)[1],
                s,
                t,
                covariance = FALSE,
                c_modif = 0)

    ## Find at which times (index of time) the ind are at risk
    span_entry <- x[, findInterval(entry, zzz$time, all.inside = FALSE) + 1]
    span_exit <- x[, findInterval(exit, zzz$time, all.inside = TRUE)]

    ## Go to C++

}
