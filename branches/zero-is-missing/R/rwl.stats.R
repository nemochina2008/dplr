`rwl.stats` <-
    function(rwl, na.rm.or.pass = TRUE, round.decimals = 3)
{
    check.flags(na.rm.or.pass)
    acf1 <- function(x, pass.na) {
        if (pass.na) {
            na.act <- na.pass
        } else {
            na.act <- na.fail
        }
        tryCatch(acf(x, na.action=na.act, lag.max=1, plot=FALSE)$acf[2],
                 error = function(...) {
                     NA_real_
                 })
    }
    skew <- function(x, na.rm) {
        if (na.rm) {
            y <- x[!is.na(x)]
        } else {
            y <- x
        }
        sum((y-mean(y))^3) / (length(y)*sd(y)^3)
    }

    rwl2 <- as.matrix(rwl)
    rwl.list <- apply(rwl2, 2,
                      function(x) {
                          idx.good <- which(!is.na(x))
                          n.good <- length(idx.good)
                          if (n.good > 0) {
                              x[idx.good[1]:idx.good[n.good]]
                          } else {
                              numeric(0)
                          }
                      })
    if (is.matrix(rwl.list)) {
        rwl.list <- as.data.frame(rwl.list)
    }

    series.stats <- data.frame(series=names(rwl))
    the.range <- vapply(rwl.list,
                        function(x) {
                            n.x <- length(x)
                            if (n.x == 0) {
                                c(NA_real_, NA_real_)
                            } else {
                                as.numeric(names(x)[c(1, n.x)])
                            }
                        },
                        numeric(2))
    series.stats$first <- the.range[1, ]
    series.stats$last <- the.range[2, ]
    series.stats$year <- series.stats$last - series.stats$first + 1
    series.stats$mean <- vapply(rwl.list, mean, 0, na.rm=na.rm.or.pass)
    series.stats$median <- vapply(rwl.list, median, 0, na.rm=na.rm.or.pass)
    series.stats$stdev <- vapply(rwl.list, sd, 0, na.rm=na.rm.or.pass)
    series.stats$skew <- vapply(rwl.list, skew, 0, na.rm=na.rm.or.pass)
    series.stats$sens1 <- vapply(rwl.list, sens1, 0, na.rm=na.rm.or.pass)
    series.stats$sens2 <- vapply(rwl.list, sens2, 0, na.rm=na.rm.or.pass)
    series.stats$gini <- vapply(rwl.list, gini.coef, 0, na.rm=na.rm.or.pass)
    series.stats$ar1 <- vapply(rwl.list, acf1, 0, pass.na=na.rm.or.pass)
    if (is.numeric(round.decimals) && length(round.decimals) > 0 &&
        is.finite(round.decimals[1]) && round.decimals[1] >= 0) {
        seq.temp <- -seq_len(4)
        series.stats[, seq.temp] <- round(series.stats[, seq.temp],
                                          round.decimals)
    }

    series.stats
}
