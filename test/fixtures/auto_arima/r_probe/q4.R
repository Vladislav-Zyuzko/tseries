options(digits = 15, warn = 1, width = 220)
suppressPackageStartupMessages({library(forecast); library(urca)})
crit <- c("0.01" = 0.739, "0.025" = 0.574, "0.05" = 0.463, "0.1" = 0.347)
al <- c(0.01, 0.025, 0.05, 0.10)
probe <- function(nm, x) {
  x <- as.numeric(x); n <- length(x)
  nd <- sapply(al, function(a) ndiffs(x, test = "kpss", alpha = a))
  st1 <- sapply(0:6, function(l) ur.kpss(x, type = "mu", use.lag = l)@teststat)
  u <- 2:4; ok1 <- sapply(0:6, function(l) all((st1[l + 1] > crit[u]) == (nd[u] >= 1)))
  line <- sprintf("%-22s n=%3d ndiffs[a=.01,.025,.05,.10]=%s | step1 lags consistent: {%s} | rule3sqrt/13=%d sqrt/4=%d l4=%d", nm, n,
    paste(nd, collapse = ","), paste((0:6)[ok1], collapse = ","), floor(3 * sqrt(n) / 13), floor(sqrt(n) / 4), floor(4 * (n / 100)^0.25))
  if (any(nd >= 1)) { y <- diff(x); m <- length(y)
    st2 <- sapply(0:6, function(l) ur.kpss(y, type = "mu", use.lag = l)@teststat)
    sel <- (nd >= 1) & c(FALSE, TRUE, TRUE, TRUE)
    ok2 <- sapply(0:6, function(l) all((st2[l + 1] > crit[sel]) == (nd[sel] >= 2)))
    line <- paste0(line, sprintf(" || step2 (T=%d) consistent: {%s} rule=%d", m, paste((0:6)[ok2], collapse = ","), floor(3 * sqrt(m) / 13)))
  }
  cat(line, "\n")
}
set.seed(20261006)
la <- log(AirPassengers)
for (n in c(70, 75, 76, 150, 168, 170)) {
  cat("## n =", n, "\n")
  for (k in 1:6) { rw <- cumsum(rnorm(n)); probe(sprintf("rw_%d", k), rw) }
  for (k in 1:4) { ar <- as.numeric(arima.sim(list(ar = 0.9), n)); probe(sprintf("ar0.9_%d", k), ar) }
  if (n <= 98) probe("lake_huron_head", LakeHuron[1:n]); if (n <= 100) { probe("nile_head", Nile[1:n]); probe("wwwusage_head", WWWusage[1:n]) }
  if (n <= 144) probe("air_log_head", la[1:n])
}
cat("## n = 100, three DGPs\n")
for (k in 1:5) { probe(sprintf("wn_%d", k), rnorm(100)); probe(sprintf("ar0.5_%d", k), arima.sim(list(ar = 0.5), 100)); probe(sprintf("ar0.95_%d", k), arima.sim(list(ar = 0.95), 100)) }
cat("## length 76, twice integrated (looking for ndiffs = 2)\n")
for (k in 1:8) probe(sprintf("i2_%d", k), cumsum(cumsum(rnorm(76))))
