options(digits = 17, warn = 1)
suppressPackageStartupMessages({library(forecast); library(urca); library(parallel)})
idx <- rbind(read.csv("/work/index_A.csv"), read.csv("/work/index_B.csv"))
ser <- rbind(read.csv("/work/series_A.csv"), read.csv("/work/series_B.csv"))
kp  <- rbind(read.csv("/work/kpss_A.csv"), read.csv("/work/kpss_B.csv"))
f17 <- function(v) formatC(v, digits = 17, format = "g")
one <- function(i) {
  r <- idx[i, ]; id <- r$series_id
  v <- ser$value[ser$series_id == id]; v <- v[order(ser$t[ser$series_id == id])]
  x <- ts(v, frequency = r$period)
  rows <- character(0)
  fitrow <- function(mode, fit) {
    o <- arimaorder(fit); g <- function(k) if (k %in% names(o)) o[[k]] else 0
    cn <- names(coef(fit))
    sprintf("%s,%s,%d,%d,%d,%d,%d,%d,%d,%d,%s", id, mode, g("p"), g("d"), g("q"), g("P"), g("D"), g("Q"),
            as.integer(any(cn %in% c("intercept", "mean"))), as.integer("drift" %in% cn), f17(fit$aicc))
  }
  safe <- function(mode, expr) tryCatch(fitrow(mode, expr), error = function(e) sprintf("%s,%s,NA,NA,NA,NA,NA,NA,NA,NA,NA", id, mode))
  if (r$period > 1) {
    rows <- c(rows, safe("step", auto.arima(x, d = r$d_star, D = r$D_star, stepwise = TRUE, approximation = FALSE)),
                    safe("ex",   auto.arima(x, d = r$d_star, D = r$D_star, stepwise = FALSE, approximation = FALSE, max.order = 5)))
  } else {
    rows <- c(rows, safe("step", auto.arima(x, d = r$d_star, stepwise = TRUE, approximation = FALSE)),
                    safe("ex",   auto.arima(x, d = r$d_star, stepwise = FALSE, approximation = FALSE, max.order = 5)))
  }
  rows <- c(rows, safe("auto", auto.arima(x, test = "kpss", stepwise = TRUE, approximation = FALSE)))
  # KPSS steps with our lag
  z <- if (r$D_star > 0) diff(x, lag = r$period, differences = r$D_star) else x
  nd <- tryCatch(ndiffs(z, test = "kpss", type = "level", alpha = 0.05), error = function(e) NA)
  ks <- kp[kp$series_id == id, ]; ks <- ks[order(ks$step), ]
  krows <- character(0)
  for (j in seq_len(nrow(ks))) {
    zk <- if (ks$step[j] > 1) diff(z, differences = ks$step[j] - 1) else z
    eta <- tryCatch(ur.kpss(as.numeric(zk), type = "mu", use.lag = ks$lag[j])@teststat, error = function(e) NA)
    krows <- c(krows, sprintf("%s,%d,%d,%d,%s,%s", id, ks$step[j], length(zk), ks$lag[j], f17(eta), as.character(nd)))
  }
  srow <- NULL
  if (r$period > 1) {
    m <- tryCatch(mstl(x), error = function(e) NULL)
    fs <- if (is.null(m)) NA else { rem <- m[, "Remainder"]; sea <- m[, grep("^Seasonal", colnames(m))[1]]; max(0, 1 - var(rem) / var(sea + rem)) }
    srow <- sprintf("%s,%s,%d", id, f17(fs), tryCatch(nsdiffs(x), error = function(e) NA))
  }
  list(o = rows, k = krows, s = srow)
}
res <- mclapply(seq_len(nrow(idx)), one, mc.cores = 4)
w <- function(path, header, lines) writeLines(c(header, lines), path)
w("/work/r_orders.csv", "series_id,mode,p,d,q,P,D,Q,has_mean,has_drift,aicc_R", unlist(lapply(res, `[[`, "o")))
w("/work/r_kpss.csv", "series_id,step,T,lag,eta_R,ndiffs_R", unlist(lapply(res, `[[`, "k")))
w("/work/r_seasonal.csv", "series_id,F_S_STL,nsdiffs_R", unlist(lapply(res, function(z) z$s)))
writeLines(c(paste("R", R.version.string), paste("forecast", as.character(packageVersion("forecast"))), paste("urca", as.character(packageVersion("urca"))),
  "image rocker/r-ver:4.4.1; run.R in this directory; date", format(Sys.time(), "%Y-%m-%d %H:%M UTC", tz = "UTC"),
  "step: auto.arima(x, d=d_star, D=D_star, stepwise=TRUE, approximation=FALSE)",
  "ex:   auto.arima(x, d=d_star, D=D_star, stepwise=FALSE, approximation=FALSE, max.order=5)",
  "auto: auto.arima(x, test='kpss', stepwise=TRUE, approximation=FALSE)",
  "kpss: ur.kpss(z_k, type='mu', use.lag=<our lag>)@teststat; ndiffs(z, test='kpss', type='level', alpha=0.05)",
  "seasonal: mstl(x) F_S = max(0, 1 - var(rem)/var(seas+rem)); nsdiffs(x)"), "/work/PROVENANCE.txt")
cat("DONE", nrow(idx), "\n")
