options(digits = 17, warn = 1)
suppressPackageStartupMessages({library(forecast); library(urca)})
idx <- read.csv("/work/index.csv"); ser <- read.csv("/work/series.csv"); kp <- read.csv("/work/ours_kpss.csv")
f17 <- function(v) formatC(v, digits = 17, format = "g")
o_rows <- c(); k_rows <- c(); s_rows <- c()
fitrow <- function(id, mode, fit) { o <- arimaorder(fit); g <- function(k) if (k %in% names(o)) o[[k]] else 0; cn <- names(coef(fit))
  sprintf("%s,%s,%d,%d,%d,%d,%d,%d,%d,%d,%s", id, mode, g("p"), g("d"), g("q"), g("P"), g("D"), g("Q"),
          as.integer(any(cn %in% c("intercept","mean"))), as.integer("drift" %in% cn), f17(fit$aicc)) }
for (i in seq_len(nrow(idx))) { r <- idx[i,]; id <- r$series_id
  v <- ser$value[ser$series_id == id]; v <- v[order(ser$t[ser$series_id == id])]; x <- ts(v, frequency = r$period); seas <- r$period > 1
  safe <- function(mode, expr) tryCatch(fitrow(id, mode, expr), error = function(e) sprintf("%s,%s,NA,NA,NA,NA,NA,NA,NA,NA,NA", id, mode))
  o_rows <- c(o_rows,
    safe("step", auto.arima(x, stepwise = TRUE,  approximation = FALSE, test = "kpss", max.p = 5, max.q = 5, max.P = 2, max.Q = 2, ic = "aicc", seasonal = seas)),
    safe("ex",   auto.arima(x, stepwise = FALSE, approximation = FALSE, test = "kpss", max.p = 5, max.q = 5, max.P = 2, max.Q = 2, ic = "aicc", seasonal = seas, max.order = 5)))
  z <- if (r$D_star > 0) diff(x, lag = r$period, differences = r$D_star) else x
  nd <- tryCatch(ndiffs(z, test = "kpss", type = "level", alpha = 0.05), error = function(e) NA)
  ks <- kp[kp$series_id == id,]; ks <- ks[order(ks$step),]
  for (j in seq_len(nrow(ks))) { zk <- if (ks$step[j] > 1) diff(z, differences = ks$step[j] - 1) else z
    eta <- tryCatch(ur.kpss(as.numeric(zk), type = "mu", use.lag = ks$lag[j])@teststat, error = function(e) NA)
    k_rows <- c(k_rows, sprintf("%s,%d,%d,%d,%s,%s", id, ks$step[j], length(zk), ks$lag[j], f17(eta), as.character(nd))) }
  if (seas) { m <- tryCatch(mstl(x), error = function(e) NULL)
    fs <- if (is.null(m)) NA else { rem <- m[, "Remainder"]; sea <- m[, grep("^Seasonal", colnames(m))[1]]; max(0, 1 - var(rem)/var(sea + rem)) }
    s_rows <- c(s_rows, sprintf("%s,%s,%d", id, f17(fs), tryCatch(nsdiffs(x), error = function(e) NA))) }
}
writeLines(c("series_id,mode,p,d,q,P,D,Q,has_mean,has_drift,aicc_R", o_rows), "/work/r_orders.csv")
writeLines(c("series_id,step,T,lag,eta_R,ndiffs_R", k_rows), "/work/r_kpss.csv")
writeLines(c("series_id,F_S_STL,nsdiffs_R", s_rows), "/work/r_seasonal.csv")
writeLines(c(paste("R", R.version.string), paste("forecast", as.character(packageVersion("forecast"))), paste("urca", as.character(packageVersion("urca"))),
  "image rocker/r-ver:4.4.1; run.R (windows); date", format(Sys.time(), "%Y-%m-%d %H:%M UTC", tz = "UTC"),
  "step: auto.arima(x, stepwise=TRUE, approximation=FALSE, test='kpss', max.p=5, max.q=5, max.P=2, max.Q=2, ic='aicc', seasonal=(s>1))",
  "ex:   same with stepwise=FALSE, max.order=5",
  "kpss: ur.kpss(z_k, type='mu', use.lag=<our lag>)@teststat; ndiffs(z, test='kpss', type='level', alpha=0.05)",
  "seasonal: mstl(x) F_S; nsdiffs(x)"), "/work/PROVENANCE.txt")
cat("DONE", nrow(idx), "\n")
