options(digits = 10, warn = 1, width = 200)
suppressPackageStartupMessages(library(forecast))
idx <- rbind(read.csv("/work/index_A.csv"), read.csv("/work/index_B.csv"))
ser <- rbind(read.csv("/work/series_A.csv"), read.csv("/work/series_B.csv"))
widx <- read.csv("/work/win_index.csv"); wser <- read.csv("/work/win_series.csv")
getx <- function(id) {
  if (id %in% idx$series_id) { r <- idx[idx$series_id == id, ]; s <- ser } else { r <- widx[widx$series_id == id, ]; s <- wser }
  v <- s$value[s$series_id == id]; v <- v[order(s$t[s$series_id == id])]
  list(x = ts(v, frequency = r$period), r = r)
}
mr <- function(co) if (length(co) == 0) NA else min(Mod(polyroot(c(1, co))))
cat("# R", R.version.string, "| forecast", as.character(packageVersion("forecast")), "\n")
cat("\n#### PART 1: refits of the R orders rejected by our root check (and ima_n200_r55 start model)\n")
cat("series,order,const,method,loglik,aicc,min_ar,min_ma,coef,se\n")
cases <- list(
  list("ar2_n200_r24", c(2,0,2), "mean"), list("ar098_n200_r59", c(1,1,1), "drift"),
  list("ar098_n200_r65", c(1,1,2), "none"), list("ar098_n500_r77", c(2,1,2), "none"),
  list("ari_n500_r89", c(2,1,1), "drift"), list("wwwusage_o6", c(2,0,3), "mean"),
  list("wwwusage_o8", c(2,0,3), "mean"), list("wwwusage_o9", c(2,0,3), "mean"),
  list("ima_n200_r55", c(2,1,2), "drift"), list("ima_n200_r55", c(3,1,0), "none"))
for (cs in cases) for (m in c("CSS-ML", "ML")) {
  g <- getx(cs[[1]]); o <- cs[[2]]; cst <- cs[[3]]
  fit <- tryCatch(Arima(g$x, order = o, include.mean = (cst == "mean"), include.drift = (cst == "drift"), method = m), error = function(e) NULL)
  if (is.null(fit)) { cat(sprintf("%s,(%d%d%d),%s,%s,ERROR,,,,,\n", cs[[1]], o[1],o[2],o[3], cst, m)); next }
  cf <- coef(fit); se <- sqrt(pmax(diag(fit$var.coef), 0)); negv <- any(diag(fit$var.coef) < 0)
  cat(sprintf("%s,(%d%d%d),%s,%s,%.6f,%.6f,%s,%s,%s,%s%s\n", cs[[1]], o[1],o[2],o[3], cst, m, fit$loglik, fit$aicc,
      format(mr(-cf[grep("^ar[0-9]", names(cf))])), format(mr(cf[grep("^ma[0-9]", names(cf))])),
      paste(names(cf), format(cf, digits = 8), sep = "=", collapse = ";"), paste(format(se, digits = 6), collapse = ";"), if (negv) ";NEGATIVE_VAR" else ""))
}
cat("\n#### PART 2: seasonal stepwise traces (airline, our d*, D*)\n")
ids <- c(paste0("airline_n200_r", 0:14), paste0("airline_n500_r", 0:14))
for (id in ids) { g <- getx(id); r <- g$r
  cat(sprintf("\n## %s  n=%d period=%d d*=%d D*=%d\n", id, length(g$x), r$period, r$d_star, r$D_star))
  f <- auto.arima(g$x, d = r$d_star, D = r$D_star, stepwise = TRUE, approximation = FALSE, trace = TRUE)
  cat("Chosen:", as.character(f), " aicc=", f$aicc, "\n") }
cat("\n#### PART 2b: generated (1,0,0)(0,1,1)[12] with nonzero drift in the seasonal difference, d=0, D=1\n")
for (sd in 1:5) { set.seed(1000 + sd)
  e <- arima.sim(list(ar = 0.5, ma = c(rep(0, 11), -0.6)), n = 228) + 0.8
  y <- ts(diffinv(e, lag = 12, xi = 50 + 5 * sin(2 * pi * (1:12) / 12)), frequency = 12)
  cat(sprintf("\n## gen_sar_d0D1_s%d  n=%d period=12 d*=0 D*=1\n", sd, length(y)))
  f <- auto.arima(y, d = 0, D = 1, stepwise = TRUE, approximation = FALSE, trace = TRUE)
  cat("Chosen:", as.character(f), " aicc=", f$aicc, "\n") }
cat("\n#### PART 3: traces for attribution\n")
for (id in c("wn_n500_r27", "arma11_n500_r14", "arma11_n500_r98", "ar098_n500_r95")) { g <- getx(id); r <- g$r
  cat(sprintf("\n## %s  n=%d period=%d d*=%d D*=%d\n", id, length(g$x), r$period, r$d_star, r$D_star))
  f <- auto.arima(g$x, d = r$d_star, stepwise = TRUE, approximation = FALSE, trace = TRUE)
  cat("Chosen:", as.character(f), " aicc=", f$aicc, "\n") }
