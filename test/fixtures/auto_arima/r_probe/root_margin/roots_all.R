options(digits = 8, warn = 1, width = 220)
suppressPackageStartupMessages(library(forecast))
idx <- rbind(read.csv("/work/index_A.csv"), read.csv("/work/index_B.csv"))
ser <- rbind(read.csv("/work/series_A.csv"), read.csv("/work/series_B.csv"))
tr <- readLines("/work/traces.txt")
cur <- NULL
minmod <- function(coefs) { co <- coefs[abs(coefs) > 1e-8 | seq_along(coefs) <= length(coefs)]; if (length(coefs) == 0) return(NA); rt <- polyroot(c(1, coefs)); if (length(rt) == 0) NA else min(Mod(rt)) }
cat("# For every model printed as Inf in the stepwise trace: refit with Arima (method CSS-ML like auto.arima) and report min |root| of AR/MA (and seasonal) polynomials\n")
cat("series,model,const,r_status,aicc_refit,min_ar,min_ma,min_sar,min_sma,converged\n")
for (ln in tr) {
  m <- regmatches(ln, regexec("^## ([^ ]+)", ln))[[1]]; if (length(m) == 2) { cur <- m[2]; next }
  if (!grepl(" : ", ln) || is.null(cur)) next; acc <- if (grepl(": Inf", ln)) "INF" else "ACC"
  g <- regmatches(ln, regexec("ARIMA[(]([0-9]),([0-9]),([0-9])[)]([(]([0-9]),([0-9]),([0-9])[)][[]([0-9]+)[]])?[ ]*(with drift|with non-zero mean|with zero mean)?", ln))[[1]]
  if (length(g) < 4) next
  r <- idx[idx$series_id == cur, ]; v <- ser$value[ser$series_id == cur]; v <- v[order(ser$t[ser$series_id == cur])]; x <- ts(v, frequency = r$period)
  o <- as.integer(g[2:4]); so <- if (nzchar(g[5])) as.integer(g[6:8]) else c(0,0,0); cst <- g[10]
  fit <- tryCatch(Arima(x, order = o, seasonal = so, include.mean = (cst == "with non-zero mean"), include.drift = (cst == "with drift"), method = "CSS-ML"), error = function(e) NULL)
  if (is.null(fit)) { cat(sprintf("%s,(%d%d%d)(%d%d%d),%s,ERROR,,,,,\n", cur, o[1],o[2],o[3],so[1],so[2],so[3], cst, acc)); next }
  cf <- coef(fit); ar <- -cf[grep("^ar[0-9]", names(cf))]; ma <- cf[grep("^ma[0-9]", names(cf))]; sar <- -cf[grep("^sar[0-9]", names(cf))]; sma <- cf[grep("^sma[0-9]", names(cf))]
  mr <- function(co) if (length(co) == 0) NA else min(Mod(polyroot(c(1, co))))
  cat(sprintf("%s,(%d%d%d)(%d%d%d),%s,%s,%.4f,%s,%s,%s,%s,%d
", cur, o[1],o[2],o[3],so[1],so[2],so[3], cst, acc, fit$aicc, format(mr(ar)), format(mr(ma)), format(mr(sar)), format(mr(sma)), fit$code))
}
