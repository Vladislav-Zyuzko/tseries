options(digits = 15, warn = 1, width = 200)
suppressPackageStartupMessages({library(forecast); library(urca)})
cat("# R", R.version.string, "| forecast", as.character(packageVersion("forecast")), "| urca", as.character(packageVersion("urca")), "\n")
k <- function(nm, x) {
  u <- ur.kpss(x, type = "mu", lags = "short")
  cat(sprintf("KPSS %s: n=%d ur.kpss(mu,short) teststat=%.12f lag=%d cval10=%.3f cval5=%.3f cval2.5=%.3f cval1=%.3f | ndiffs(kpss,alpha=0.05)=%d ndiffs(alpha=0.10)=%d\n",
    nm, length(x), u@teststat, u@lag, u@cval[1], u@cval[2], u@cval[3], u@cval[4], ndiffs(x, test = "kpss", alpha = 0.05), ndiffs(x, test = "kpss", alpha = 0.10)))
}
la <- log(AirPassengers)
k("air_log", la); k("air_log_sdiff12", diff(la, 12)); k("air_log_diff1", diff(la)); k("wwwusage", WWWusage); k("wwwusage_diff1", diff(WWWusage))
k("lake_huron", LakeHuron); k("nile", Nile); k("lynx", lynx); k("us_acc_deaths", USAccDeaths); k("us_acc_sdiff12", diff(USAccDeaths, 12))
cat("# formals(ndiffs):", paste(names(formals(ndiffs)), sapply(formals(ndiffs), function(v) paste(deparse(v), collapse="")), sep="=", collapse="; "), "\n")
cat("# auto.arima air_log d chosen:", arimaorder(auto.arima(la, approximation = FALSE))[["d"]], " nsdiffs:", nsdiffs(la), "\n")
cat("\n## LakeHuron trace (d=1, stepwise, approximation=FALSE)\n")
f <- auto.arima(LakeHuron, d = 1, stepwise = TRUE, approximation = FALSE, trace = TRUE)
cat("\n## LakeHuron fixed orders\n")
for (o in list(c(0,1,0), c(0,1,1), c(1,1,0), c(1,1,2), c(2,1,1))) for (dr in c(FALSE, TRUE)) {
  m <- tryCatch(Arima(LakeHuron, order = o, include.drift = dr, method = "ML"), error = function(e) NULL)
  if (!is.null(m)) cat(sprintf("Arima(%d,%d,%d) drift=%s loglik=%.9f aicc=%.9f sigma2=%.9f coef=%s\n", o[1], o[2], o[3], dr, as.numeric(logLik(m)), m$aicc, m$sigma2, paste(names(coef(m)), sprintf("%.9f", coef(m)), sep="=", collapse=",")))
}
