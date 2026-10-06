# Produces test/fixtures/auto_arima/reference_r_forecast.txt and the series CSVs.
# Run outside the package build, e.g.:
#   docker run --rm -v "$PWD:/work" -w /work rocker/r-ver:4.4.1 \n#     bash -c "R -q -e 'install.packages(\"forecast\")' && Rscript auto_arima_reference.R"
# R and forecast (GPL) are used only as an external tool to obtain numbers.
options(digits = 15, warn = 1)
suppressPackageStartupMessages(library(forecast))
sets <- list(
  air_passengers_log = log(AirPassengers), wwwusage = WWWusage, lynx = lynx,
  lake_huron = LakeHuron, nile = Nile, us_acc_deaths = USAccDeaths)
out <- file("reference.txt", "w")
w <- function(...) writeLines(paste0(...), out)
w("# R ", R.version.string, " | forecast ", as.character(packageVersion("forecast")), " | ", format(Sys.time(), "%Y-%m-%d", tz = "UTC"))
w("# modes: stepwise = auto.arima(y, stepwise=TRUE, approximation=FALSE); exhaustive = auto.arima(y, stepwise=FALSE, approximation=FALSE, max.order=5)")
w("# name|mode|p|d|q|P|D|Q|period|constant|loglik|aic|aicc|bic|sigma2|coef")
for (nm in names(sets)) {
  y <- sets[[nm]]
  write.table(data.frame(value = as.numeric(y)), file = paste0("series_", nm, ".csv"),
              row.names = FALSE, col.names = FALSE)
  w("# ", nm, ": n=", length(y), " frequency=", frequency(y),
    " ndiffs_kpss=", ndiffs(y, test = "kpss"),
    if (frequency(y) > 1) paste0(" nsdiffs_seas=", nsdiffs(y, test = "seas")) else "")
  for (mode in c("stepwise", "exhaustive")) {
    fit <- if (mode == "stepwise") auto.arima(y, stepwise = TRUE, approximation = FALSE) else auto.arima(y, stepwise = FALSE, approximation = FALSE, max.order = 5)
    o <- arimaorder(fit); g <- function(k) if (k %in% names(o)) o[[k]] else 0
    cn <- names(coef(fit)); const <- if ("drift" %in% cn) "drift" else if (any(cn %in% c("intercept", "mean"))) "mean" else "none"
    w(nm, "|", mode, "|", g("p"), "|", g("d"), "|", g("q"), "|", g("P"), "|", g("D"), "|", g("Q"), "|",
      if ("Frequency" %in% names(o)) o[["Frequency"]] else 1, "|", const, "|",
      format(as.numeric(logLik(fit)), digits = 15), "|", format(fit$aic, digits = 15), "|",
      format(fit$aicc, digits = 15), "|", format(fit$bic, digits = 15), "|", format(fit$sigma2, digits = 15), "|",
      paste(cn, format(as.numeric(coef(fit)), digits = 12), sep = "=", collapse = ","))
  }
}
close(out)
cat(readLines("reference.txt"), sep = "\n")
