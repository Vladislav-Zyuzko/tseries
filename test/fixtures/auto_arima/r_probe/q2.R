options(digits = 15, warn = 1, width = 200)
suppressPackageStartupMessages({library(forecast); library(urca)})
la <- log(AirPassengers)
sets <- list(air_log = la, air_log_sdiff12 = diff(la, 12), wwwusage = WWWusage, wwwusage_diff1 = diff(WWWusage), lake_huron = LakeHuron,
             nile = Nile, lynx = lynx, us_acc_deaths = USAccDeaths, us_acc_sdiff12 = diff(USAccDeaths, 12))
for (nm in names(sets)) { x <- sets[[nm]]
  st <- sapply(0:8, function(l) ur.kpss(x, type = "mu", use.lag = l)@teststat)
  cat(sprintf("%s n=%d ndiffs05=%d | teststat by use.lag 0..8: %s\n", nm, length(x), ndiffs(x, test = "kpss", alpha = 0.05), paste(sprintf("%.6f", st), collapse = " ")))
}
cat("# ndiffs at alpha grid (wwwusage / air_log_sdiff12):\n")
for (a in c(0.01, 0.025, 0.05, 0.1)) cat(sprintf("alpha=%.3f wwwusage=%d air_sdiff12=%d\n", a, ndiffs(WWWusage, alpha = a), ndiffs(diff(la, 12), alpha = a)))
