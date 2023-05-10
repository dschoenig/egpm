mod.names <-
  list(c("binary_imbalance_high", "egp_som25"),
       c("binary_imbalance_high", "match"),
       c("binary_imbalance_low", "egp_som25"),
       c("binary_imbalance_low", "match"),
       c("imbalance_high", "egp_som10"),
       c("imbalance_high", "egp_som25"),
       c("imbalance_high", "egp_som25_noaprox"),
       c("imbalance_high", "egp_som50"),
       c("imbalance_high", "egp_som25_sam005"),
       c("imbalance_high", "egp_som25_sam02"),
       c("imbalance_high", "match"),
       c("imbalance_high", "match_sam005"),
       c("imbalance_high", "match_sam02"),
       c("imbalance_low", "egp_som25"),
       c("imbalance_low", "match"))

sim.res <- logical(0)
for(i in seq_along(mod.names)) {
  path.log <- paste0("../models/", mod.names[[i]][1], "/", mod.names[[i]][2], ".log")
  sim.log <- read.table(path.log)[,1]
  sim.res[i] <- all(1:1000 %in% sim.log)
}

all(sim.res)
sim.res
