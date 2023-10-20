



library(data.table)
library(gt)

source("utilities.R")

path.base <- "../"
path.results <- paste0(path.base, "results/")
path.comp <- paste0(path.results, "comparison/")
file.egp.match.boot <- paste0(path.comp, "egp.match.boot.rds")

egp.match.boot <- readRDS(file.egp.match.boot)


# DATA PREPARATION ##################################################



egp.match.boot[ls.response == "all" & ls.imbalance != "all"] |>
format_comparisons(by.response = FALSE)

egp.match.boot[ls.response != "all" & ls.imbalance != "all"] |>
format_comparisons(by.response = TRUE)

