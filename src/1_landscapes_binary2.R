args <- commandArgs(trailingOnly = TRUE)

# SETUP #############################################################

source("utilities.R")

overwrite <- TRUE

ls.original <- args[1]
ls.binary <- args[2]
task.id <- as.integer(args[3])
task.count <- as.integer(args[4])
parallel <- as.integer(args[5])

# ls.original <- "imbalance_low"
# ls.binary <- "binary_imbalance_low"
# task.id <- 1
# task.count <- 1

paste0("Settings:\n",
       "  Original landscape ", ls.original, "\n",
       "  Binary landscape ", ls.binary) |>
message()

path.base <- "../"

path.ls.o <- paste0(path.base, "landscapes/", ls.original, "/")
path.ls.o.data <- paste0(path.ls.o, "data/")

path.ls.b <- paste0(path.base, "landscapes/", ls.binary, "/")
path.ls.b.data <- paste0(path.ls.b, "data/")
if(!dir.exists(path.ls.b.data)) dir.create(path.ls.b.data, recursive = TRUE)

file.par.o <- paste0(path.ls.o, "parameters.rds")
file.par.b <- paste0(path.ls.b, "parameters.rds")

file.log <- paste0(path.ls.b, "simulate.log")


# PARAMETERS ########################################################

par.o <- readRDS(file.par.o)
n <- nrow(par.o)

set.seed(18820125) # Virginia Woolfe
int.code <- as.integer(rawToBits(charToRaw(paste0(ls.original, ls.binary))))
ls.seeds <- round(runif(n, 0, sum(int.code)/length(int.code)) * 1e8)
par.b <-
  data.table(id = par.o$id,
             seed = ls.seeds,
             ls.original = ls.original,
             p.ref = 0.3,
             p.trt = 0.5,
             e.df = 1,
             opt.prec = sqrt(.Machine$double.eps),
             opt.grid = 1024
             )

par.b[, file.name := stri_pad_left(id, ceiling(log10(n))+1, 0)]
par.b[, file.path := paste0(path.ls.b.data, file.name, ".rds")]

if(!file.exists(file.par.b) | task.id == 1) saveRDS(par.b, file.par.b)

if(!overwrite & file.exists(file.log)) {
  simulated <- as.integer(read.table(file.log)[,1])
} else {
  simulated <- integer(0)
}

row.chunks <- chunk_seq(1,
                        nrow(par.b),
                        ceiling(nrow(par.b) / task.count))
chunk <- row.chunks$from[task.id]:row.chunks$to[task.id]

chunk <- setdiff(chunk, simulated)

# low <- read.table("sum_bin_low.txt")[,1]
# chunk <- which(low > 1e-4)

for(i in chunk) {

  message(paste0("Generating binary landscape ", i, " / ", nrow(par.b), " …"))

  ta <- Sys.time()

  ls.o <- readRDS(par.o$file.path[i])

  ls.par <- 
    as.list(par.b[id == i,]) |>
    lapply(unlist, recursive = FALSE)
  ls.par$ls <- ls.o$landscape
  ls.par$verbose <- 2
  ls.par$e.var <- par.o[id == i, sum(e.exp.var, e.nug.var)]
  file.ls.b <- ls.par$file.path
  ls.b <- NULL
  while(is.null(ls.b)) {
    try({ls.b <- do.call(generate_landscape_4cov_nl_binary, ls.par)})
    if(is.null(ls)) message("Simulation failed. Trying again …")
  }

  ls.b <-
    c(ls.b,
      ls.o[c("fun", "optim")])

  message(paste0("Saving to ", file.ls.b, " …"))
  saveRDS(ls.b, file.ls.b)

  rm(ls.b, ls.o)
  gc()

  tb <- Sys.time()
  te <- tb-ta
  print(te)

  system(paste0('echo "', i, '" >> ', file.log), intern = TRUE)

}

