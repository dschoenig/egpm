chunk <- seq(1, 10, 2)

files.res <- paste0("./test_",
                    stri_pad_left(chunk, 4, 0),
                    ".rds")

results <- list()
for(i in chunk) {
  results[[i]] <- rnorm(10)
}

for(i in seq_along(results)) {
  saveRDS(object = results[[i]], file = file.res[i])
  }
 
