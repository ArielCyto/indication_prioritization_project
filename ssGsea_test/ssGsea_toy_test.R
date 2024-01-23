library(cytoreason.gx)
library(Matrix, lib.loc = "/usr/lib/R/library")
library(checkmate)

# toy data (20 samples)
x <- matrix(runif(20 * 1000), ncol = 20,
            dimnames = list(paste0("g", 1:1000), NULL))
# gene-set collection
gs <- setNames(
  replicate(10, paste0("g", sample.int(1000, 5)), simplify = FALSE),
  LETTERS[1:10]
)

# fit
res <- service_ssgsea(x, gs)
