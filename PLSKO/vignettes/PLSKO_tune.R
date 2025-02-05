## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = FALSE,
  comment = "#>"
)

## ----case1tune----------------------------------------------------------------
library(PLSKO)
data(cfRNA_placenta)
str(cfRNA_placenta)

# Basic input
X <- cfRNA_placenta$counts
n_ko <- 10
q <- 0.05
p_s <- 10

# Parameters to test
ncomp <- c(3, 5, 7)
threshold.abs <- c(0, 0.3, 0.5, 0.8)
sparsity <- 1
# We test on:
expand.grid(ncomp = ncomp, threshold.abs = threshold.abs, sparsity = sparsity)

# Run the tuning
start.time <- Sys.time()

cfRNA_tune <- plsko_tuning(X, n_ko = n_ko, p_s = p_s, q = q, ncomp = ncomp, threshold.abs = threshold.abs, sparsity = sparsity)

end.time <- Sys.time()

# Print the results
plot(cfRNA_tune)

# Print the time
end.time - start.time

## ----case1real----------------------------------------------------------------
y <- cfRNA_placenta$metadata$PE
ncomp <- cfRNA_tune$optimal$ncomp
threshold.abs <- cfRNA_tune$optimal$threshold.abs
sparsity <- cfRNA_tune$optimal$sparsity

cfRNA_result <- plsAKO(X, y, n_ko = 25, ncomp = ncomp, threshold.abs = threshold.abs, sparsity = sparsity)
print(cfRNA_result)

# Or if you like using the configuration with the lowest average FDP
user_optimal <- cfRNA_tune$mean[which.min(cfRNA_tune$mean$mean.fdp), ]
cfRNA_result_alt <- plsAKO(X, y, n_ko = 25, ncomp = user_optimal$ncomp, threshold.abs = user_optimal$threshold.abs, sparsity = user_optimal$sparsity)
print(cfRNA_result_alt)

## -----------------------------------------------------------------------------
sessionInfo()

