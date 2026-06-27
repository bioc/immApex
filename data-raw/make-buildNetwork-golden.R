# Regenerate the buildNetwork regression golden snapshot.
# Run deliberately (NOT as part of the rewrite) only when an intended change to
# buildNetwork's output set is being accepted. From the package root:
#   Rscript data-raw/make-buildNetwork-golden.R
suppressMessages(devtools::load_all(".", quiet = TRUE))
source("tests/testthat/helper-buildNetwork-regression.R")
sc <- .bn_scenarios()
extra <- if ("expand" %in% names(formals(buildNetwork))) list(expand = "clique") else list()
golden <- lapply(sc$grid, function(scen) {
  set <- suppressMessages(.bn_run(sc$datasets, scen, extra))
  list(data = scen$data, args = scen$args, set = set)
})
saveRDS(golden, "tests/testthat/fixtures/buildNetwork_golden.rds")
cat("Wrote", length(golden), "scenarios.\n")
