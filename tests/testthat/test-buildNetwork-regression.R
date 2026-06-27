# Regression: the candidate-indexing rewrite must reproduce, as a SET, the
# exact edges + distances of the original brute-force implementation across the
# full parameter matrix. `expand = "clique"` is the bit-identical mode (star
# changes edge multiplicity by design and is checked separately).
#
# The golden snapshot in fixtures/buildNetwork_golden.rds was generated from the
# pre-rewrite implementation. Regenerate ONLY with intent via
# data-raw/make-buildNetwork-golden.R (never blindly).

test_that("buildNetwork reproduces frozen golden edge sets (clique mode)", {
  golden_path <- test_path("fixtures", "buildNetwork_golden.rds")
  skip_if_not(file.exists(golden_path), "golden snapshot missing")
  golden <- readRDS(golden_path)
  sc <- .bn_scenarios()

  # Pass expand = "clique" only once the argument exists (post-rewrite).
  extra <- if ("expand" %in% names(formals(buildNetwork))) list(expand = "clique") else list()

  fails <- character(0)
  for (i in seq_along(golden)) {
    g <- golden[[i]]
    scen <- sc$grid[[i]]
    got <- tryCatch(
      suppressMessages(.bn_run(sc$datasets, scen, extra)),
      error = function(e) paste0("ERROR:", conditionMessage(e)))
    # Compare as SETS. The contract is the edge set, not row order, so a
    # difference in collation order (e.g. C locale under R CMD check vs the
    # locale the golden was built in) must not register as a failure.
    only_new <- setdiff(got, g$set)
    only_old <- setdiff(g$set, got)
    if (length(only_new) || length(only_old)) {
      lab <- sprintf("[%d] %s {%s}", i, g$data,
                     paste(names(scen$args), unlist(scen$args), sep = "=", collapse = ", "))
      fails <- c(fails, sprintf("%s  +%d/-%d", lab, length(only_new), length(only_old)))
    }
  }
  expect_identical(fails, character(0),
                   info = paste0("\n", paste(utils::head(fails, 20), collapse = "\n")))
})
