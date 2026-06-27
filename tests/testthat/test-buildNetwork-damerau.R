# Guards the SymSpell radius trap: a pure adjacent transposition is Damerau
# distance 1 but Levenshtein distance 2, so the deletion index for Damerau must
# use radius 2*k. A radius-k index would silently miss these pairs.

test_that("Damerau finds adjacent transpositions that Levenshtein does not", {
  d <- data.frame(
    sequence_id = c("a", "b", "c", "d"),
    junction_aa = c("ABCDEFGHIK", "BACDEFGHIK",   # swap 1-2
                    "MNPQRSTVWY", "MNPQRSTVYW"),   # swap 9-10
    v_call = "V1", j_call = "J1", stringsAsFactors = FALSE)

  e_dam <- suppressMessages(buildNetwork(d, seq_col = "junction_aa", ids = d$sequence_id,
                                         dist_type = "damerau", threshold = 1))
  e_lev <- suppressMessages(buildNetwork(d, seq_col = "junction_aa", ids = d$sequence_id,
                                         dist_type = "levenshtein", threshold = 1))

  key <- function(df) if (nrow(df)) sort(paste(pmin(df$from, df$to), pmax(df$from, df$to))) else character(0)
  expect_identical(key(e_dam), c("a b", "c d"))
  expect_true(all(e_dam$dist == 1))
  expect_identical(key(e_lev), character(0))   # swaps are Levenshtein distance 2
})
