# Direct correctness check for the banded Levenshtein engine, independent of the
# golden snapshot. The banded routine previously missed pairs whose only optimal
# alignment rides the band edge (true distance == |length difference| == thr),
# including identical sequences at a tight normalized threshold. This compares
# buildNetwork's edge set against brute-force utils::adist for several
# thresholds and a length-diverse sequence set (which guarantees band edges are
# exercised).

test_that("Levenshtein edges match brute-force adist (incl. band edges)", {
  set.seed(2024)
  aa <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
  # wide length spread so |len diff| frequently equals the threshold
  seqs <- vapply(1:120, function(i) paste(sample(aa, sample(8:22, 1), replace = TRUE), collapse = ""),
                 character(1))
  ids <- paste0("s", seq_along(seqs))
  D <- as.matrix(utils::adist(seqs))

  for (thr in 0:4) {
    # brute-force truth: all i<j with edit distance <= thr
    truth <- which(D <= thr & upper.tri(D), arr.ind = TRUE)
    truth_key <- if (nrow(truth)) sort(paste(ids[pmin(truth[,1], truth[,2])],
                                              ids[pmax(truth[,1], truth[,2])])) else character(0)

    e <- buildNetwork(input.sequences = seqs, ids = ids,
                      threshold = thr, dist_type = "levenshtein")
    got_key <- if (nrow(e)) sort(paste(pmin(e$from, e$to), pmax(e$from, e$to))) else character(0)

    expect_identical(got_key, truth_key, info = paste("threshold", thr))
    # distances reported must equal the true edit distance
    if (nrow(e)) {
      td <- D[cbind(match(e$from, ids), match(e$to, ids))]
      expect_identical(as.integer(e$dist), as.integer(td), info = paste("dist values, threshold", thr))
    }
  }
})

test_that("identical sequences always edge, even at maxd 0", {
  d <- data.frame(sequence_id = c("a", "b", "c"),
                  junction_aa = c("ACDEFGHIK", "ACDEFGHIK", "ACDEFGHIK"),
                  v_call = "V1", stringsAsFactors = FALSE)
  # normalize maxlen with a tiny threshold -> maxd = floor(0.05 * 9) = 0
  e <- buildNetwork(d, seq_col = "junction_aa", ids = d$sequence_id,
                    threshold = 0.05, normalize = "maxlen", filter.v = TRUE)
  expect_equal(nrow(e), 3L)          # all three identical pairs
  expect_true(all(e$dist == 0))
})
