# ---------------------------------------------------------------------------
# Regression harness for buildNetwork()
#
# Goal: freeze the *set* of edges (and their distances) that the current
# implementation produces, so that the candidate-indexing rewrite can be
# proven bit-identical in `expand = "clique"` mode across the full parameter
# matrix. Row order is intentionally ignored (OpenMP merge order is not
# deterministic); we compare canonical sets.
#
# Both the golden-snapshot generator (data-raw style script) and the
# regression test source these helpers, so the scenario grid lives in one
# place.
# ---------------------------------------------------------------------------

# Canonical, order-independent representation of a buildNetwork() result.
# Handles both the edge-list data.frame and the sparse dgCMatrix output.
.bn_canon <- function(res) {
  if (inherits(res, "data.frame")) {
    if (nrow(res) == 0L) return(character(0))
    a <- pmin(res$from, res$to)
    b <- pmax(res$from, res$to)
    # %.12g round-trips the integer and normalized doubles exactly while
    # staying robust to trivial formatting noise.
    key <- paste(a, b, sprintf("%.12g", res$dist), sep = "\031")
    sort(key)
  } else if (inherits(res, "dgCMatrix")) {
    res <- methods::as(res, "TsparseMatrix")
    nm <- rownames(res)
    fr <- nm[res@i + 1L]
    to <- nm[res@j + 1L]
    a <- pmin(fr, to)
    b <- pmax(fr, to)
    key <- paste(a, b, sprintf("%.12g", res@x), sep = "\031")
    sort(unique(key))
  } else {
    stop("Unexpected buildNetwork result class: ", class(res)[1])
  }
}

# Deterministic data generators (fixed seeds) covering the tricky cases.
.bn_data <- function() {
  aa <- c("A","C","D","E","F","G","H","I","K","L",
          "M","N","P","Q","R","S","T","V","W","Y")
  rseq <- function(len) paste(sample(aa, len, replace = TRUE), collapse = "")

  out <- list()

  # 1. Duplicate-heavy: many identical (seq,V,J) tuples across distinct ids,
  #    plus near-neighbors. This is the dedup/expansion stress case.
  set.seed(101)
  base <- replicate(15, rseq(sample(12:16, 1)))
  seqs <- sample(base, 120, replace = TRUE)          # heavy repetition
  # inject a handful of 1-2 edit neighbors of base[1]
  mut <- base[1]
  near <- vapply(1:6, function(k) {
    s <- strsplit(mut, "")[[1]]
    p <- sample(seq_along(s), 1); s[p] <- sample(aa, 1)
    paste(s, collapse = "")
  }, character(1))
  seqs <- c(seqs, near)
  n <- length(seqs)
  out$dup_heavy <- data.frame(
    sequence_id = paste0("dup_", seq_len(n)),
    junction_aa = seqs,
    v_call = paste0("IGHV", sample(1:3, n, replace = TRUE)),
    j_call = paste0("IGHJ", sample(1:2, n, replace = TRUE)),
    stringsAsFactors = FALSE)

  # 2. Mostly unique sequences, mixed lengths (length-blocking stress).
  set.seed(202)
  m <- 90
  out$unique_mixed <- data.frame(
    sequence_id = paste0("u_", seq_len(m)),
    junction_aa = vapply(seq_len(m), function(i) rseq(sample(10:20, 1)), character(1)),
    v_call = paste0("IGHV", sample(1:4, m, replace = TRUE)),
    j_call = paste0("IGHJ", sample(1:3, m, replace = TRUE)),
    stringsAsFactors = FALSE)

  # 3. Empty strings + single-member groups + short sequences.
  out$edge_cases <- data.frame(
    sequence_id = paste0("e_", 1:12),
    junction_aa = c("", "", "AC", "ACD", "ACDE", "ACDF",
                    "ACDEFGHIK", "ACDEFGHIL", "ACDEFGHIK",
                    "LONELYSEQ1", "WXYV", "WXYW"),
    v_call = c("V1","V1","V1","V1","V1","V1",
               "V2","V2","V2","V9","V3","V3"),
    j_call = c("J1","J1","J1","J1","J1","J1",
               "J1","J1","J1","J9","J1","J1"),
    stringsAsFactors = FALSE)

  # 4. Equal-length block for Hamming (with duplicates + 1/2 substitutions).
  set.seed(303)
  L <- 14L; q <- 60
  eq <- replicate(8, rseq(L))
  hs <- sample(eq, q, replace = TRUE)
  out$equal_len <- data.frame(
    sequence_id = paste0("h_", seq_len(q)),
    junction_aa = hs,
    v_call = paste0("IGHV", sample(1:2, q, replace = TRUE)),
    j_call = "IGHJ1",
    stringsAsFactors = FALSE)

  # 5. Pure transposition pairs (Damerau radius-2k trap). seq vs adjacent swap.
  out$transpose <- data.frame(
    sequence_id = paste0("t_", 1:8),
    junction_aa = c("ABCDEFGHIK", "BACDEFGHIK",   # swap pos 1-2  (Dam 1, Lev 2)
                    "ABCDEFGHIK", "ABDCEFGHIK",   # swap pos 3-4
                    "MNPQRSTVWY", "MNPQRSTVYW",    # swap pos 9-10
                    "MNPQRSTVWY", "NMPQRSTVWY"),   # swap pos 1-2
    v_call = "IGHV1",
    j_call = "IGHJ1",
    stringsAsFactors = FALSE)

  out
}

# The parameter grid. Each scenario is list(data=<name>, args=<list>).
# `args` are passed to buildNetwork() along with seq_col/ids.
.bn_scenarios <- function() {
  datasets <- .bn_data()
  grid <- list()
  add <- function(data_name, ...) {
    grid[[length(grid) + 1L]] <<- list(data = data_name, args = list(...))
  }

  filters <- list(c(FALSE, FALSE), c(TRUE, FALSE),
                  c(FALSE, TRUE), c(TRUE, TRUE))

  # Levenshtein / Damerau: absolute + normalized thresholds, all filters.
  for (dn in c("dup_heavy", "unique_mixed", "edge_cases")) {
    for (mt in c("levenshtein", "damerau")) {
      for (f in filters) {
        add(dn, dist_type = mt, threshold = 1, filter.v = f[1], filter.j = f[2])
        add(dn, dist_type = mt, threshold = 2, filter.v = f[1], filter.j = f[2])
        add(dn, dist_type = mt, threshold = 3, filter.v = f[1], filter.j = f[2])
        add(dn, dist_type = mt, threshold = 0.1, normalize = "maxlen",
            filter.v = f[1], filter.j = f[2])
        add(dn, dist_type = mt, threshold = 0.15, normalize = "length",
            filter.v = f[1], filter.j = f[2])
        # the loose-then-refilter pattern scRepertoire uses today
        add(dn, dist_type = mt, threshold = 0.85, normalize = "length",
            filter.v = f[1], filter.j = f[2])
      }
    }
  }

  # Hamming: equal-length and edge-case data, absolute + normalized.
  for (dn in c("equal_len", "edge_cases", "dup_heavy")) {
    for (f in filters) {
      add(dn, dist_type = "hamming", threshold = 1, filter.v = f[1], filter.j = f[2])
      add(dn, dist_type = "hamming", threshold = 2, filter.v = f[1], filter.j = f[2])
      add(dn, dist_type = "hamming", threshold = 3, filter.v = f[1], filter.j = f[2])
      add(dn, dist_type = "hamming", threshold = 0.15, normalize = "maxlen",
          filter.v = f[1], filter.j = f[2])
    }
  }

  # Damerau transposition stress.
  add("transpose", dist_type = "damerau", threshold = 1)
  add("transpose", dist_type = "damerau", threshold = 2)
  add("transpose", dist_type = "levenshtein", threshold = 1)
  add("transpose", dist_type = "levenshtein", threshold = 2)

  # NW / SW (no pruning path): a few thresholds.
  for (dn in c("unique_mixed", "equal_len")) {
    add(dn, dist_type = "nw", dist_mat = "BLOSUM62", threshold = 5, filter.v = TRUE)
    add(dn, dist_type = "sw", dist_mat = "BLOSUM62", threshold = 5, filter.v = TRUE)
    add(dn, dist_type = "nw", dist_mat = "PAM30", threshold = 8, filter.v = TRUE)
  }

  # A couple of sparse-output scenarios to lock the matrix path.
  add("dup_heavy", dist_type = "levenshtein", threshold = 2,
      filter.v = TRUE, output = "sparse", weight = "dist")
  add("dup_heavy", dist_type = "levenshtein", threshold = 2,
      filter.v = TRUE, output = "sparse", weight = "binary")

  list(datasets = datasets, grid = grid)
}

# Run one scenario through buildNetwork() and return its canonical set.
.bn_run <- function(datasets, scen, extra = list()) {
  d <- datasets[[scen$data]]
  args <- modifyList(scen$args, extra)
  args$input.data <- d
  args$seq_col <- "junction_aa"
  args$v_col <- "v_call"
  args$j_col <- "j_call"
  args$ids <- d$sequence_id
  res <- do.call(buildNetwork, args)
  .bn_canon(res)
}
