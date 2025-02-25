# test script for tokenizeSequences.R - testcases are NOT comprehensive!

# 1. Check module availability once.
keras_installed <- reticulate::py_module_available("keras")
numpy_installed <- reticulate::py_module_available("numpy")

# 2. If not installed, skip everything:
if (!keras_installed || !numpy_installed) {
  test_that("Skipping tokenizeSequences tests", {
    skip("Required Python modules (Keras, NumPy) are not available.")
  })
} else {
  
  test_that("tokenizeSequences works", {
    
    sequences <- getdata("generateSequences", "generateSequences_T1")
  
    token.default <- tokenizeSequences(sequences)
    
    expect_equal(
      token.default,
      getdata("tokenizeSequences", "tokenizeSequences_default")
    )
    
    token.noAdded <- tokenizeSequences(sequences,
                                        add.startstop = FALSE)
    
    expect_equal(
      token.noAdded,
      getdata("tokenizeSequences", "tokenizeSequences_noStartStop")
    )
    
    token.matrix <- tokenizeSequences(sequences,
                                       convert.to.matrix = FALSE)
    
    expect_equal(
      token.matrix,
      getdata("tokenizeSequences", "tokenizeSequences_matrix")
    )
    
    
  })
}
