# test script for generateSequences.R - testcases are NOT comprehensive!

# 1. Check module availability once.
keras_installed <- reticulate::py_module_available("keras")
numpy_installed <- reticulate::py_module_available("numpy")

# 2. If not installed, skip everything:
if (!keras_installed || !numpy_installed) {
  test_that("Skipping generateSequences tests", {
    skip("Required Python modules (Keras, NumPy) are not available.")
  })
} else {
  
  # 3. If installed, run tests as normal:
  test_that("generateSequences works", {
    
    set.seed(42)
    t1.sequences <- generateSequences(prefix.motif = "CAS",
                                      suffix.motif = "YF",
                                      number.of.sequences = 1000,
                                      min.length = 8,
                                      max.length = 16)
    
    expect_equal(
      t1.sequences,
      getdata("generateSequences", "generateSequences_T1")
    )
    
    set.seed(42)
    t2.sequences <- generateSequences(number.of.sequences = 1000,
                                      min.length = 2,
                                      max.length = 16, 
                                      sequence.dictionary = c("A", "C", "T", "G"))
    
    expect_equal(
      t2.sequences,
      getdata("generateSequences", "generateSequences_T2")
    )
  })
}