# test script for onehotEncoder.R - testcases are NOT comprehensive!

# 1. Check module availability once.
keras_installed <- reticulate::py_module_available("keras")
numpy_installed <- reticulate::py_module_available("numpy")

# 2. If not installed, skip everything:
if (!keras_installed || !numpy_installed) {
  test_that("Skipping generateSequences tests", {
    skip("Required Python modules (Keras, NumPy) are not available.")
  })
} else {
  
  test_that("onehotEncoder works", {
   
    sequences <- getdata("generateSequences", "generateSequences_T1")
    if(reticulate::py_module_available("numpy")) {
      ohe.default <- onehotEncoder(sequences)
      
      expect_equal(
        ohe.default,
        getdata("ohe.encoder", "onehotEncoder_matrix")
      )
      
      ohe.2mer <- onehotEncoder(sequences, 
                                  motif.length = 2)
      expect_equal(
        ohe.2mer,
        getdata("ohe.encoder", "onehotEncoder_2mer.matrix")
      )
      
      ohe.padded <- onehotEncoder(sequences,
                                    max.length = 40)
      
      expect_equal(
        ohe.padded,
        getdata("ohe.encoder", "onehotEncoder_padded.matrix")
      )
      
      ohe.array <- onehotEncoder(sequences,
                                   convert.to.matrix = FALSE)
      
      expect_equal(
        ohe.array,
        getdata("ohe.encoder", "onehotEncoder_array")
      )
      
      ohe.padded.array <- onehotEncoder(sequences,
                                   max.length = 40,
                                   convert.to.matrix = FALSE)
      
      expect_equal(
        ohe.padded.array,
        getdata("ohe.encoder", "onehotEncoder_padded.array")
      )
      
      nt.sequences <- getdata("generateSequences", "generateSequences_T2")
      
      ohe.nt <- onehotEncoder(nt.sequences,
                                sequence.dictionary = c("A", "C", "T", "G"))
      
      expect_equal(
        ohe.nt,
        getdata("ohe.encoder", "onehotEncoder_nucleotide.matrix")
      )
    }
  })
}
