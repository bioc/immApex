# test script for variationalSequences.R - testcases are NOT comprehensive!
# 1. Check module availability once.
keras_installed <- reticulate::py_module_available("keras")
numpy_installed <- reticulate::py_module_available("numpy")

# 2. If not installed, skip everything:
if (!keras_installed || !numpy_installed) {
  test_that("Skipping generateSequences tests", {
    skip("Required Python modules (Keras, NumPy) are not available.")
  })
} else {
  test_that("variationalSequences works", {
      set.seed(42)
      test_that("Test for valid input sequence length", {
        expect_error(variationalSequences(input.sequences = character(0)),
                     "input.sequences must have at least one sequence.")
      })
      
      test_that("Test for valid encoder function", {
        sequences <- c("CASGY", "CASDY", "CASTY")
        expect_error(variationalSequences(input.sequences = sequences,
                                          encoder.function = "invalidEncoder"),
                     "Invalid encoder provided.")
      })
      
      test_that("Test for valid optimizer", {
        sequences <- c("CASGY", "CASDY", "CASTY")
        expect_error(variationalSequences(input.sequences = sequences,
                                          optimizer = "invalidOptimizer"),
                     "Please select a compatible optimizer function in the Keras R implementation.")
      })
      
      test_that("Test for correct output type", {
        sequences <- getdata("generateSequences", "generateSequences_T1")[1:100]
        result <- variationalSequences(input.sequences = sequences, 
                                       sequence.dictionary = amino.acids[1:20])
        expect_true(is.vector(result))
      })
      
      test_that("Test for correct sequence generation", {
        sequences <- getdata("generateSequences", "generateSequences_T1")[1:100]
        number_of_sequences <- 5
        result <- variationalSequences(input.sequences = sequences, 
                                       number.of.sequences = number_of_sequences)
        expect_equal(length(result), number_of_sequences)
      })
  })
}