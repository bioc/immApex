# test script for sequenceDecoder.R - testcases are NOT comprehensive!

# 1. Check module availability once.
keras_installed <- reticulate::py_module_available("keras")
numpy_installed <- reticulate::py_module_available("numpy")

# 2. If not installed, skip everything:
if (!keras_installed || !numpy_installed) {
  test_that("Skipping generateSequences tests", {
    skip("Required Python modules (Keras, NumPy) are not available.")
  })
} else {
  
test_that("sequenceDecoder works", {
    
    sequence.matrix <- getdata("ohe.encoder", "onehotEncoder_matrix")
    
    decoded.ohe <- sequenceDecoder(sequence.matrix)
      
    expect_equal(
        decoded.ohe,
        getdata("sequenceDecoder", "sequenceDecoder_ohe"))
      
    sequence.matrix <- getdata("propertyEncoder", "propertyEncoder_AtchleyFactors_matrix")
     
    decoded.AF <- sequenceDecoder(sequence.matrix,
                                  encoder = "propertyEncoder",
                                  aa.method.to.use = "atchleyFactors", 
                                  call.threshold = 2)
      
    expect_equal(decoded.AF,
                  getdata("sequenceDecoder", "sequenceDecoder_AF"))
      
    sequence.array <- getdata("propertyEncoder", "propertyEncoder_KideraFactors_array")
      
    decoded.KF <- sequenceDecoder(sequence.array,
                                  encoder = "propertyEncoder",
                                  aa.method.to.use = "kideraFactors",
                                  call.threshold =1)
      
    expect_equal(decoded.KF,
                 getdata("sequenceDecoder", "sequenceDecoder_KF"))
  })
}
