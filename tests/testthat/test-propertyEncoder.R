# test script for propertyEncoder.R - testcases are NOT comprehensive!

# 1. Check module availability once.
keras_installed <- reticulate::py_module_available("keras")
numpy_installed <- reticulate::py_module_available("numpy")

# 2. If not installed, skip everything:
if (!keras_installed || !numpy_installed) {
  test_that("Skipping propertyEncoder tests", {
    skip("Required Python modules (Keras, NumPy) are not available.")
  })
} else {
    
  test_that("propertyEncoder works", {
    
    sequences <- getdata("generateSequences", "generateSequences_T1")
      if(reticulate::py_module_available("numpy")) {
      #Return Matrix
      af.matrix <- propertyEncoder(sequences,
                                    method.to.use = "atchleyFactors",
                                    convert.to.matrix = TRUE)
      
      expect_equal(
        af.matrix,
        getdata("propertyEncoder", "propertyEncoder_AtchleyFactors_matrix")
      )
      
      #Return Array
      kf.array <- propertyEncoder(sequences,
                                   method.to.use = "kideraFactors",
                                   convert.to.matrix = FALSE)
      
      expect_equal(
        kf.array,
        getdata("propertyEncoder", "propertyEncoder_KideraFactors_array")
      )
      
      #Padded Matrix
      fasgai.matrix <- propertyEncoder(sequences,
                                        max.length = 40,
                                        method.to.use = "FASGAI",
                                        convert.to.matrix = FALSE)
      
      expect_equal(
        fasgai.matrix,
        getdata("propertyEncoder", "propertyEncoder_FASGAI_matrix")
      )
      
      #Padded Array
      vhse.array <- propertyEncoder(sequences,
                                     max.length = 40,
                                     method.to.use = "VHSE",
                                     convert.to.matrix = FALSE)
      
      expect_equal(
        vhse.array,
        getdata("propertyEncoder", "propertyEncoder_VSHE_array")
      )
      
      #Return Summary Matrix
      median.matrix <- propertyEncoder(sequences,
                                        method.to.use = "atchleyFactors",
                                        summary.function = "median")
      
      expect_equal(
        median.matrix,
        getdata("propertyEncoder", "propertyEncoder_AtchleyFactors_median.matrix")
      )
      
      #Return multiple properties
      multi.matrix <- propertyEncoder(sequences,
                                       method.to.use = c("atchleyFactors", "kideraFactors"))
      
      expect_equal(
        multi.matrix,
        getdata("propertyEncoder", "propertyEncoder_multi_matrix")
      )
      
      multi.array <- propertyEncoder(sequences,
                                      method.to.use = c("VHSE", "tScales"),
                                      convert.to.matrix = FALSE)
      
      expect_equal(
        multi.array,
        getdata("propertyEncoder", "propertyEncoder_multi_array")
      )
    }
  })
}