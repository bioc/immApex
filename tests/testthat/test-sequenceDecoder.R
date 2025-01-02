# test script for sequenceDecoder.R - testcases are NOT comprehensive!

test_that("sequenceDecoder works", {
  
  sequence.matrix <- getdata("ohe.encoder", "onehotEncoder_matrix")
  
  if(reticulate::py_module_available("numpy")) {
    decoded.ohe <- sequenceDecoder(sequence.matrix)
    
    expect_equal(
      decoded.ohe,
      getdata("sequenceDecoder", "sequenceDecoder_ohe")
    )
    
    sequence.matrix <- getdata("propertyEncoder", "propertyEncoder_AtchleyFactors_matrix")
   
    decoded.AF <- sequenceDecoder(sequence.matrix,
                                  encoder = "propertyEncoder",
                                  aa.method.to.use = "atchleyFactors", 
                                  call.threshold = 2)
    
    expect_equal(
      decoded.AF,
      getdata("sequenceDecoder", "sequenceDecoder_AF")
    )
    
    sequence.array <- getdata("propertyEncoder", "propertyEncoder_KideraFactors_array")
    
    decoded.KF <- sequenceDecoder(sequence.array,
                                  encoder = "propertyEncoder",
                                  aa.method.to.use = "kideraFactors",
                                  call.threshold =1)
    
    expect_equal(
      decoded.KF,
      getdata("sequenceDecoder", "sequenceDecoder_KF")
    )
  }
})
