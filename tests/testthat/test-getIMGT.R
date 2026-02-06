# test script for getIMGT.R - testcases are NOT comprehensive!

test_that("getIMGT returns expected structure for amino acid sequences", {
  TRBV_human_aa <- getIMGT(species = "human",
                            chain = "TRB",
                            region = "v",
                            sequence.type = "aa")

  expect_type(TRBV_human_aa, "list")
  expect_true("sequences" %in% names(TRBV_human_aa))
  expect_true("misc" %in% names(TRBV_human_aa))
  expect_true(length(TRBV_human_aa$sequences) > 0)
  expect_equal(TRBV_human_aa$misc$chain, "TRB")
  expect_equal(TRBV_human_aa$misc$region, "v")
  expect_equal(TRBV_human_aa$misc$sequence.type, "aa")
})

test_that("getIMGT works with different regions and species", {
  TRBJ_mouse_nt <- getIMGT(species = "mouse",
                            chain = "TRB",
                            region = "j",
                            sequence.type = "nt")

  expect_type(TRBJ_mouse_nt, "list")
  expect_true(length(TRBJ_mouse_nt$sequences) > 0)
  expect_equal(TRBJ_mouse_nt$misc$region, "j")
})

test_that("getIMGT works with nucleotide sequences", {
  IGHV_rat_nt <- getIMGT(species = "rat",
                          chain = "IGH",
                          region = "v",
                          sequence.type = "nt")

  expect_type(IGHV_rat_nt, "list")
  expect_true(length(IGHV_rat_nt$sequences) > 0)
  expect_equal(IGHV_rat_nt$misc$sequence.type, "nt")
})

test_that("getIMGT works with additional species", {
  TRBV_rabbit_aa <- getIMGT(species = "rabbit",
                             chain = "TRB",
                             region = "v",
                             sequence.type = "aa")

  expect_type(TRBV_rabbit_aa, "list")
  expect_true(length(TRBV_rabbit_aa$sequences) > 0)
})

test_that("getIMGT input validation works", {
  expect_error(getIMGT(region = "x"), "Invalid region")
  expect_error(getIMGT(sequence.type = "protein"), "Invalid sequence.type")
})

test_that(".mapSpecies maps correctly", {
  expect_equal(.mapSpecies("human"), "human")
  expect_equal(.mapSpecies("Human"), "human")
  expect_equal(.mapSpecies("rhesus monkey"), "rhesus_monkey")
  expect_equal(.mapSpecies("Rhesus monkey"), "rhesus_monkey")
  expect_equal(.mapSpecies("rhesus Monkey"), "rhesus_monkey")
  expect_error(.mapSpecies("platypus"), "Invalid species")
})
