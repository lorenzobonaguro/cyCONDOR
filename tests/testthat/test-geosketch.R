library(testthat)
library(devtools)
devtools::load_all()
condor<-readRDS(here::here("../data/bonn_covid_test_data/condor_bonn_monocytes_factor.rds"))



test_that("Checking if reacts to non existent slot", {
  expect_error(subsample_geosketch(condor,"bla",10), "not exist")
})

test_that("Checking if reacts to wrong n_sub", {
  expect_error(subsample_geosketch(condor,"orig",nrow(condor$expr$orig)+1), "is higher")
  expect_error(subsample_geosketch(condor,"orig",(-10)), "an integer")
  expect_error(subsample_geosketch(condor,"orig","test"), "be an integer")
  expect_error(subsample_geosketch(condor,"orig",1.5), "fraction between")
})



condor_filter <- filter_fcd(fcd = condor,
                            cell_ids = rownames(condor$expr$orig)[1:20000])

test_that("Checking if object has expected dimensions after subset", {
  expect_equal(dim(subsample_geosketch(condor,"orig",20000)$expr$orig), dim(condor_filter$expr$orig))
})
