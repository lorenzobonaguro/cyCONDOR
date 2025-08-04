library(testthat)
library(devtools)
devtools::load_all()

condor_raw <- prep_fcd(data_path = "/home/data/cyCONDOR_demo/fcs/", # Select the folder where the compensated fcs files are located
                   max_cell = 10000,  # Select the maximum number of cells to analyse for each file
                   useCSV = FALSE, # If csv files are used instead of fcs
                   transformation = "auto_logi", # Type of trasformation "a" = autologicle
                   remove_param = c("FSC-H", "SSC-H", "FSC-W", "SSC-W", "Time", "InFile"), # parameter to exclude
                   anno_table = "/home/data/cyCONDOR_demo/metadata.csv", # Metadata file for the annotation
                   filename_col = "filename", # Column of the metadata file containing the file name
                   seed = 91)

condor<-readRDS("/home/data/bonn_covid_test_data/condor_bonn_monocytes_factor.rds")


test_that("Checking if reacts to non existent slot", {
  expect_error(harmonize_PCA(condor,data_slot="bla",batch_var="sex"), "not found")
})

test_that("Checking if reacts to missing PCA", {
  expect_error(harmonize_PCA(condor_raw,batch_var="sex"),"not computed")
})

test_that("Checking if reacts to non existent meta column", {
  expect_error(harmonize_PCA(condor,batch_var="bla"), "not found in meta table")
})

test_that("Checking if switches to GPU", {
  expect_message(harmonize_PCA(condor,batch_var="sex",GPU=T), "rapids")
})

test_that("Checking if reacts to missing rapids", {
  expect_error(harmonize_PCA(condor,batch_var="sex",GPU=T,rapids_dir="bla"),"not found")
})



test_that("Checking if object has harmonised PCA", {
  expect_equal(names(harmonize_PCA(condor,batch_var="sex",GPU=T)$pca), c("orig","norm"))
})
