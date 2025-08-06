library(testthat)
library(devtools)
suppressWarnings(devtools::load_all())

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
  expect_error(cluster_GPU(condor,data_slot="bla"), "not exist")
})

test_that("Checking if reacts to missing PCA", {
  expect_error(cluster_GPU(condor_raw),"not computed")
})


test_that("Checking if reacts to missing rapids", {
  expect_error(cluster_GPU(condor,rapids_dir="bla"),"not found")
})


test_that("Checking if reacts to illegal values", {
  #expect_error(cluster_GPU(condor,GPU_device =1.5), "not integer") ##device allocation currently not functional
  expect_error(cluster_GPU(condor,n_neighbor=1.5), "a null pointer")
  expect_error(cluster_GPU(condor,n_neighbor=-5), "a positive integers")
  expect_error(cluster_GPU(condor,res="a"), "be real number")

})

test_that("Checking if reacts to illegal values in n_neighbor", {
  expect_snapshot_error(cluster_GPU(condor,n_neighbor="a"))

})

test_that("Checking if reacts to illegal values in seed", {
expect_snapshot_error(cluster_GPU(condor,seed="a"))
})



test_that("Checking if object has the correct metadata column", {
  expect_true("leiden_0.6"%in%names(cluster_GPU(condor)$clustering))
})
w
