require(shiny)
library(debrowser)
library(DESeq2)
library(testthat)

load(paste0(system.file("extdata/demo/", package="debrowser"), "demodata.Rda"))
columns <- c("exper_rep1","exper_rep2","exper_rep3",
             "control_rep1","control_rep2","control_rep3")
conds <- factor( c("Control","Control", "Control",
                   "Treat", "Treat","Treat") )
data <- data.frame(demodata[, columns])

test_that("Able to run DESeq2", {
  deseqrun <- runDESeq(data, columns, conds)
  
  expect_true(exists("deseqrun"))
  expect_equal(deseqrun[[2]][[1]], 0.1641255385)
})
