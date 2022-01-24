SRC_= file.choose("Escoge el archivo a probar")
setwd(file.path(dirname(SRC_)))
source(SRC_)
library(testthat)


test_that("DSLR", {
  
  library("xlsx")
  data_  <- read.xlsx("./rain_test.xls", sheetIndex = 1)
  RT = as.double(data_$Rain)
  expect_equal(1,DSLR_(1,RT,197))
  expect_equal(8,DSLR_(8,RT,197))
  expect_equal(5,DSLR_(93,RT,197))
})

test_that("INWS", {
  expect_equal(0,INSW(-1,0,1))
  expect_equal(1,INSW(1,0,1))
})

test_that("LIMIT", {
  expect_equal(2,LIMIT(1,2,3))
  expect_equal(3,LIMIT(3,2,1))
  expect_equal(2,LIMIT(1,3,2))
})

test_that("INTGRL <- function(IC,X,T,C)",{
        expect_equal(7,INTGRL(1,2,3,4))
})
  
test_that("AMOD <-function(X,P)",{
  expect_equal(1,AMOD(3,2))
  expect_equal(197,AMOD(197,365))
  
})
  
test_that("SUERRM <- function(MNR,X,XMIN,XMAX,NUNIT)",{
  expect_silent(SUERRM(0,2,1,3,6.))
  expect_error(SUERRM(1,0,1,3.,6.))
})

        
test_that("FUWCHK <-function(CKWFL,CKWIN,TIME,IDATE)",{
  
  expect_silent(FUWCHK(1,1,1,1))
  expect_output(FUWCHK(2,1,1,1))
  
})

test_that("Porcentaje de error",{
  expect_lt(porc_ , 5)
})







                    