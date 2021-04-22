library(testthat)
library(devtools)
library(flintyR)

# either of these should work just ensure that you are in the proper working
# directory

#devtools::test(pkg="flintyR")
testthat::test_package("flintyR")
