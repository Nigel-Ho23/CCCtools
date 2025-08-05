## code to prepare `DATASET` dataset goes here

set.seed(123)
testdata <- rnorm(1000)
usethis::use_data(testdata, overwrite = TRUE)

my_iris <- iris[1:50, ]
usethis::use_data(my_iris, overwrite = TRUE)
