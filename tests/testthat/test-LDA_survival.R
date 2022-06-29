test_that("form", {
  data(LDAdata)
  cell.line <- unique(LDAdata$name)[1]
  x <- subset.data.frame(LDAdata, subset = (name==cell.line))
  sf <- LDA_survival(x[,1:4])
  expect_equal(class(sf),"list")
  expect_error(LDA_survival(list(x[,1:4])))
  expect_error(LDA_survival(x))
})
