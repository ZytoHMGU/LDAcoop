test_that("form", {
  data(LDAdata)
  cell.line <- unique(LDAdata$name)[1]
  x <- subset.data.frame(LDAdata, subset = (name==cell.line))
  act <- LDA_activity(x)
  expect_equal(class(act), "LDA_activity_list")
  expect_error(LDA_activity(list(x)))
  expect_warning(LDA_activity(x[,1:3]))
})
