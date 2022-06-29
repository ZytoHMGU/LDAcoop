test_that("plot", {
  data(LDAdata)
  cell.line <- unique(LDAdata$name)[1]
  x <- subset.data.frame(LDAdata, subset = (name==cell.line))
  act <- LDA_activity(x[,1:4])
  expect_error(LDA_plot(x[,1:4]))
  LDA_plot(act)
  x <- subset.data.frame(LDAdata, subset = (name==cell.line) & (Group==0))
  act <- LDA_activity_single(x[,1:3])
  LDA_plot(act)
})
