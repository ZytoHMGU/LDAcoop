test_that("table", {
  data(LDAdata)
  cell.line <- unique(LDAdata$name)[1]
  x <- subset.data.frame(LDAdata, subset = (name==cell.line))
  expect_equal(class(LDA_table(x[,1:3])),"LDA_activity_object")
  x$`S-value` <- as.character(x$`S-value`)
  expect_error(LDA_table(x[,1:4]))
  x <- subset.data.frame(LDAdata, subset = (name==cell.line))
  LDA_table(x[,1:4])
})
test_that("strange things", {
  data(LDAdata)
  cell.line <- unique(LDAdata$name)[1]
  x <- subset.data.frame(LDAdata, subset = (name==cell.line))
  x.shuff <- rbind(x[x$Group!=0,],x[x$Group==0,])
  expect_warning(LDA_table(x.shuff[,1:4]))
})
