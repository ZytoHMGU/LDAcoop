#' @title LDA_activity
#'
#' @description calculation of activity in a table of LDA data
#'   (i.e. dose, number of wells, number of positive wells, group).
#'
#' @param x numeric data.frame or matrix with three columns (dose, number of
#'   wells, number of positive wells, group (optional))
#' @param name optional: experiment name (e.g. name of cell line)
#'
#' @return list object with LDA-activities as returned by LDA_activity_single
#'
#' @examples
#' x <- data.frame("dose" = c(10,50,100,250,10,50,100,250),
#'                 "wells" = rep(25,8),
#'                 "positive" = c(2,5,10,20,1,2,6,11),
#'                 "group" = c(rep("A",4),rep("B",4)))
#' act <- LDA_activity(x)
#' @export
#'
LDA_activity <- function(x,name = "LDA cells"){
  if (!(class(x)[1] %in% c("data.frame","matrix"))){
    stop("error: x must be of class data.frame or matrix")
  }
  if (ncol(x) == 3){
    act <- LDA_activity_single(x,name)
  }
  if (ncol(x) > 3){
    x <- x[,1:4]
    colnames(x) <- c("dose","wells","positive","group")
    groups <- unique(x$group)
    act <- vector(mode = "list",length = length(groups))
    for (i in 1:length(groups)){
      act[[i]] <- LDA_activity_single(
        x = subset.data.frame(x = x,
                              subset = x$group == groups[i],
                              select = c("dose","wells","positive"),
                              drop = T),name,treat = groups[i])

    }
  }
  return(act)
}
