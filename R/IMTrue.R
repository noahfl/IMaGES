#' IMaGES Test Data
#'
#' Test data to be used with IMaGES. Included is
#' a list called \code{IMTrue}, a list containing the original
#' graph structures for each dataset in \code{IMTrue}.
#'
#' @docType data
#'
#' @usage data(IMData)
#'
#' @format Objects of class \code{"list"}.
#'
#' @keywords datasets
#'
#'
#' @examples
#' data(IMData)
#' data(IMTrue)
#' 
#' result <- IMaGES(matrices=IMTrue)
#' 
#' plotAll(result)
#' 
#' for (i in 1:length(IMTrue)) {
#'   plot(IMTrue[[i]])
#' }
#' 
#' 
"IMTrue"