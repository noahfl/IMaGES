#' IMaGES Test Data
#'
#' Test data to be used with IMaGES. Included is
#' a list called \code{IMData}, a list containing three
#' data files that can be passed into IMaGES.
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
#' result <- IMaGES(matrices=IMData)
#' 
#' plotAll(result)
#' 
#' for (i in 1:length(IMTrue)) {
#'   plot(IMTrue[[i]])
#' }
#' 
#' 
"IMData"