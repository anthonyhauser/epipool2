#' Simulated pooled test results
#'
#' @description A list of 5 datasets containing a simulated prevalence over time, as well as simulated
#' pool test results in 5 different regions.
#' The datasets were simulated assuming a sensitivity of 83% and specificity of 99.5%.
#' Each dataset contains the following variables:
#'
#' \itemize{
#'   \item time. Time (1--30)
#'   \item region. Regions (1--5)
#'   \item pool.size. Pool size (number of samples by pool, fixed at 10)
#'   \item n.pools. Total number of pools at given time point and region (fixed at 10)
#'   \item n.pos.pools. Number of positive pools (should be between 0 and 10).
#'   \item prev. True (simulated) overall prevalence for a given time.
#'   \item prev_region. True (simulated) prevalence for a given region and a given time.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name epipool_data1
#' @usage data(epipool_data1)
#' @format A list of 5 data frame with 150 rows and 7 variables
#' @examples
#' #load data
#' data(epipool_data1)
NULL

#' Simulated pooled test results
#'
#' @description A list of 5 datasets containing a simulated prevalence over time, as well as simulated
#' pool test results in 5 different regions.
#' The datasets were simulated assuming a sensitivity of 83% and specificity of 99.5%.
#' Each dataset contains the following variables:
#'
#' \itemize{
#'   \item time. Time (1--30)
#'   \item region. Regions (1--5)
#'   \item pool.size. Pool size (number of samples by pool, fixed at 10)
#'   \item n.pools. Total number of pools at given time point and region (fixed at 100)
#'   \item n.pos.pools. Number of positive pools (should be between 0 and 10).
#'   \item prev. True (simulated) overall prevalence for a given time.
#'   \item prev_region. True (simulated) prevalence for a given region and a given time.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name epipool_data2
#' @usage data(epipool_data2)
#' @format A list of 5 data frame with 150 rows and 7 variables
#' @examples
#' #load data
#' data(epipool_data2)
NULL

#' Simulated pooled test results
#'
#' @description A list of 5 datasets containing a simulated prevalence over time, as well as simulated
#' pool test results in 5 different regions.
#' The datasets were simulated assuming a sensitivity of 83% and specificity of 99.5%.
#' Each dataset contains the following variables:
#'
#' \itemize{
#'   \item time. Time (1--30)
#'   \item region. Regions (1--5)
#'   \item pool.size. Pool size (number of samples by pool, fixed at 10)
#'   \item n.pools. Total number of pools at given time point and region (fixed at 20)
#'   \item n.pos.pools. Number of positive pools (should be between 0 and 10).
#'   \item prev. True (simulated) overall prevalence for a given time.
#'   \item prev_region. True (simulated) prevalence for a given region and a given time.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name epipool_data3
#' @usage data(epipool_data3)
#' @format A list of 5 data frame with 150 rows and 7 variables
#' @examples
#' #load data
#' data(epipool_data3)
NULL
