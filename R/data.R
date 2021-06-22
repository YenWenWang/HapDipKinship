#' @title Pedigrees
#'
#' @description A data set with parents of individuals.
#'
#' @format A data frame with 42 rows and 4 variables:
#' \describe{
#'   \item{id}{individuals' id.}
#'   \item{father}{individual's father.}
#'   \item{mother}{individual's mother}
#'   \item{sex}{individuals' sex.}
#' }
#' @source <https://www.github.com/YenWenWang/HapDipKinship>
"pedigree"

#' @title Relatedness of interest
#'
#' @description A dataset denoting which relatives are of interest (for simulation).
#'
#' @format A data frame with 164 rows and 4 variables:
#' \describe{
#'   \item{ind1}{individual 1.}
#'   \item{ind1}{individual 2.}
#'   \item{relationship}{manually coded relationship.}
#'   \item{pairtype}{ploidy of ind1 and ind2.}
#' }
#' @source <https://www.github.com/YenWenWang/HapDipKinship>
"RelOfInterest"
