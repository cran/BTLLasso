#' Create response object for BTLLasso
#' 
#' Create a response object for \code{BTLLasso} and \code{cv.BTLLasso}
#' 
#' 
#' @param response Vector containing results (binary or ordinal) of single paired
#' comparisons.
#' @param first.object Vector (character or factor, same length as response) indicating the first
#' object of the respective paired comparison from response.
#' @param second.object Vector (character or factor, same length as response) indicating the second
#' object of the respective paired comparison from response.
#' @param subject Vector (character, same length as response) indicating the subject that
#' generated the respective paired comparison from response.
#' @return Object of class \code{response.BTLLasso}
#' @author Gunther Schauberger\cr \email{gunther@@stat.uni-muenchen.de}\cr
#' \url{http://www.statistik.lmu.de/~schauberger/}
#' @seealso \code{\link{BTLLasso}}, \code{\link{cv.BTLLasso}}
#' @references Schauberger, Gunther and Tutz, Gerhard (2015): Modelling
#' Heterogeneity in Paired Comparison Data - an L1 Penalty Approach with an
#' Application to Party Preference Data, \emph{Department of Statistics, LMU
#' Munich}, Technical Report 183
#' 
#' Schauberger, Gunther, Groll Andreas and Tutz, Gerhard (2016): Modelling 
#' Football Results in the German Bundesliga Using Match-specific Covariates, 
#' \emph{Department of Statistics, LMU Munich}, Technical Report 197
#' @export response.BTLLasso
response.BTLLasso <- function(response, first.object, second.object, 
  subject = NULL) {
  
  withS <- FALSE
  if (!is.null(subject)) {
    withS <- TRUE
    if (!is.character(subject)) 
      stop("Argument subject has to be a character vector")
  }
  
  if (!withS) {
    subject <- 1:length(response)
  }
  
  ly <- length(response)
  lo1 <- length(first.object)
  lo2 <- length(second.object)
  ls <- length(subject)
  
  if (!all(sapply(list(lo1, lo2, ls), identical, ly))) 
    stop("The arguments response, first.object, second.object and (if specified) subject
     have to be of the same length")
  
  
  all.objects <- as.factor(as.character(unlist(list(first.object, 
    second.object))))
  object.names <- levels(all.objects)
  
  first.object <- as.numeric(all.objects[1:ly])
  second.object <- as.numeric(all.objects[(ly + 1):(2 * ly)])
  
  m <- length(object.names)
  
  ## make response ordered
  response <- as.ordered(response)
  
  # number of response categories
  q <- length(levels(response)) - 1
  k <- q + 1
  
  
  ## everything about the subjects
  subject.names <- levels(as.factor(subject))
  n <- length(subject.names)
  
  
  RET <- list(response = response, first.object = first.object, 
    second.object = second.object, subject = subject, withS = withS, 
    subject.names = subject.names, object.names = object.names, 
    n = n, m = m, k = k, q = q)
  
  class(RET) <- "responseBTLLasso"
  
  RET
}
