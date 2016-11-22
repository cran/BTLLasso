

#' BTLLasso
#' 
#' Performs BTLLasso, a method to model heterogeneity in paired comparison
#' data. Different types of covariates are allowd to have an influence on the
#' attractivity/strength of the objects. Covariates can be subject-specific, 
#' object-specific or subject-object-specific. L1 penalties are used to reduce the 
#' complexiy of the model by enforcing clusters of equal effects or by elimination of irrelevant
#' covariates.  Several additional functions are provided, such as
#' cross-validation, bootstraped confidence intervals, and plot
#' functions.
#' 
#' 
#' @name BTLLasso-package
#' @docType package
#' @author Gunther Schauberger\cr \email{gunther@@stat.uni-muenchen.de}\cr
#' \url{http://www.statistik.lmu.de/~schauberger/}
#' @seealso \code{\link{BTLLasso}}, \code{\link{cv.BTLLasso}}
#' @references Schauberger, Gunther and Tutz, Gerhard (2015): Modelling
#' Heterogeneity in Paired Comparison Data - an L1 Penalty Approach with an
#' Application to Party Preference Data, \emph{Department of Statistics, LMU
#' Munich}, Technical Report 183
#' @keywords package BTL Bradley-Terry BTLLasso
#' @examples
#' 
#' \dontrun{
#' ##### Example with simulated data set containing X, Z1 and Z2
#' data(SimData)
#' 
#' ## Specify tuning parameters
#' lambda <- exp(seq(log(151),log(1.05),length=30))-1
#' 
#' ## Specify control argument, allow for object-specific order effects and penalize intercepts
#' ctrl <- ctrl.BTLLasso(penalize.intercepts = TRUE, object.order.effect = TRUE,
#'                       penalize.order.effect.diffs = TRUE)
#' 
#' ## Simple BTLLasso model for tuning parameters lambda
#' m.sim <- BTLLasso(Y = SimData$Y, X = SimData$X, Z1 = SimData$Z1, 
#'                   Z2 = SimData$Z2, lambda = lambda, control = ctrl)
#' 
#' singlepaths(m.sim, x.axis = "loglambda")
#' 
#' ## Cross-validate BTLLasso model for tuning parameters lambda
#' set.seed(5)
#' m.sim.cv <- cv.BTLLasso(Y = SimData$Y, X = SimData$X, Z1 = SimData$Z1, 
#'                         Z2 = SimData$Z2, lambda = lambda, control = ctrl)
#' 
#' 
#' singlepaths(m.sim.cv, x.axis = "loglambda", plot.order.effect = FALSE, plot.intercepts = FALSE, 
#'             plot.Z2 = FALSE)
#' paths(m.sim.cv, y.axis="L2")
#' 
#' ## Example for bootstrap confidence intervals for illustration only
#' ## Don't calculate bootstrap confidence intervals with B = 10
#' set.seed(5)
#' m.sim.boot <- boot.BTLLasso(m.sim.cv, B = 10, cores = 10)
#' ci.BTLLasso(m.sim.boot)
#' 
#' ##### Example with small version from GLES data set
#' data(GLESsmall)
#' 
#' ## vector of subtitles, containing the coding of the X covariates
#' subs.X <- c("","female (1); male (0)")
#' 
#' ## vector of tuning parameters
#' lambda <- exp(seq(log(61),log(1),length=30))-1
#' 
#' 
#' ## compute BTLLasso model 
#' m.gles <- BTLLasso(Y = GLESsmall$Y, X = GLESsmall$X, Z1 = GLESsmall$Z1, lambda = lambda)
#' 
#' singlepaths(m.gles, x.axis = "loglambda", subs.X = subs.X)
#' paths(m.gles, y.axis="L2")
#' 
#' ## Cross-validate BTLLasso model 
#' m.gles.cv <- cv.BTLLasso(Y = GLESsmall$Y, X = GLESsmall$X, Z1 = GLESsmall$Z1, lambda = lambda)
#' 
#' singlepaths(m.gles.cv, x.axis = "loglambda", subs.X = subs.X)
#' }
#' 
NULL


#' Bundesliga Data 2015/16 (Buli1516)
#' 
#' Data from the German Bundesliga from the season 2015/16. 
#' The data contain all 306 matches of the season treated as paired comparisons with 5 different 
#' response categories. Additionally, different match-specific covariates are given as, for example, 
#' the percentage of ball possession or the total distance ran per team and per match.
#' 
#' @name Buli1516
#' @docType data
#' @format A list containing data from the German Bundesliga with 306 observations. 
#' The list contains both information on the response (paired comparisons) and different covariates.
#' \describe{ 
#' \item{Y}{A response.BTLLasso object for the Buli1516 data including
#' \itemize{
#' \item{response: Ordinal paired comparison response vector} 
#' \item{first.object: Vector containing the first-named team per paired comparison (home team)}
#' \item{second.object: Vector containing the second-named team per paired comparison (away team)}
#' \item{subject: Vector containing a match-day identifier per paired comparison}
#' }}
#' \item{Z1}{Matrix containing all team-match-specific covariates
#' \itemize{
#' \item{Distance: Total amount of km run} 
#' \item{BallPossession: Percentage of ball possession}
#' \item{TacklingRate: Rate of won tacklings}
#' \item{ShotsonGoal: Total number of shots on goal} 
#' \item{CompletionRate: Percentage of passes reaching teammates} 
#' \item{FoulsSuffered: Number of fouls suffered} 
#' \item{Offside: Number of offsides (in attack)}
#' }
#' }
#' \item{Z2}{Matrix containing all the average market values of the teams as a team-specific covariate} 
#' }
#' @references 
#' Schauberger, Gunther and Tutz, Gerhard (2015): Modelling Heterogeneity in
#' Paired Comparison Data - an L1 Penalty Approach with an Application to Party
#' Preference Data, \emph{Department of Statistics, LMU Munich}, Technical
#' Report 183
#' @source
#' \url{http://www.kicker.de/}
#' @keywords datasets
#' @examples
#' 
#' data(Buli1516)
#' 
NULL




#' German Longitudinal Election Study (GLES)
#' 
#' Data from the German Longitudinal Election Study (GLES), see Rattinger et
#' al. (2014). The GLES is a long-term study of the German electoral process.
#' It collects pre- and post-election data for several federal elections, the
#' data used here originate from the pre-election study for 2013.
#' 
#' @name GLES
#' @docType data
#' @format A list containing data from the German Longitudinal Election Study with 2003 
#' (partly incomplete) observations. 
#' The list contains both information on the response (paired comparisons) and different covariates.
#' \describe{ 
#' \item{Y}{A response.BTLLasso object for the GLES data including
#' \itemize{
#' \item{response: Ordinal paired comparison response vector} 
#' \item{first.object: Vector containing the first-named party per paired comparison}
#' \item{second.object: Vector containing the second-named party per paired comparison}
#' \item{subject: Vector containing a person identifier per paired comparison}
#' }}
#' \item{X}{Matrix containing all eight person-specific covariates
#' \itemize{
#' \item{Age: Age in years} 
#' \item{Gender (0: male, 1: female)}
#' \item{EastWest (0: West Germany, 1: East Germany)}
#' \item{PersEcon: Personal economic situation, 1: good or very good,
#' 0: else} 
#' \item{Abitur: School leaving certificate, 1: Abitur/A
#' levels, 0: else} 
#' \item{Unemployment: 1: currently unemployed, 0:
#' else} 
#' \item{Church: Frequency of attendence in a
#' church/synagogue/mosque/..., 1: at least once a month, 0: else}
#' \item{Migration: Have you been a German citizen since birth? 1: yes,
#' 0: no} 
#' }
#' }
#' \item{Z1}{Matrix containing all four person-party-specific covariates
#' \itemize{
#' \item{Climate: Self-perceived distance of each person to all five parties with respect to 
#' ones attitude towards climate change.}
#' \item{SocEc: Self-perceived distance of each person to all five parties with respect to 
#' ones attitude towards socio-economic issues.}
#' \item{Immigration: Self-perceived distance of each person to all five parties with respect to 
#' ones attitude towards immigration.}
#' }
#' }
#' }
#' @references Rattinger, H., S. Rossteutscher, R. Schmitt-Beck, B. Wessels,
#' and C. Wolf (2014): Pre-election cross section (GLES 2013). \emph{GESIS Data
#' Archive, Cologne ZA5700 Data file Version 2.0.0.}
#' 
#' Schauberger, Gunther and Tutz, Gerhard (2015): Modelling Heterogeneity in
#' Paired Comparison Data - an L1 Penalty Approach with an Application to Party
#' Preference Data, \emph{Department of Statistics, LMU Munich}, Technical
#' Report 183
#' @source
#' \url{http://gles.eu/wordpress/english/}
#' @keywords datasets
#' @examples
#' 
#' data(GLES)
#' 
NULL





#' Subset of the GLES data set with 200 observations and 4 covariates.
#' 
#' This is a subset of the \code{\link{GLES}} data set from the German
#' Longitudinal Election Study (GLES), see Rattinger et al. (2014). The subset contains 
#' only 200 of the 2003 observations and only  a small part of the covariates. The GLES is
#' a long-term study of the German electoral process. It collects pre- and
#' post-election data for several federal elections, the data used here
#' originate from the pre-election study for 2013.
#' 
#' @name GLESsmall
#' @docType data
#' @format A list containing data from the German Longitudinal Election Study with 200 observations. 
#' The list contains both information on the response (paired comparisons) and different covariates.
#' \describe{ 
#' \item{Y}{A response.BTLLasso object for the GLES data including
#' \itemize{
#' \item{response: Ordinal paired comparison response vector} 
#' \item{first.object: Vector containing the first-named party per paired comparison}
#' \item{second.object: Vector containing the second-named party per paired comparison}
#' \item{subject: Vector containing a person identifier per paired comparison}
#' }}
#' \item{X}{Matrix containing all eight person-specific covariates
#' \itemize{
#' \item{Age: Age in years} 
#' \item{Gender (0: male, 1: female)}
#' }
#' }
#' \item{Z1}{Matrix containing all four person-party-specific covariates
#' \itemize{
#' \item{Climate: Self-perceived distance of each person to all five parties with respect to 
#' ones attitude towards climate change.}
#' \item{Immigration: Self-perceived distance of each person to all five parties with respect to 
#' ones attitude towards immigration.}
#' }
#' }
#' }
#' @references Rattinger, H., S. Rossteutscher, R. Schmitt-Beck, B. Wessels,
#' and C. Wolf (2014): Pre-election cross section (GLES 2013). \emph{GESIS Data
#' Archive, Cologne ZA5700 Data file Version 2.0.0.}
#' 
#' Schauberger, Gunther and Tutz, Gerhard (2015): Modelling Heterogeneity in
#' Paired Comparison Data - an L1 Penalty Approach with an Application to Party
#' Preference Data, \emph{Department of Statistics, LMU Munich}, Technical
#' Report 183
#' @source
#' \url{http://gles.eu/wordpress/english/}
#' @keywords datasets
#' @examples
#' 
#' \dontrun{
#' data(GLESsmall)
#' 
#' ## vector of subtitles, containing the coding of the X covariates
#' subs.X <- c("","female (1); male (0)")
#' 
#' ## vector of tuning parameters
#' lambda <- exp(seq(log(61),log(1),length=30))-1
#' 
#' 
#' ## compute BTLLasso model 
#' m.gles <- BTLLasso(Y = GLESsmall$Y, X = GLESsmall$X, Z1 = GLESsmall$Z1, lambda = lambda)
#' 
#' singlepaths(m.gles, x.axis = "loglambda", subs.X = subs.X)
#' paths(m.gles, y.axis="L2")
#' 
#' ## Cross-validate BTLLasso model 
#' m.gles.cv <- cv.BTLLasso(Y = GLESsmall$Y, X = GLESsmall$X, Z1 = GLESsmall$Z1, lambda = lambda)
#' 
#' singlepaths(m.gles.cv, x.axis = "loglambda", subs.X = subs.X)
#' }
#' 
NULL

#' Simulated Data Set for illustration
#' 
#' This data set is a simulated data set including all possible types of covariates (X, Z1 and Z2)
#' and is intended to serve for illustration purpose. The data set contains paired comparisons between
#' four objects with five different response categories from 200 subjects.
#' 
#' @name SimData
#' @docType data
#' @format A list containing  simulated data for 200 observations. 
#' The list contains both information on the response (paired comparisons) and different covariates.
#' \describe{ 
#' \item{Y}{A response.BTLLasso object with simulated responses including
#' \itemize{
#' \item{response: Ordinal paired comparison response vector} 
#' \item{first.object: Vector containing the first-named object per paired comparison}
#' \item{second.object: Vector containing the second-named object per paired comparison}
#' \item{subject: Vector containing a subject identifier per paired comparison}
#' }}
#' \item{X}{Matrix containing both subject-specific covariates
#' \itemize{
#' \item{X_var1} 
#' \item{X_var2}
#' }
#' }
#' \item{Z1}{Matrix containing both subject-object-specific covariates
#' \itemize{
#' \item{Z1_var1}
#' \item{Z1_var2}
#' }
#' }
#' \item{Z2}{Matrix containing both object-specific covariates
#' \itemize{
#' \item{Z2_var1}
#' \item{Z2_var2}
#' }
#' }
#' }
#' @keywords datasets
#' @examples
#' 
#' \dontrun{
#' ##### Example with simulated data set containing X, Z1 and Z2
#' data(SimData)
#' 
#' ## Specify tuning parameters
#' lambda <- exp(seq(log(151),log(1.05),length=30))-1
#' 
#' ## Specify control argument, allow for object-specific order effects and penalize intercepts
#' ctrl <- ctrl.BTLLasso(penalize.intercepts = TRUE, object.order.effect = TRUE,
#'                       penalize.order.effect.diffs = TRUE)
#' 
#' ## Simple BTLLasso model for tuning parameters lambda
#' m.sim <- BTLLasso(Y = SimData$Y, X = SimData$X, Z1 = SimData$Z1, 
#'                   Z2 = SimData$Z2, lambda = lambda, control = ctrl)
#' 
#' singlepaths(m.sim, x.axis = "loglambda")
#' 
#' ## Cross-validate BTLLasso model for tuning parameters lambda
#' set.seed(5)
#' m.sim.cv <- cv.BTLLasso(Y = SimData$Y, X = SimData$X, Z1 = SimData$Z1, 
#'                         Z2 = SimData$Z2, lambda = lambda, control = ctrl)
#' 
#' 
#' singlepaths(m.sim.cv, x.axis = "loglambda", plot.order.effect = FALSE, plot.intercepts = FALSE, 
#'             plot.Z2 = FALSE)
#' paths(m.sim.cv, y.axis="L2")
#' }
#' 
NULL


