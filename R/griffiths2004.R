# Implementation of the method, described in:
#
# Griffiths, T. L. & Steyvers, M. Finding scientific topics Proceedings of the
# National Academy of Sciences, National Academy of Sciences, 2004, 101,
# 5228-5235
#
# Ponweiser, M. Latent Dirichlet Allocation in R. PhdThesis. Institute for
# Statistics and Mathematics, University of Economics and Business, 2012
#
#
# ToDo:
# 1. make Rmpfr library optional (looks like precision is enough)
# 2. dtm/tdm detection
# 3. It seems normalization of cm1 can be omitted
# 4. Hmisc

#' @keywords internal
griffiths2004harmonic <- function(logLikelihoods, precision=2000L) {
  # code is a little tricky, see explanation in Ponweiser2012
  # ToDo: add variant without "Rmpfr"
  llMed <- median(logLikelihoods)
  as.double( llMed - log( Rmpfr::mean( exp(
    -Rmpfr::mpfr(logLikelihoods, prec=precision) + llMed
  ))))
}

#' griffiths2004
#'
#' Calculates griffiths2004 metric for a series of LDA models to estimate the
#' most preferable number of topics.
#'
#' @inheritParams arun2010
#'
#' @return Data-frame with numbers of topics and corresponding values of metric
#'   (higher is better). Can be directly used by \code{\link{griffiths2004plot}}
#'   to draw a plot.
#'
#' @section References:
#'
#'   Griffiths, T. L. & Steyvers, M. Finding scientific topics (2004). In
#'   \emph{Proceedings of the National Academy of Sciences}, National Academy of
#'   Sciences, 101, 5228-5235. DOI:
#'   \href{http://dx.doi.org/10.1073/pnas.0307752101}{10.1073/pnas.0307752101}.
#'
#'   Ponweiser, M. Latent Dirichlet Allocation in R (2012). Diploma Thesis.
#'   Institute for Statistics and Mathematics, University of Economics and
#'   Business, Vienna. URL: \url{http://epub.wu.ac.at/id/eprint/3558}.
#'
#' @examples
#' library(topicmodels)
#' data("AssociatedPress", package="topicmodels")
#' dtm <- AssociatedPress[1:10, ]
#' griffiths2004(dtm, topics = 2:10, mc.cores = 2L)
#'
#' @export
griffiths2004 <- function(dtm, topics=seq(10, 40, by = 10), control=NULL,
                          mc.cores=2L, verbose=FALSE) {
  # check parameters
  if (length(topics[topics < 2]) != 0) {
    if (verbose) cat("warning: topics count can't to be less than 2, incorrect values was removed.\n")
    topics <- topics[topics >= 2]
  }
  if (is.null(control)) control = list()
  if (!"burnin" %in% names(control)) control$burnin <- 1000
  if (!"iter"   %in% names(control)) control$iter   <- 1000
  if (!"keep"   %in% names(control)) control$keep   <- 50
  # fit models
  models <- parallel::mclapply(
    topics,
    FUN = function(x) topicmodels::LDA(dtm, k=x, method="Gibbs", control=control),
    mc.cores = mc.cores
  )
  # log-likelihoods
  logLiks <- lapply(models, function(model) {
    model@logLiks[-c(1:(control$burnin / control$keep))]
  })
  # harmonic means
  harmony <- sapply(logLiks, function(x) griffiths2004harmonic(x))
  return(data.frame(topics, harmony))
}

#' griffiths2004plot
#'
#' Is a supplementary function to illustrate results of
#' \code{\link{griffiths2004}} function.
#'
#' @param results A data.frame with first column "topics" and second (or many)
#'   contains correspond harmonic mean values (griffiths2004 metrics).
#' @param show.max If true (default), find maximum value of metric and point it
#'   on the plot.
#' @param verbose If true, produce an additional information.
#' @return Returns \code{\link[ggplot2]{ggplot}} object, which can be print or
#'   used for subsequent tuning by the user desire.
#'
#' @examples
#' library(topicmodels)
#' data("AssociatedPress", package="topicmodels")
#' dtm <- AssociatedPress[1:10, ]
#' metrics <- griffiths2004(dtm, topics = 2:10, mc.cores = 2L)
#' griffiths2004plot(metrics)
#'
#' @export
#' @importFrom Hmisc smedian.hilow
#' @importFrom ggplot2 median_hilow
griffiths2004plot <- function(results, show.max = TRUE, verbose = FALSE) {
  # Function smedian.hilow (Hmisc) is wrapped by ggplot2, the above used
  # @importFrom works well, except for knitr (vignettes).
  if (!requireNamespace("reshape2", quietly = TRUE)) {
    stop("function 'griffiths2004plot' requires 'reshape2' package.")
  }
  p <- ggplot2::ggplot(
    reshape2::melt(results, id.vars = "topics"),
    ggplot2::aes_string(x = "topics", y = "value")
  )
  p <- p + ggplot2::geom_point()
  p <- p + ggplot2::labs(
    title = "Number of topics (by Griffiths2004)",
    x = "number of topics",
    y = "harmonic means"
  )
  p <- p + ggplot2::stat_summary(fun.data="median_hilow", geom="line", colour="black")
  p <- p + ggplot2::stat_summary(fun.data="median_hilow", geom="ribbon", alpha = 0.1)
  p <- p + ggplot2::scale_x_continuous(breaks=results$topics)

  # select the best number of topics
  if (show.max) {
    if (dim(results)[2] == 2) {
      # normal output of griffiths2004
      bestTopic <- results[ which.max(t(results[2])),"topics"]
    } else {
      # data is a result of multi-fitting
      median <- apply(results[,-1], 1, function(x) Hmisc::smedian.hilow(x))
      bestTopic <- which.max(median["Median",])
      bestTopic <- results[bestTopic,"topics"]
    }
    if (verbose) cat(sprintf("info: the recommend number of topics is %d.\n", bestTopic))
    p <- p + ggplot2::geom_vline(xintercept = bestTopic, linetype = "longdash")
  }
  p + ggplot2::theme_bw()
}
