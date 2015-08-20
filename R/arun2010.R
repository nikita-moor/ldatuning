# Implementation of the method, described in:
#
# Arun, R.; Suresh, V.; Veni Madhavan, C. E. & Narasimha Murthy, M. N. On
# Finding the Natural Number of Topics with Latent Dirichlet Allocation: Some
# Observations. Zaki, M. J.; Yu, J. X.; Ravindran, B. & Pudi, V. (Eds.) In
# Advances in Knowledge Discovery and Data Mining, Springer Berlin Heidelberg,
# 2010, 6118, 391-402
#
# ToDo:
# 1. make SLAM library optional (as optimization for sparse-matrix)
#    Is slam used by tm and topicmodels?
# 2. add "estimate.beta = TRUE" in control
# 3. dtm/tdm detection
# 4. It seems normalization of cm1 can be omitted

#' arun2010divergence
#'
#' symmetric Kullback-Leibler divergence
#' @keywords internal
arun2010divergence <- function(p,q) {
  entropy::KL.plugin(p,q) + entropy::KL.plugin(q,p)
}

#' arun2010
#'
#' Calculates arun2010 metric for a series of LDA models to estimate the most
#' preferable number of topics.
#'
#' @param dtm An object of class "\link[tm]{DocumentTermMatrix}" with
#'   term-frequency weighting or an object coercible to a
#'   "\link[slam]{simple_triplet_matrix}" with integer entries.
#' @param topics A list with number of topics to compare different models.
#' @param control A named list of the control parameters for estimation or an
#'   object of class "\linkS4class{LDAcontrol}".
#' @param mc.cores Integer; The number of CPU cores to processes models
#'   simultaneously (using \code{mclapply}).
#' @param verbose If false (default), supress all warnings and additional
#'   information.
#'
#' @return Data-frame with numbers of topics and corresponding values of metric
#'   (lower is better). Can be directly used by \code{\link{arun2010plot}} to
#'   draw a plot.
#'
#' @section References:
#'
#'   Arun, R., Suresh, V., Veni Madhavan, C. E. & Narasimha Murthy, M. N.
#'   (2010). On Finding the Natural Number of Topics with Latent Dirichlet
#'   Allocation: Some Observations. Zaki, M. J., Yu, J. X., Ravindran, B. &
#'   Pudi, V. (Eds.) In \emph{Advances in Knowledge Discovery and Data Mining},
#'   Springer Berlin Heidelberg, 6118, 391-402. DOI:
#'   \href{http://dx.doi.org/10.1073/10.1007/978-3-642-13657-3_43}{10.1007/978-3-642-13657-3_43}.
#'
#' @examples
#' library(topicmodels)
#' data("AssociatedPress", package="topicmodels")
#' dtm <- AssociatedPress[1:10, ]
#' arun2010(dtm, topics = 2:10, mc.cores = 2L)
#'
#' @export
arun2010 <- function(dtm, topics=seq(10, 40, by = 10), control=NULL,
                     mc.cores=2L, verbose=FALSE) {
  # check parameters
  if (length(topics[topics < 2]) != 0) {
    if (verbose) cat("warning: topics count can't to be less than 2, incorrect values was removed.\n")
    topics <- topics[topics >= 2]
  }
  # length of documents in words
  len <- slam::row_sums(dtm, dim=1)
  # fit models
  models <- parallel::mclapply(
    topics,
    FUN = function(x) topicmodels::LDA(dtm, k=x, method="Gibbs", control=control),
    mc.cores = mc.cores
  )
  # evaluate metrics
  metrics <- parallel::mclapply(
    models,
    FUN = function(model) {
      # matrix M1 topic-word
      m1 <- exp(model@beta) # rowSums(m1) == 1
      m1.svd <- svd(m1)
      cm1 <- as.matrix(m1.svd$d)
      # matrix M2 document-topic
      m2   <- model@gamma   # rowSums(m2) == 1
      cm2  <- len %*% m2    # crossprod(len, m2)
      norm <- norm(as.matrix(len), type="F")
      cm2  <- as.vector(cm2 / norm)
      return ( arun2010divergence(cm1, cm2) )
    },
    mc.cores = mc.cores
  )
  result <- data.frame(topics, as.numeric(metrics))
  names(result) <- c("topics", "arun2010")
  return(result)
}

#' arun2010plot
#'
#' Is a supplementary function to illustrate results of \code{\link{arun2010}}
#' function.
#'
#' @param results A data.frame with first column "topics" and second "arun2010".
#' @return Returns \code{\link[ggplot2]{ggplot}} object, which can be print or
#'   used for subsequent tuning by the user desire.
#'
#' @examples
#' library(topicmodels)
#' data("AssociatedPress", package="topicmodels")
#' dtm <- AssociatedPress[1:10, ]
#' metrics <- arun2010(dtm, topics = 2:10, mc.cores = 2L)
#' arun2010plot(metrics)
#'
#' @export
arun2010plot <- function(results) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("function 'arun2010plot' requires 'ggplot2' library.")
  }
  p <- ggplot2::ggplot(data=results, ggplot2::aes_string(x="topics", y="arun2010"))
  p <- p + ggplot2::geom_line()
  p <- p + ggplot2::labs(
    title = "Number of topics (by Arun2010)",
    x = "number of topics",
    y = "symmetric Kullback-Leiber divergence"
  )
  p <- p + ggplot2::scale_x_continuous(breaks=results$topics)
  p + ggplot2::theme_bw()
}
