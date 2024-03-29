# Copyright (c) 2017  Nikita Murzintcev

# ToDo:
# 1. add conversion: TermDocumentMatrix > DocumentTermMatrix
#    if (is(dtm, "TermDocumentMatrix")) dtm <- t(dtm)
# 2. CaoJuan2009: check with lsa::cosine - http://stackoverflow.com/questions/2535234/find-cosine-similarity-in-r


#' FindTopicsNumber
#'
#' Calculates different metrics to estimate the most preferable number of topics
#' for LDA model.
#'
#' @param dtm An object of class "\link[tm]{DocumentTermMatrix}" with
#'   term-frequency weighting or an object coercible to a
#'   "\link[slam]{simple_triplet_matrix}" with integer entries.
#' @param topics Vector with number of topics to compare different models.
#' @param metrics String or vector of possible metrics: "Griffiths2004",
#'   "CaoJuan2009", "Arun2010", "Deveaud2014".
#' @param method The method to be used for fitting; see \link[topicmodels]{LDA}.
#' @param control A named list of the control parameters for estimation or an
#'   object of class "\linkS4class{LDAcontrol}".
#' @param mc.cores NA, integer or, cluster; the number of CPU cores to process models
#'   simultaneously. If an integer, create a cluster on the local machine. If a
#'   cluster, use but don't destroy it (allows multiple-node clusters). Defaults to
#'   NA, which triggers auto-detection of number of cores on the local machine.
#' @param return_models Whether or not to return the model objects of class
#'   "\link[topicmodels]{LDA}. Defaults to false. Setting to true requires the tibble package.
#' @param verbose If false (default), suppress all warnings and additional
#'   information.
#' @param libpath Path to R packages (use only if your R installation can't find
#'   'topicmodels' package, [issue #3](https://github.com/nikita-moor/ldatuning/issues/3).
#'   For example: "C:/Program Files/R/R-2.15.2/library" (Windows),
#'                "/home/user/R/x86_64-pc-linux-gnu-library/3.2" (Linux)
#'
#' @return Data-frame with one or more metrics.  numbers of topics and
#'   corresponding values of metric. Can be directly used by
#'   \code{\link{FindTopicsNumber_plot}} to draw a plot.
#'
#' @examples
#' \dontrun{
#'
#' library(topicmodels)
#' data("AssociatedPress", package="topicmodels")
#' dtm <- AssociatedPress[1:10, ]
#' FindTopicsNumber(dtm, topics = 2:10, metrics = "Arun2010", mc.cores = 1L)
#' }
#'
#' @export
#' @import topicmodels
FindTopicsNumber <- function(dtm, topics = seq(10, 40, by = 10),
                             metrics = "Griffiths2004",
                             method = "Gibbs", control = list(),
                             mc.cores = NA, return_models = FALSE,
                             verbose = FALSE, libpath = NULL) {
  # check parameters
  if (length(topics[topics < 2]) != 0) {
    if (verbose) cat("warning: topics count can't to be less than 2, incorrect values was removed.\n")
    topics <- topics[topics >= 2]
  }
  topics <- sort(topics, decreasing = TRUE)

  if ("Griffiths2004" %in% metrics) {
    if (method == "VEM") {
      # memory allocation error
      if (verbose) cat("'Griffiths2004' is incompatible with 'VEM' method, excluded.\n")
      metrics <- setdiff(metrics, "Griffiths2004")
    } else {
      # save log-likelihood when generating model
      if (!"keep" %in% names(control)) control <- c(control, keep = 50)
    }
  }
  if ("Arun2010" %in% metrics) {
    if ( max(topics) > ncol(dtm) ) {
      if (verbose) cat("To use \'Arun2010\', the number of columns must be greater than or equal to the number of topics.\n")
      metrics <- setdiff(metrics, "Arun2010")
    }
  }

  # fit models
  if (verbose) cat("fit models...")

  # Parallel setup
  if (any(class(mc.cores) == "cluster")) {
    cl <- mc.cores
  } else if (isTRUE(class(mc.cores) == "integer")) {
    cl <- parallel::makeCluster(mc.cores)
  } else {
    cl <- parallel::makeCluster(parallel::detectCores())
  }
  parallel::setDefaultCluster(cl)
  parallel::clusterExport(varlist = c("dtm", "method", "control"),
                          envir = environment())
  models <- parallel::parLapply(cl = cl, X = topics, fun = function(x) {
  if (is.null(libpath) == FALSE) { .libPaths(libpath) }
    topicmodels::LDA(dtm, k = x, method = method, control = control)
  })
  if (! any(class(mc.cores) == "cluster")) {
    parallel::stopCluster(cl)
  }
  if (verbose) cat(" done.\n")

  # calculate metrics
  if (verbose) cat("calculate metrics:\n")

  if (return_models &
      requireNamespace("tibble", quietly = TRUE)
  ) {
    result <- cbind(topics, tibble::enframe(models, value = "LDA_model"))
    result$name <- NULL
  } else {
    if (return_models) {
      message("The tibble package is required for returning models. Returning results only.")
    }
    result <- data.frame(topics)
  }
  for(m in metrics) {
    if (verbose) cat(sprintf("  %s...", m))
    if (! m %in% c("Griffiths2004", "CaoJuan2009", "Arun2010", "Deveaud2014")) {
      cat(" unknown!\n")
    } else {
      result[m] <- switch(m,
        "Griffiths2004" = Griffiths2004(models, control),
        "CaoJuan2009"   = CaoJuan2009(models),
        "Arun2010"      = Arun2010(models, dtm),
        "Deveaud2014"   = Deveaud2014(models),
        NaN
      )
      if (verbose) cat(" done.\n")
    }
  }

  return(result)
}

#' Griffiths2004
#'
#' Implement scoring algorithm. In order to use this algorithm, the LDA model MUST
#' be generated using the keep control parameter >0 (defaults to 50) so that the
#' logLiks vector is retained.
#' @param models An object of class "\link[topicmodels]{LDA}
#' @param control A named list of the control parameters for estimation or an
#'   object of class "\linkS4class{LDAcontrol}".
#' @return A scalar LDA model score
#'
#' @export
#'
Griffiths2004 <- function(models, control) {
  # log-likelihoods (remove first burnin stage)
  burnin  <- ifelse("burnin" %in% names(control), control$burnin, 0)

  logLiks <- lapply(models, function(model) {
    # Check to make sure logLiks were kept; if not, value is NaN
    if (length(model@logLiks) == 0) {
      message("No logLiks were kept, which is required to use this scoring algorithm. Please regenerate the model using the keep control parameter set to a reasonable value (default = 50).")
      NaN
    } else {
      utils::tail(model@logLiks, n = length(model@logLiks) - burnin/control$keep)
      # model@logLiks[-(1 : (control$burnin/control$keep))]
    }
  })

  # harmonic means for every model
  metrics <- sapply(logLiks, function(x) {
    # code is a little tricky, see explanation in [Ponweiser2012 p. 36]
    # ToDo: add variant without "Rmpfr"
    llMed <- stats::median(x)
    metric <- as.double(
      llMed - log( Rmpfr::mean( exp( -Rmpfr::mpfr(x, prec=2000L) + llMed )))
    )
    return(metric)
  })
  return(metrics)
}

#' CaoJuan2009
#'
#' Implement scoring algorithm
#' @param models An object of class "\link[topicmodels]{LDA}
#' @return A scalar LDA model score
#'
#' @export
#'
CaoJuan2009 <- function(models) {
  metrics <- sapply(models, function(model) {
    # topic-word matrix
    m1 <- exp(model@beta)
    # pair-wise cosine distance
    pairs <- utils::combn(nrow(m1), 2)
    cos.dist <- apply(pairs, 2, function(pair) {
      x <- m1[pair[1], ]
      y <- m1[pair[2], ]
      # dist <- lsa::cosine(x, y)
      dist <- crossprod(x, y) / sqrt(crossprod(x) * crossprod(y))
      return(dist)
    })
    # metric
    metric <- sum(cos.dist) / (model@k*(model@k-1)/2)
    return(metric)
  })
  return(metrics)
}

#' Arun2010
#'
#' Implement scoring algorithm
#' @param models An object of class "\link[topicmodels]{LDA}
#' @param dtm An object of class "\link[tm]{DocumentTermMatrix}" with
#'   term-frequency weighting or an object coercible to a
#'   "\link[slam]{simple_triplet_matrix}" with integer entries.
#' @return A scalar LDA model score
#'
#' @export
#'
Arun2010 <- function(models, dtm) {
  # length of documents (count of words)
  len <- slam::row_sums(dtm)
  # evaluate metrics
  metrics <- sapply(models, FUN = function(model) {
    # matrix M1 topic-word
    m1 <- exp(model@beta) # rowSums(m1) == 1
    m1.svd <- svd(m1)
    cm1 <- as.matrix(m1.svd$d)
    # matrix M2 document-topic
    m2   <- model@gamma   # rowSums(m2) == 1
    cm2  <- len %*% m2    # crossprod(len, m2)
    norm <- norm(as.matrix(len), type="m")
    cm2  <- as.vector(cm2 / norm)
    # symmetric Kullback-Leibler divergence
    divergence <- sum(cm1*log(cm1/cm2)) + sum(cm2*log(cm2/cm1))
    return ( divergence )
  })
  return(metrics)
}

#' Deveaud2014
#'
#' Implement scoring algorithm
#' @param models An object of class "\link[topicmodels]{LDA}
#' @return A scalar LDA model score
#'
#' @export
#'
Deveaud2014 <- function(models) {
  metrics <- sapply(models, function(model) {
    ### original version
    # topic-word matrix
    m1 <- exp(model@beta)
    # prevent NaN
    if (any(m1 == 0)) { m1 <- m1 + .Machine$double.xmin }
    # pair-wise Jensen-Shannon divergence
    pairs  <- utils::combn(nrow(m1), 2)
    jsd <- apply(pairs, 2, function(pair) {
      x <- m1[pair[1], ]
      y <- m1[pair[2], ]
      ### standard Jensen-Shannon divergence
      # m <- (x + y) / 2
      # jsd <- 0.5 * sum(x*log(x/m)) + 0.5 * sum(y*log(y/m))
      ### divergence by Deveaud2014
      jsd <- 0.5 * sum(x*log(x/y)) + 0.5 * sum(y*log(y/x))
      return(jsd)
    })

#     ### optimized version
#     m1   <- model@beta
#     m1.e <- exp(model@beta)
#     pairs  <- utils::combn(nrow(m1), 2)
#     jsd <- apply(pairs, 2, function(pair) {
#       x   <- m1[pair[1], ]
#       y   <- m1[pair[2], ]
#       x.e <- m1.e[pair[1], ]
#       y.e <- m1.e[pair[2], ]
#       jsd <- ( sum(x.e*(x-y)) + sum(y.e*(y-x)) ) / 2
#       return(jsd)
#     })

    # metric
    metric <- sum(jsd) / (model@k*(model@k-1))
    return(metric)
  })
  return(metrics)
}


#' FindTopicsNumber_plot
#'
#' Support function to analyze optimal topic number. Use output of the
#' \code{\link{FindTopicsNumber}} function.
#'
#' @param values Data-frame with first column named `topics` and other columns
#'   are values of metrics.
#'
#' @examples
#' \dontrun{
#'
#' library(topicmodels)
#' data("AssociatedPress", package="topicmodels")
#' dtm <- AssociatedPress[1:10, ]
#' optimal.topics <- FindTopicsNumber(dtm, topics = 2:10,
#'   metrics = c("Arun2010", "CaoJuan2009", "Griffiths2004")
#' )
#' FindTopicsNumber_plot(optimal.topics)
#' }
#'
#' @export
#' @import ggplot2
FindTopicsNumber_plot <- function(values) {
  # Drop models if present, as they won't rescale
  if ("LDA_model" %in% names(values)) {
    values <- values[!names(values) %in% c("LDA_model")]
  }
  # normalize to [0,1]
  columns <- base::subset(values, select = 2:ncol(values))
  values <- base::data.frame(
    values["topics"],
    base::apply(columns, 2, function(column) {
      scales::rescale(column, to = c(0, 1), from = range(column))
    })
  )

  # melt
  values <- reshape2::melt(values, id.vars = "topics", na.rm = TRUE)

  # separate max-arg & min-arg metrics
  values$group <- values$variable %in% c("Griffiths2004", "Deveaud2014")
  values$group <- base::factor(
    values$group,
    levels = c(FALSE, TRUE),
    labels = c("minimize", "maximize")
  )

  # standard plot
  p <- ggplot(values, aes_string(x = "topics", y = "value", group = "variable"))
  p <- p + geom_line()
  p <- p + geom_point(aes_string(shape = "variable"), size = 3)
  p <- p + guides(size = FALSE, shape = guide_legend(title = "metrics:"))
  p <- p + scale_x_continuous(breaks = values$topics)
  p <- p + labs(x = "number of topics", y = NULL)

  # separate in two parts
  p <- p + facet_grid(group ~ .)

  # style
  # p <- p + theme_bw(base_size = 14, base_family = "") %+replace% theme(
  p <- p + theme_bw() %+replace% theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_line(colour = "grey70"),
    panel.grid.minor.x = element_blank(),
    legend.key = element_blank(),
    strip.text.y = element_text(angle = 90)
  )

  # move strip block to left side
  g <- ggplotGrob(p)
  g$layout[g$layout$name == "strip-right", c("l", "r")] <- 3
  grid::grid.newpage()
  grid::grid.draw(g)

  # return(p)
}
