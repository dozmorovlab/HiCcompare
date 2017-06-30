# Functions that are called within main functions of HiCLoess

# loess with Automatic Smoothing Parameter Selection adjusted possible
# range of smoothing originally from fANCOVA package
.loess.as <- function(x, y, degree = 1, criterion = c("aicc", "gcv"),
                      family = c("gaussian",
                                 "symmetric"), user.span = NULL, plot = FALSE, ...) {
  criterion <- match.arg(criterion)
  family <- match.arg(family)
  x <- as.matrix(x)

  if ((ncol(x) != 1) & (ncol(x) != 2))
    stop("The predictor 'x' should be one or two dimensional!!")
  if (!is.numeric(x))
    stop("argument 'x' must be numeric!")
  if (!is.numeric(y))
    stop("argument 'y' must be numeric!")
  if (any(is.na(x)))
    stop("'x' contains missing values!")
  if (any(is.na(y)))
    stop("'y' contains missing values!")
  if (!is.null(user.span) && (length(user.span) != 1 || !is.numeric(user.span)))
    stop("argument 'user.span' must be a numerical number!")
  if (nrow(x) != length(y))
    stop("'x' and 'y' have different lengths!")
  if (length(y) < 3)
    stop("not enough observations!")

  data.bind <- data.frame(x = x, y = y)
  if (ncol(x) == 1) {
    names(data.bind) <- c("x", "y")
  } else {
    names(data.bind) <- c("x1", "x2", "y")
  }

  opt.span <- function(model, criterion = c("aicc", "gcv"), span.range = c(0.01,
                                                                           0.9)) {
    as.crit <- function(x) {
      span <- x$pars$span
      traceL <- x$trace.hat
      sigma2 <- sum(x$residuals^2)/(x$n - 1)
      aicc <- log(sigma2) + 1 + 2 * (2 * (traceL + 1))/(x$n - traceL -
                                                          2)
      gcv <- x$n * sigma2/(x$n - traceL)^2
      result <- list(span = span, aicc = aicc, gcv = gcv)
      return(result)
    }
    criterion <- match.arg(criterion)
    fn <- function(span) {
      mod <- update(model, span = span)
      as.crit(mod)[[criterion]]
    }
    result <- optimize(fn, span.range)
    return(list(span = result$minimum, criterion = result$objective))
  }

  if (ncol(x) == 1) {
    if (is.null(user.span)) {
      fit0 <- loess(y ~ x, degree = degree, family = family, data = data.bind,
                    ...)
      span1 <- opt.span(fit0, criterion = criterion)$span
    } else {
      span1 <- user.span
    }
    fit <- loess(y ~ x, degree = degree, span = span1, family = family,
                 data = data.bind, ...)
  } else {
    if (is.null(user.span)) {
      fit0 <- loess(y ~ x1 + x2, degree = degree, family = family,
                    data.bind, ...)
      span1 <- opt.span(fit0, criterion = criterion)$span
    } else {
      span1 <- user.span
    }
    fit <- loess(y ~ x1 + x2, degree = degree, span = span1, family = family,
                 data = data.bind, ...)
  }
  if (plot) {
    if (ncol(x) == 1) {
      m <- 100
      x.new <- seq(min(x), max(x), length.out = m)
      fit.new <- predict(fit, data.frame(x = x.new))
      plot(x, y, col = "lightgrey", xlab = "x", ylab = "m(x)", ...)
      lines(x.new, fit.new, lwd = 1.5, ...)
    } else {
      m <- 50
      x1 <- seq(min(data.bind$x1), max(data.bind$x1), len = m)
      x2 <- seq(min(data.bind$x2), max(data.bind$x2), len = m)
      x.new <- expand.grid(x1 = x1, x2 = x2)
      fit.new <- matrix(predict(fit, x.new), m, m)
      persp(x1, x2, fit.new, theta = 40, phi = 30, ticktype = "detailed",
            xlab = "x1", ylab = "x2", zlab = "y", col = "lightblue",
            expand = 0.6)
    }
  }
  return(fit)
}


# function to perform loess normalization on a hic.table called from
# within hic_loess main function
.loess.matrix <- function(hic.table, degree = 1, Plot = FALSE, span = NA,
                          loess.criterion = "gcv") {
  # perform loess on data
  if (is.na(span)) {
    l <- .loess.as(x = hic.table$D, y = hic.table$M, degree = degree,
                   criterion = loess.criterion,
                   control = loess.control(surface = "interpolate",
                                           statistics = "approximate", trace.hat = "approximate"))
    # calculate gcv and AIC
    traceL <- l$trace.hat
    sigma2 <- sum(l$residuals^2)/(l$n - 1)
    aicc <- log(sigma2) + 1 + 2 * (2 * (traceL + 1))/(l$n - traceL -
                                                        2)
    gcv <- l$n * sigma2/(l$n - traceL)^2
  } else {
    l <- loess(hic.table$M ~ hic.table$D, degree = degree, span = span,
               control = loess.control(surface = "interpolate", statistics = "approximate",
                                       trace.hat = "approximate"))
    # calculate gcv and AIC
    traceL <- l$trace.hat
    sigma2 <- sum(l$residuals^2)/(l$n - 1)
    aicc <- log(sigma2) + 1 + 2 * (2 * (traceL + 1))/(l$n - traceL -
                                                        2)
    gcv <- l$n * sigma2/(l$n - traceL)^2
  }
  # print the span picked by gcv
  message(paste("Span for loess: ", l$pars$span, sep = ""))
  message(paste("GCV for loess: ", gcv, sep = ""))
  message(paste("AIC for loess: ", aicc, sep = ""))
  # get the correction factor
  mc <- predict(l, hic.table$D)
  if (Plot) {
    MD.plot1(hic.table$M, hic.table$D, mc)
  }
  # create mhat matrix using mc/2 which will be subtracted/added to the
  # original matrices to produce the loess normalized matrices
  mhat <- mc/2
  # normalize original matrices
  if (sum(hic.table$IF1 == 0) + sum(hic.table$IF2 == 0) > 0) {
    hic.table[, `:=`(adj.IF1, 2^(log2(IF1 + 1) + mhat))]
    hic.table[, `:=`(adj.IF2, 2^(log2(IF2 + 1) - mhat))]
    hic.table[, `:=`(adj.M, log2((adj.IF2)/(adj.IF1)))]
    hic.table[, `:=`(mc, mc)]
  } else {
    hic.table[, `:=`(adj.IF1, 2^(log2(IF1) + mhat))]
    hic.table[, `:=`(adj.IF2, 2^(log2(IF2) - mhat))]
    hic.table[, `:=`(adj.M, log2((adj.IF2)/(adj.IF1)))]
    hic.table[, `:=`(mc, mc)]
  }
  # check for negative values in normalized matrices
  if (sum(hic.table$adj.IF1 < 0, na.rm = TRUE) > 0) {
    message("negatives introduced in normalization")
  }
  if (sum(c(hic.table$adj.IF2) < 0, na.rm = TRUE) > 0) {
    message("negatives introduced in normalization")
  }
  return(hic.table)
}


# Permutation test function called from hic_loess function or hic_diff
# function
.perm.test <- function(data, iterations) {
  n <- length(data)
  numerator <- sapply(data, function(x) {
    # to ignore any NAs in the data
    if (!is.finite(x))
      return(NA) else {
        perm.data <- sample(data, size = iterations, replace = TRUE)
        test.stat <- ifelse(abs(perm.data) >= abs(x), 1, 0)
        return(sum(test.stat, na.rm = TRUE))
      }
  })
  p.value <- (numerator + 1)/(iterations + 1)
  return(p.value)
}


# function to calculate a difference threshold based on the
# distribution of M will produce a difference threshold of 2 * SD(M)
.calc.diff.thresh <- function(hic.table) {
  sd_M <- sd(hic.table$adj.M)
  diff.thresh <- 2 * sd_M
  return(diff.thresh)
}

# Fucntion to calculate p-values based on distance Called from within
# hic_loess or hic_diff functions uses perm.test function
.calc.pval <- function(hic.table, diff.thresh = NA, p.adj.method = "fdr",
                       Plot = TRUE, iterations = 10000) {
  temp <- list()
  for (dist in 0:ceiling(0.85 * max(hic.table$D))) {
    temp[[dist + 1]] <- subset(hic.table, D == dist)
    p.temp <- .perm.test(temp[[dist + 1]]$adj.M, iterations = iterations)
    temp[[dist + 1]][, `:=`(p.value, p.temp)]
    temp[[dist + 1]][, `:=`(p.adj, p.adjust(p.temp, method = p.adj.method))]
    # method to check for significant calls when the actual difference
    # between the two values is very small
    if (!is.na(diff.thresh)) {
      # M specifies the log2 fold change between IF1 and IF2. Want to call
      # differences less than user set diff.thresh fold change not clinically
      # significant
      temp[[dist + 1]][, `:=`(p.value, ifelse(p.value < 0.05 & abs(adj.M) <
                                                diff.thresh, 0.5, p.value))]
    }
  }
  # for permutation to work need to combine top distances together into
  # one group
  temp[[dist + 2]] <- subset(hic.table, D > dist)
  p.temp <- .perm.test(temp[[dist + 2]]$adj.M, iterations = iterations)
  temp[[dist + 2]][, `:=`(p.value, p.temp)]
  temp[[dist + 2]][, `:=`(p.adj, p.adjust(p.temp, method = p.adj.method))]
  # method to check for significant calls when the actual difference
  # between the two values is very small
  if (!is.na(diff.thresh)) {
    temp[[dist + 2]][, `:=`(p.value, ifelse(p.value < 0.05 & abs(adj.M) <
                                              diff.thresh, 0.5, p.value))]
  }
  hic.table <- rbindlist(temp)
  ## Dont need p.adj so remove it from hic.table before returning
  hic.table[, `:=`(p.adj, NULL)]
  # add fold change column
  hic.table[, `:=`(fold.change, adj.IF2/adj.IF1)]
  if (Plot) {
    mdplot <- MD.plot2(M = hic.table$adj.M, D = hic.table$D, p.val = hic.table$p.value,
                       diff.thresh = diff.thresh)
    print(mdplot)
    return(hic.table)
  } else {
    return(hic.table)
  }
}
