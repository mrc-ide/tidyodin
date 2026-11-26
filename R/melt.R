## Core reshaping: everything in the package ends up here.
##
## `x` is an array whose first dimension is packed state (or packed
## parameters), with one name per row in `names`. `trailing` describes
## every other dimension, in array order, as a named list of the values
## each dimension takes. The special trailing name ".draw" is expanded
## into `.chain`, `.iteration` and `.draw` using `n_samples`.

column_order <- c(".chain", ".iteration", ".draw", "group", "particle",
                  "time")

melt_state <- function(x, names, trailing, vars = NULL, dims = NULL,
                       labels = NULL, name_col = "variable",
                       n_samples = NULL, call = rlang::caller_env()) {
  expected <- c(length(names), lengths(trailing, use.names = FALSE))
  actual <- dim(x) %||% length(x)
  if (!identical(as.integer(actual), as.integer(expected))) {
    cli::cli_abort(
      c("Unexpected array dimensions",
        x = "Expected {.val {paste(expected, collapse = ' x ')}} ({paste(c(name_col, names(trailing)), collapse = ' x ')})",
        x = "Found {.val {paste(actual, collapse = ' x ')}}"),
      call = call)
  }

  info <- parse_names(names)
  n_rest <- prod(expected[-1])

  if (!is.null(vars)) {
    unknown <- setdiff(vars, info$variable)
    if (length(unknown) > 0) {
      cli::cli_abort(
        c("Unknown {name_col}{?s} in {.arg vars}: {.val {unknown}}",
          i = "Available: {.val {unique(info$variable)}}"),
        call = call)
    }
    keep <- info$variable %in% vars
    x <- matrix(x, length(names), n_rest)[keep, , drop = FALSE]
    info <- lapply(info, function(el) {
      if (is.matrix(el)) el[keep, , drop = FALSE] else el[keep]
    })
  }

  n_state <- length(info$variable)
  ret <- list()

  ## Trailing dimensions: state varies fastest, the last dimension slowest
  n_before <- n_state
  for (i in seq_along(trailing)) {
    values <- trailing[[i]]
    n_after <- prod(lengths(trailing)[-seq_len(i)])
    ret[[names(trailing)[[i]]]] <-
      rep(rep(values, each = n_before), times = n_after)
    n_before <- n_before * length(values)
  }

  if (!is.null(ret$.draw)) {
    draw <- ret$.draw
    ret$.chain <- (draw - 1L) %/% n_samples + 1L
    ret$.iteration <- (draw - 1L) %% n_samples + 1L
  }
  if (!is.null(ret$group)) {
    ret$group <- apply_labels(ret$group, labels$group, "group", call)
  }
  ret <- ret[intersect(column_order, names(ret))]

  ret[[name_col]] <- rep(info$variable, times = n_rest)
  index <- index_columns(info, dims, labels, call)
  for (nm in names(index)) {
    ret[[nm]] <- rep(index[[nm]], times = n_rest)
  }
  if (!all(is.na(info$group))) {
    ret$group <- rep(info$group, times = n_rest)
  }
  ret$value <- as.vector(x)

  tibble::new_tibble(ret, nrow = length(ret$value))
}


## Split flat packer names into their parts:
##
##   "S"         -> variable "S"
##   "x[2,3]"    -> variable "x", index 2, 3
##   "beta<a>"   -> variable "beta", group "a"
##   "x[1]<a>"   -> variable "x", index 1, group "a"
parse_names <- function(names) {
  has_group <- grepl("<[^<>]*>$", names)
  group <- rep(NA_character_, length(names))
  group[has_group] <- sub("^.*<([^<>]*)>$", "\\1", names[has_group])
  base <- sub("<[^<>]*>$", "", names)

  has_index <- grepl("\\[[0-9, ]+\\]$", base)
  variable <- base
  variable[has_index] <- sub("\\[[0-9, ]+\\]$", "", base[has_index])
  index_str <- sub("^.*\\[([0-9, ]+)\\]$", "\\1", base[has_index])
  parts <- strsplit(index_str, ",", fixed = TRUE)
  rank <- max(0L, lengths(parts))
  index <- matrix(NA_integer_, length(names), rank)
  for (i in seq_along(parts)) {
    p <- as.integer(parts[[i]])
    index[which(has_index)[[i]], seq_along(p)] <- p
  }
  list(variable = variable, index = index, group = group)
}


## Build the index columns for array variables. `dims` names the
## dimensions of each variable and `labels` gives labels for a named
## dimension.
index_columns <- function(info, dims, labels, call) {
  rank_max <- ncol(info$index)
  if (rank_max == 0) {
    return(list())
  }
  variables <- unique(info$variable)
  rank <- vapply(variables, function(v) {
    sum(!is.na(info$index[match(v, info$variable), ]))
  }, integer(1))

  if (!is.null(dims) && !is.list(dims)) {
    if (!is.character(dims) || !is.null(names(dims))) {
      cli::cli_abort(
        "{.arg dims} must be an unnamed character vector or a named list",
        call = call)
    }
    if (length(dims) < rank_max) {
      cli::cli_abort(
        "{.arg dims} names {length(dims)} dimension{?s} but some variables have {rank_max}",
        call = call)
    }
    dims <- lapply(rank[rank > 0], function(r) dims[seq_len(r)])
  }

  ret <- list()
  for (v in variables[rank > 0]) {
    r <- rank[[v]]
    nms <- dims[[v]] %||% paste0("index", seq_len(r))
    if (length(nms) != r) {
      cli::cli_abort(
        "{.arg dims} for {.val {v}} has length {length(nms)} but {.val {v}} has {r} dimension{?s}",
        call = call)
    }
    rows <- info$variable == v
    for (j in seq_len(r)) {
      if (is.null(ret[[nms[[j]]]])) {
        ret[[nms[[j]]]] <- rep(NA_integer_, length(info$variable))
      }
      ret[[nms[[j]]]][rows] <- info$index[rows, j]
    }
  }
  for (nm in names(ret)) {
    ret[[nm]] <- apply_labels(ret[[nm]], labels[[nm]], nm, call)
  }
  ret
}


apply_labels <- function(index, labels, name, call) {
  if (is.null(labels)) {
    return(index)
  }
  if (!is.numeric(index)) {
    return(factor(index, levels = labels))
  }
  n <- max(index, na.rm = TRUE)
  if (length(labels) < n) {
    cli::cli_abort(
      "{.arg labels} for {.val {name}} has {length(labels)} element{?s} but {n} are needed",
      call = call)
  }
  factor(labels[index], levels = labels)
}


`%||%` <- function(x, y) if (is.null(x)) y else x
