trim <- function(x) gsub("^\\s+|\\s+$", "", x)
rm.1level.col = function(mt) mt[, apply(mt, 2, function(x) length(unique(x)) != 1),
                                drop = F]
anysubset = function(mt, mt2col, fill = 0) {
  diffcol = setdiff(mt2col, dimnames(mt)[[2]])
  diffmt = matrix(fill, nrow(mt), length(diffcol), dimnames = list(NULL, diffcol))
  mt2 = cbind(mt, diffmt)[, mt2col]
  return(mt2)
}
Anysubset = function(mt, mt2col, fill = 0) {
  diffcol = setdiff(mt2col, dimnames(mt)[[2]])
  diffmt = Matrix::Matrix(fill, nrow(mt), length(diffcol), dimnames = list(NULL, diffcol))
  mt2 = cbind(mt, diffmt)[, mt2col]
  return(mt2)
}
poly.coef = function(rformula, data) {
  f = as.character(rformula)[2]
  fp = regmatches(f, gregexpr("poly\\([^)]+\\)", f))
  fpu = unique(fp[[1]])
  fpo = sub("poly\\((.*?),.*", "\\1", fpu)
  getpolycoef2 = function(var, df = data) {
    txt = paste("with (df, ", var, ")", sep = "")
    a = epp(txt)
    ac = attr(a, "coefs")
    return(ac)
  }
  poly.coef = lapply(fpu, getpolycoef2)
  names(poly.coef) = fpo
  return(poly.coef)
}
epp = function(x) eval.parent(parse(text = x))

#' Create model.matrix and save the fitted parameters
#'
#' This function create model.matrix and save the fitted parameters
#'
#' @param rformula a formula, e.g. output of formula()
#' @param data a data.frame
#' @param sparse bool, if True return a sparse matrix.
#' @param center bool, if center the output
#' @param scale  bool, if scale the output
#' @param remove_1st_dummy bool, if remove the first dummy variable in one hot key transformation
#' @param remove_same_value_column bool, if remove column with one unique value
#' @return A ModelMatrixModel class,which includes the transformed matrix and  the fitted parameters.
#' @export
ModelMatrixModel = function(rformula, data, sparse = T, center = FALSE, scale = FALSE,
                              remove_1st_dummy = F, remove_same_value_column = F,verbose=T) {
  f = as.character(rformula)[2]
  f2 = trim(unlist(strsplit(f, "\\+")))
  f2 = f2[f2 != ""]
  fdf = data.frame(nm = f2, bt = (1:length(f2)))
  nbt = length(f2)
  center.attr = NULL
  scale.attr = NULL
  factor.levelses = NULL
  for (nm in names(data)) {
    if (grepl(nm, f)) {
      if (class(data[[nm]]) == "character")
        data[[nm]] = as.factor(data[[nm]])
      if (class(data[[nm]]) == "factor") {
        le = list(levels(data[[nm]]))
        names(le) = nm
        factor.levelses = c(factor.levelses, le)
      }
    }
  }
  for (n in 1:nbt) {
    if (verbose)print(paste(n, ' out of ',nbt))
    fnc = paste(fdf[fdf$bt == n, "nm"], collapse = "+")
    if (remove_1st_dummy) {
      fn = formula(paste("~", fnc))
      xn = model.matrix(fn, data = data)
      xn = xn[, -1, drop = F]
    }
    else {
      fn = formula(paste("~0+", fnc))
      xn = model.matrix(fn, data = data)
    }
    if (remove_same_value_column)
      xn = rm.1level.col(xn)
    xn = scale(xn, center = center, scale = scale)
    if (center == T) {
      xns = attr(xn, "scaled:center")
      center.attr = c(center.attr, xns)
    }
    if (scale == T) {
      xns = attr(xn, "scaled:scale")
      scale.attr = c(scale.attr, xns)
    }
    if (sparse == T) {
      xn = Matrix::Matrix(xn, sparse = sparse)
      if (n == 1) {
        x = xn
      }
      else {
        x = cbind(x, xn)
      }
    }
    else {
      if (n == 1) {
        x = xn
      }
      else {
        x = cbind(x, xn)
      }
    }
  }
  colnames(x) = gsub("[^[:alnum:]]", "_", sub(":", "_X_", colnames(x)))
  if (center == T)
    names(center.attr) = gsub("[^[:alnum:]]", "_", names(center.attr))
  if (scale == T)
    names(scale.attr) = gsub("[^[:alnum:]]", "_", names(scale.attr))
  poly.coef = poly.coef(rformula, data = data)
  mx = list(rformula = rformula, x = x, poly.coef = poly.coef, scale.attr = scale.attr,
            center.attr = center.attr, x.colnames = colnames(x), sparse = sparse, factor.levelses = factor.levelses,
            remove_1st_dummy = remove_1st_dummy)
  class(mx) = "ModelMatrixModel"
  return(mx)
}



#' Create model.matrix and save the fitted parameters
#'
#' This function transform new data based on information from a ModelMatrixModel object
#'
#' @param object a ModelMatrixModel object
#' @param data a data.frame
#' @param handleInvalid a string,'keep' or 'error'.  In dummy variable transformation, if categorical variable has a factor level that is unseen before, 'keep' will keep the record, output dummy variables will be all zero.
#' @param ... other parameters
#' @return A ModelMatrixModel class,which includes the transformed matrix and  the fitted parameters copied from input object.
#' @export
predict.ModelMatrixModel = function(object, data, handleInvalid = "keep",...) {
  f = tail(as.character(object$rformula), 1)
  f2 = trim(unlist(strsplit(f, "\\+")))
  f2 = f2[f2 != ""]
  fdf = data.frame(nm = f2, bt = (1:length(f2)))
  fdf$nm2 = gsub("(poly\\((.*?),.*)\\)", "\\1,coefs=object$poly.coef$\\2)", fdf$nm)
  nbt = length(f2)
  for (nn in names(data)) {
    if (grepl(paste0("\\b", nn, "\\b"), f)) {
      if (class(data[[nn]]) == "character")
        data[[nn]] = as.factor(data[[nn]])
      if (class(data[[nn]]) == "factor") {
        if (length(setdiff(levels(data[[nn]]), object$factor.levelses[[nn]])) !=
            0 & handleInvalid != "keep") {
          stop("invalid level(s): \"", setdiff(levels(data[[nn]]), object$factor.levelses[[nn]]),
               "\", in column   \"", nn, "\" in test data  but not training data\n to avoid this error set handleInvalid=\"keep\"")
        }
        else if (length(levels(data[[nn]])) == 1 | (object$remove_1st_dummy &
                                                    !isTRUE(all.equal(object$factor.levelses[[nn]], levels(data[[nn]]))))) {
          data[[nn]] = factor(data[[nn]], levels = union(object$factor.levelses[[nn]],
                                                         levels(data[[nn]])))
        }
      }
    }
  }
  for (n in 1:nbt) {
    if (verbose)print(paste(n, ' out of ',nbt))
    fnc = paste(fdf[fdf$bt == n, "nm2"], collapse = "+")
    fn = formula(paste("~0+", fnc))
    xn = model.matrix(fn, data = data)
    colnames(xn) = gsub("(poly\\(.*),.+\\)", "\\1)", colnames(xn))
    colnames(xn) = gsub("[^[:alnum:]]", "_", sub(":", "_X_", colnames(xn)))
    if (object$sparse == T) {
      xn = Matrix::Matrix(xn, sparse = object$sparse)
      if (n == 1) {
        x = xn
      }
      else {
        x = cbind(x, xn)
      }
    }
    else {
      if (n == 1) {
        x = xn
      }
      else {
        x = cbind(x, xn)
      }
    }
  }
  if (is.null(object$center.attr))
    center = F
  else center = object$center.attr[colnames(x)]
  if (is.null(object$scale.attr))
    scale = F
  else scale = object$scale.attr[colnames(x)]
  x = scale(x, center = center, scale = scale)
  if (object$sparse == T)
    x = Anysubset(x, object$x.colnames)
  else x = anysubset(x, object$x.colnames)
  object$x = x
  return(object)
}


