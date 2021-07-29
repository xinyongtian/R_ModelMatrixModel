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
poly.coef = function(form, data) {
  f = as.character(form)[2]
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
model.matrix.model = function(form, data, sparse = T, center = FALSE, scale = FALSE,
                              remove_1st_dummy = F, remove_same_value_column = F) {
  f = as.character(form)[2]
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
    print(n)
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
  poly.coef = poly.coef(form, data = data)
  mx = list(form = form, x = x, poly.coef = poly.coef, scale.attr = scale.attr,
            center.attr = center.attr, x.colnames = colnames(x), sparse = sparse, factor.levelses = factor.levelses,
            remove_1st_dummy = remove_1st_dummy)
  class(mx) = "model.matrix.model"
  return(mx)
}
predict.model.matrix.model = function(obj, data, handleInvalid = "keep") {
  f = tail(as.character(obj$form), 1)
  f2 = trim(unlist(strsplit(f, "\\+")))
  f2 = f2[f2 != ""]
  fdf = data.frame(nm = f2, bt = (1:length(f2)))
  fdf$nm2 = gsub("(poly\\((.*?),.*)\\)", "\\1,coefs=obj$poly.coef$\\2)", fdf$nm)
  nbt = length(f2)
  for (nn in names(data)) {
    if (grepl(paste0("\\b", nn, "\\b"), f)) {
      if (class(data[[nn]]) == "character")
        data[[nn]] = as.factor(data[[nn]])
      if (class(data[[nn]]) == "factor") {
        if (length(setdiff(levels(data[[nn]]), obj$factor.levelses[[nn]])) !=
            0 & handleInvalid != "keep") {
          stop("invalid level(s): \"", setdiff(levels(data[[nn]]), obj$factor.levelses[[nn]]),
               "\", in column   \"", nn, "\" in test data  but not training data\n to avoid this error set handleInvalid=\"keep\"")
        }
        else if (length(levels(data[[nn]])) == 1 | (obj$remove_1st_dummy &
                                                    !isTRUE(all.equal(obj$factor.levelses[[nn]], levels(data[[nn]]))))) {
          data[[nn]] = factor(data[[nn]], levels = union(obj$factor.levelses[[nn]],
                                                         levels(data[[nn]])))
        }
      }
    }
  }
  for (n in 1:nbt) {
    print(n)
    fnc = paste(fdf[fdf$bt == n, "nm2"], collapse = "+")
    fn = formula(paste("~0+", fnc))
    xn = model.matrix(fn, data = data)
    colnames(xn) = gsub("(poly\\(.*),.+\\)", "\\1)", colnames(xn))
    colnames(xn) = gsub("[^[:alnum:]]", "_", sub(":", "_X_", colnames(xn)))
    if (obj$sparse == T) {
      xn = Matrix::Matrix(xn, sparse = obj$sparse)
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
  if (is.null(obj$center.attr))
    center = F
  else center = obj$center.attr[colnames(x)]
  if (is.null(obj$scale.attr))
    scale = F
  else scale = obj$scale.attr[colnames(x)]
  x = scale(x, center = center, scale = scale)
  if (obj$sparse == T)
    x = Anysubset(x, obj$x.colnames)
  else x = anysubset(x, obj$x.colnames)
  obj$x = x
  return(obj)
}
