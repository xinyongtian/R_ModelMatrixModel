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
extract_poly_coef = function(rformula, data) {
  rformula_str = as.character(rformula)[2]
  fp = regmatches(rformula_str, gregexpr("poly\\([^)]+\\)", rformula_str))
  fpu = unique(fp[[1]])
  fpo = sub("poly\\((.*?),.*", "\\1", fpu)
  getpolycoef2 = function(var, df = data) {
    txt = paste("with (df, ", var, ")", sep = "")
    a = eval_parent(txt)
    ac = attr(a, "coefs")
    return(ac)
  }
  extract_poly_coef = lapply(fpu, getpolycoef2)
  names(extract_poly_coef) = fpo
  return(extract_poly_coef)
}
eval_parent = function(x) eval.parent(parse(text = x))

#' ModelMatrixModel() function
#'
#' This function  transforms a data.frame to matrix  based on a r formula. It is similar to model.matrix() function. But it output a class stored with the transformed matrix and the transforming parameters.
#'
#' @param rformula a formula, e.g. rformula=formula("~ 1+x1+x2").Note the interpreting of the formula might be different slightly from model.matrix function. In model.matrix(),intercept column will be included in output matrix with or without "1" in the formula. But in ModelMatrixModel(),intercept column will  be included in output matrix only when  "1" is present. Moreover "0" in the formula will be ignored.
#' @param data a data.frame
#' @param sparse bool, if True return a sparse matrix, i.e. a
#' @param center bool, if center the output
#' @param scale  bool, if scale the output
#' @param remove_1st_dummy bool, if remove the first dummy variable in one hot key transformation
#' @param verbose bool, if print out progress
#' @return A ModelMatrixModel class,which includes the transformed matrix and  the transforming parameters.
#' @export
#' @example
#' see vignettes
ModelMatrixModel = function(rformula, data, sparse = T, center = F, scale = F,
                              remove_1st_dummy = F,verbose=F) {
  rformula_str = as.character(rformula)[2]
  rformula_items = trim(unlist(strsplit(rformula_str, "[+-]")))
  rformula_items =rformula_items[!rformula_items%in%c("0","")]
  rformula_df = data.frame(item_name = rformula_items, item_number = (1:length(rformula_items)))
  rformula_items_length = length(rformula_items)
  center.attr = NULL
  scale.attr = NULL
  factor.levelses = NULL
  for (item_name in names(data)) {
    if (grepl(item_name, rformula_str)) {
      if (class(data[[item_name]]) == "character")
        data[[item_name]] = as.factor(data[[item_name]])
      if (class(data[[item_name]]) == "factor") {
        le = list(levels(data[[item_name]]))
        names(le) = item_name
        factor.levelses = c(factor.levelses, le)
      }
    }
  }
  for (n in 1:rformula_items_length) {
    if (verbose)print(paste(n, ' out of ',rformula_items_length))
    nth_item_formula_str =rformula_df[rformula_df$item_number == n, "item_name"]
    if (remove_1st_dummy) {
      nth_item_formula = formula(paste("~", nth_item_formula_str))
      nth_item_modelmatrix = model.matrix(nth_item_formula, data = data)
      if (nth_item_formula_str!="1") nth_item_modelmatrix = nth_item_modelmatrix[, -1, drop = F]
    }
    else {
      nth_item_formula = formula(paste("~0+", nth_item_formula_str))
      nth_item_modelmatrix = model.matrix(nth_item_formula, data = data)
    }
    nth_item_modelmatrix=scale(nth_item_modelmatrix, center = center, scale = scale)
    if (center == T) {
      xns = attr(nth_item_modelmatrix, "scaled:center")
      center.attr = c(center.attr, xns)
    }
    if (scale == T) {
      xns = attr(nth_item_modelmatrix, "scaled:scale")
      scale.attr = c(scale.attr, xns)
    }
    if (sparse == T) {
      nth_item_modelmatrix = Matrix::Matrix(nth_item_modelmatrix, sparse = sparse)
      if (n == 1) {
        x = nth_item_modelmatrix
      }
      else {
        x = cbind(x, nth_item_modelmatrix)
      }
    }
    else {
      if (n == 1) {
        x = nth_item_modelmatrix
      }
      else {
        x = cbind(x, nth_item_modelmatrix)
      }
    }
  }
  colnames(x) = gsub("[^[:alnum:]]", "_", sub(":", "_X_", colnames(x)))
  if (center == T)
    names(center.attr) = gsub("[^[:alnum:]]", "_", names(center.attr))
  if (scale == T)
    names(scale.attr) = gsub("[^[:alnum:]]", "_", names(scale.attr))
  extract_poly_coef = extract_poly_coef(rformula, data = data)
  mx = list(rformula = rformula, x = x, extract_poly_coef = extract_poly_coef, scale.attr = scale.attr,
            center.attr = center.attr, x.colnames = colnames(x), sparse = sparse, factor.levelses = factor.levelses,
            remove_1st_dummy = remove_1st_dummy)
  class(mx) = "ModelMatrixModel"
  return(mx)
}



#' predict() function
#'
#' This function transform new data based on transforming parameters from a ModelMatrixModel object
#'
#' @param object a ModelMatrixModel object
#' @param data a data.frame
#' @param handleInvalid a string,'keep' or 'error'.  In dummy variable transformation, if categorical variable has a factor level that is unseen before, 'keep' will keep the record, output dummy variables will be all zero.
#' @param ... other parameters
#' @param verbose bool, if print out progress
#' @return A ModelMatrixModel class,which includes the transformed matrix and  the necessary transforming parameters copied from input object.
#' @export
#' @examples
#' see vignettes
predict.ModelMatrixModel = function(object, data, handleInvalid = "keep",verbose=F, ...) {
  rformula_str = tail(as.character(object$rformula), 1)
  rformula_items = trim(unlist(strsplit(rformula_str, "[+-]")))
  rformula_items =rformula_items[!rformula_items%in%c("0","")]
  rformula_df = data.frame(item_name = rformula_items, item_number = (1:length(rformula_items)))
  rformula_df$item_name_modified = gsub("(poly\\((.*?),.*)\\)", "\\1,coefs=object$extract_poly_coef$\\2)", rformula_df$item_name)
  rformula_items_length = length(rformula_items)
  for (col in names(data)) {
    if (grepl(paste0("\\b", col, "\\b"), rformula_str)) {
      if (class(data[[col]]) == "character")
        data[[col]] = as.factor(data[[col]])
      if (class(data[[col]]) == "factor") {
        if (length(setdiff(levels(data[[col]]), object$factor.levelses[[col]])) !=
            0 & handleInvalid != "keep") {
          stop("invalid level(s): \"", setdiff(levels(data[[col]]), object$factor.levelses[[col]]),
               "\", in column \"", col, "\" in test data  but not training data\n to avoid this error set handleInvalid=\"keep\"")
        }
        else if (length(levels(data[[col]])) == 1 | (object$remove_1st_dummy &
                                                    !isTRUE(all.equal(object$factor.levelses[[col]], levels(data[[col]]))))) {
          data[[col]] = factor(data[[col]], levels = union(object$factor.levelses[[col]],
                                                         levels(data[[col]])))
        }
      }
    }
  }
  for (n in 1:rformula_items_length) {
    if (verbose)print(paste(n, ' out of ',rformula_items_length))
    nth_item_formula_str = rformula_df[rformula_df$item_number == n, "item_name_modified"]
    nth_item_formula = formula(paste("~0+", nth_item_formula_str))
    nth_item_modelmatrix = model.matrix(nth_item_formula, data = data)
    colnames(nth_item_modelmatrix) = gsub("(poly\\(.*),.+\\)", "\\1)", colnames(nth_item_modelmatrix))
    colnames(nth_item_modelmatrix) = gsub("[^[:alnum:]]", "_", sub(":", "_X_", colnames(nth_item_modelmatrix)))
    if (object$sparse == T) {
      nth_item_modelmatrix = Matrix::Matrix(nth_item_modelmatrix, sparse = object$sparse)
      if (n == 1) {
        x = nth_item_modelmatrix
      }
      else {
        x = cbind(x, nth_item_modelmatrix)
      }
    }
    else {
      if (n == 1) {
        x = nth_item_modelmatrix
      }
      else {
        x = cbind(x, nth_item_modelmatrix)
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


