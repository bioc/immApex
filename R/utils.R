amino.acids <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")

"%!in%" <- Negate("%in%")

.check.sequences <- function(sequences, sequence.dictionary) {
  if(length(sequences) >= 10) {
    any(unlist(strsplit(sequences[seq(1,10)], "")) %!in% sequence.dictionary)
  } else {
    any(unlist(strsplit(sequences[seq_len(length(sequences))], "")) %!in% sequence.dictionary)
  }
}
  


#Add additional sequence padding to max length
.padded.strings <- function(strings, 
                            max.length,
                            padded.token = NULL,
                            concatenate = TRUE) {
  max_length <- max.length
  
  #Produce combined padded strings
  if(concatenate) {
    x <- lapply(strings, function(str) {
      str_len <- nchar(str)
      if (str_len < max_length) {
        str <- paste(c(str, rep(padded.token, max_length - str_len)), collapse = "")
      } else {
        str
      }
    })
  #Produce a list of strings with seperate vectors
  } else {
    x <- lapply(strings, function(str) {
      str_len <- length(str)
      if (str_len < max_length) {
        str <- c(str, rep(padded.token, max_length - str_len))
      } else {
        str
      }
    })
  }
  return(x)
  
}

substring.extractor <- function(strings, motif.length) {
  lapply(strings, function(x) {
    string_length <- nchar(x)
    num_substrings <- string_length - motif.length + 1
    if (num_substrings > 0) {
      # Generate all substrings of the specified length
      substrings <- vapply(seq_len(num_substrings), 
                           function(j) substring(x, j, j + motif.length - 1), 
                           FUN.VALUE = character(1))
    } else {
      # Return NA if the string is too short
      substrings <- NA
    }
    substrings
  })
}

.min.max.normalize <- function(x){
  (x- min(x)) /(max(x)-min(x))
}

#' @importFrom matrixStats colMedians colMeans2 colSums2 colVars colMads
.get.stat.function <- function(method) {
  statFunc <- switch(method,
                     "median" = colMedians,
                     "mean"  = colMeans2,
                     "sum"      = colSums2,
                     "vars" = colVars,
                     "mads"  = colMads, 
                     stop("Invalid summary.function provided"))
  return(statFunc) 
}

.is_seurat_object <- function(obj) inherits(obj, "Seurat")
.is_se_object <- function(obj) inherits(obj, "SummarizedExperiment")
.is_seurat_or_se_object <- function(obj) {
  .is_seurat_object(obj) || .is_se_object(obj)
}

.get.genes.updated <- function(input.data, technology, region) {
  if (.is_seurat_or_se_object(input.data)) {
    genes.updated <- region
  } else {
    if (is.null(technology)) technology <- NA
    if (technology %in% c("TenX", "Adaptive")) {
      potential_col <- paste0(region, "_gene")
      if (!(potential_col %in% colnames(input.data))) {
        genes.updated <- paste0(region, "GeneName")
      } else {
        genes.updated <- potential_col
      }
    } else if (technology %in% c("AIRR")) {
      genes.updated <- paste0(region, "_call")
    } else {
      genes.updated <- region
    }
  }
  return(genes.updated)
}

#' @importFrom stringr str_sort
.array.dimnamer <- function(array) {
  combinations <- expand.grid(dimnames(array)[[2]], dimnames(array)[[3]], stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
  combinations[,1] <- str_sort(combinations[,1], numeric = TRUE)
  combinations[,2] <- dimnames(array)[[3]]
  combined_strings <- apply(combinations, 1, function(x) paste0(x[1], "_", x[2]))
  return(combined_strings)
}

.array.reshape <- function(x, dim, order = c("C", "F")) {
  
  order <- match.arg(order)
  dim   <- as.integer(dim)
  
  # Checks 
  if (!is.atomic(x))
    stop("`x` must be an atomic vector, matrix, or array.", call. = FALSE)
  if (any(dim <= 0))
    stop("All elements of `dim` must be positive.", call. = FALSE)
  if (length(x) != prod(dim))
    stop(sprintf("Input of length %d cannot be reshaped to %s (product = %d).",
                 length(x), paste(dim, collapse = "x"), prod(dim)),
         call. = FALSE)
  
  # row-major flatten 
  row_flatten <- function(a) {
    if (is.null(dim(a))) {
      as.vector(a)                                 # already 1D
    } else {
      as.vector(aperm(a, rev(seq_along(dim(a)))))  # reverse dims then flatten
    }
  }
  
  #  build result 
  if (order == "F") {
    array(as.vector(x), dim = dim)
    
  } else {  
    vals_C <- row_flatten(x)
    tmp    <- array(vals_C, dim = rev(dim))
    aperm(tmp, rev(seq_along(dim)))
  }
}


