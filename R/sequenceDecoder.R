#' One Hot Decoder from One Hot Encoded Matrix or 3D Array
#' 
#' Use this to transform one hot encoded sequences back
#' into amino acid or nucleotide sequences.
#' 
#' @examples
#' new.sequences <- generateSequences(prefix.motif = "CAS",
#'                                    suffix.motif = "YF",
#'                                    number.of.sequences = 100,
#'                                    min.length = 8,
#'                                    max.length = 16)
#'                           
#' sequence.matrix <- onehotEncoder(new.sequences, 
#'                                  convert.to.matrix = TRUE)
#'                                  
#' decoded.sequences <- sequenceDecoder(sequence.matrix,
#'                                      padding.symbol = ".")
#' 
#' @param sequence.matrix The encoded sequences to decode in an array or matrix
#' @param encoder The method to prepare the sequencing information - 
#' "onehotEncoder" or "propertyEncoder"
#' @param aa.method.to.use The method or approach to use for the conversion:
#' \itemize{
#'   \item{Individual sets: atchleyFactors, crucianiProperties, FASGAI, kideraFactors, MSWHIM,
#'   ProtFP, stScales, tScales, VHSE, zScales"}
#'   \item{Multiple Sets: c("atchleyFactors", "VHSE") }
#' } 
#' @param call.threshold The relative strictness of sequence calling with higher values being more
#' stringent
#' @param sequence.dictionary The letters to use in sequence generation 
#' (default are all amino acids)
#' @param padding.symbol Symbol to use for padding at the end of sequences
#' @param remove.padding Remove the additional symbol from the end of decoded sequences
#' @importFrom utils data
#' 
#' @export
#' @return Decoded amino acid or nucleotide sequences
sequenceDecoder <- function(sequence.matrix,
                            encoder = "onehotEncoder",
                            aa.method.to.use = NULL,
                            call.threshold = 0.5,
                            sequence.dictionary = amino.acids,
                            padding.symbol = ".", 
                            remove.padding = TRUE) {
  if(call.threshold <= 0) {
    stop("Please select number > 0 for call.threshold")
  }
  
  if(encoder %!in% c("onehotEncoder", "propertyEncoder")) {
    stop("Invalid encoder provided, please select either 'onehotEncoder' or 'propertyEncoder'.")
  }
  if(encoder == "onehotEncoder") {
    decoded_sequences <- .onehotDecoder(sequence.matrix,
                                        sequence.dictionary,
                                        padding.symbol,
                                        call.threshold)
    
  } else if (encoder == "propertyEncoder") {
    data("immapex_AA.data", package = "immApex", envir = environment())
    if(any(aa.method.to.use %!in% names(immapex_AA.data))) {
      stop(paste0("Please select one of the following for aa.method.to.use: ", paste(sort(names(immapex_AA.data)), collapse = ", ")))
    }
    decoded_sequences <- .propertyDecoder(sequence.matrix,
                                          aa.method.to.use,
                                          padding.symbol,
                                          call.threshold)
  }
  if(remove.padding) {
    remove_repetitive_end <- function(x) {
      gsub(paste0("(\\", padding.symbol, "*)$"), "", x)
    }
    decoded_sequences <- decoded_sequences <- vapply(decoded_sequences, 
                                                     remove_repetitive_end, 
                                                     FUN.VALUE = character(1), 
                                                     USE.NAMES = FALSE)
    
  }
  
  return(decoded_sequences)
}

.euclidean.distance <- function(vec1, vec2) {
  sqrt(sum((vec1 - vec2)^2))
}

#' @importFrom reticulate array_reshape
#' @importFrom utils data
.propertyDecoder <- function(sequence.matrix,
                             aa.method.to.use,
                             padding.symbol,
                             call.threshold) {
  
  data("immapex_AA.data", package = "immApex", envir = environment())
  call.threshold <- 1/call.threshold
  vectors <- immapex_AA.data[aa.method.to.use]
  vector.names <- as.vector(unlist(lapply(vectors, names)))
  vectors <- do.call(c, vectors)
  names(vectors) <- vector.names
  vectors <- lapply(vectors, .min.max.normalize)
  
  
  if (inherits(sequence.matrix, "matrix")) {
    num_sequences <- nrow(sequence.matrix)
    sequence_length <- ncol(sequence.matrix) / length(vectors)
    sequence.matrix <- array_reshape(sequence.matrix, 
                                     c(num_sequences, sequence_length, length(vectors)))
  } else {
    num_sequences <- dim(sequence.matrix)[1]
    sequence_length <- dim(sequence.matrix)[2]
  }
  vectors <- do.call(rbind, vectors)
  vectors <- cbind(vectors, c(rep(0, nrow(vectors))))
  colnames(vectors)[21] <- padding.symbol
  decoded_sequences <- character(num_sequences)
  
  for (i in seq_len(num_sequences)) {
    sequence <- ""
    for (j in seq_len(sequence_length)) {
      distances <- apply(vectors, 2, function(col) .euclidean.distance(sequence.matrix[i, j, ], col))
      if(min(distances) < call.threshold) {
        index <- names(sort(distances)[1])
        sequence <- paste0(sequence, index)
      } else {
        sequence <- paste0(sequence, padding.symbol)
      }
    }
    decoded_sequences[i] <- sequence
  }
  return(decoded_sequences)
}

#TODO Add testthat 

#' @importFrom reticulate array_reshape
.onehotDecoder <- function(sequence.matrix,
                           sequence.dictionary,
                           padding.symbol,
                           call.threshold) {
  if(call.threshold > 1) {
    call.threshold <- 1
  }
  
  if (inherits(sequence.matrix, "matrix")) {
    num_sequences <- nrow(sequence.matrix)
    sequence_length <- ncol(sequence.matrix) / (length(sequence.dictionary) + 1)
    sequence.matrix <- array_reshape(sequence.matrix, 
                                     c(num_sequences, sequence_length, (length(sequence.dictionary) + 1)))
  } else {
    num_sequences <- dim(sequence.matrix)[1]
    sequence_length <- dim(sequence.matrix)[2]
  }
  
  char_set <- c(sequence.dictionary, padding.symbol)
  
  decoded_sequences <- character(num_sequences)
  
  
  for (i in seq_len(num_sequences)) {
    sequence <- ""
    for (j in seq_len(sequence_length)) {
      index <- which(sequence.matrix[i, j, ] == max(sequence.matrix[i, j, ]))
      if (length(index) == 1 & max(sequence.matrix[i, j, ]) >= call.threshold) {
        sequence <- paste0(sequence, char_set[index])
      } else {
        sequence <- paste0(sequence, padding.symbol)
      }
    }
    decoded_sequences[i] <- sequence
  }
  return(decoded_sequences)
}

