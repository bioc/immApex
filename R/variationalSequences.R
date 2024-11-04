#' Generate Similar Sequences using Variational Autoencoder
#' 
#' Use this to simulate sequences using a variational 
#' autoencoder (VAE) and perturbation of the probability distributions.
#' 
#' @examples
#' \dontrun{
#' sequences <- generateSequences(prefix.motif = "CAS",
#'                                suffix.motif = "YF",
#'                                number.of.sequences = 100,
#'                                min.length = 8,
#'                                max.length = 16)
#'                                 
#' new.sequences <- variationalSequences(sequences, 
#'                                       encoder = "onehotEncoder",
#'                                       encoder.hidden.dim = c(256, 128),
#'                                       latent.dim = 16,
#'                                       batch.size = 16)
#' } 
#' 
#' @param input.sequences The amino acid or nucleotide sequences to use
#' @param encoder.function The method to prepare the sequencing information - 
#' "onehotEncoder" or "propertyEncoder"
#' @param aa.method.to.use The method or approach to use for the conversion:
#' \itemize{
#'   \item{Individual sets: atchleyFactors, crucianiProperties, FASGAI, kideraFactors, MSWHIM,
#'   ProtFP, stScales, tScales, VHSE, zScales"}
#'   \item{Multiple Sets: c("atchleyFactors", "VHSE") }
#' } 
#' @param number.of.sequences Number of sequences to generate
#' @param encoder.hidden.dim A vector of the neurons to use in the hidden layers
#' for the encoder portion of the model
#' @param decoder.hidden.dim A vector of the neurons to use in the hidden layers 
#' for the decoder portion of the model. If NULL assumes symmetric autoencoder
#' @param latent.dim The size of the latent dimensions
#' @param batch.size The batch size to use for VAE training
#' @param epochs The number of epochs to use in VAE training
#' @param learning.rate The learning rate to use in VAE training
#' @param epsilon.std The epsilon to use in VAE training
#' @param call.threshold The relative strictness of sequence calling
#'  with higher values being more stringent
#' @param activation.function The activation for the dense connected layers
#' @param optimizer The optimizer to use in VAE training
#' @param disable.eager.execution Disable the eager execution parameter for
#' tensorflow.
#' @param sequence.dictionary The letters to use in sequence mutation
#' (default are all amino acids)
#' @param verbose Print messages corresponding to the processing step
#' 
#' @importFrom keras3 layer_dense layer_lambda  keras_model compile 
#' fit layer_input loss_binary_crossentropy optimizer_adadelta 
#' optimizer_adagrad optimizer_adam optimizer_adamax optimizer_ftrl 
#' optimizer_nadam optimizer_rmsprop optimizer_sgd  
#' callback_early_stopping layer_normalization 
#' @importFrom magrittr %>%
#' @importFrom stats predict runif
#' @importFrom tensorflow tf
#' @export 
#' @return A vector of mutated sequences

variationalSequences <- function(input.sequences,
                                 encoder.function = "onehotEncoder",
                                 aa.method.to.use = NULL,
                                 number.of.sequences = 100,
                                 encoder.hidden.dim = c(128,64),
                                 decoder.hidden.dim = NULL,
                                 latent.dim = 16,
                                 batch.size = 16,
                                 epochs = 50,
                                 learning.rate = 0.001,
                                 epsilon.std = 1,
                                 call.threshold = 0.2,
                                 activation.function = "relu",
                                 optimizer = "adam",
                                 disable.eager.execution = FALSE,
                                 sequence.dictionary = amino.acids,
                                 verbose = TRUE) {
  
  n_train <- floor(length(input.sequences) * 0.8)  # Default to 80% for training
  
  # Input validation
  if(length(input.sequences) < 1) stop("input.sequences must have at least one sequence.")

  
  if (disable.eager.execution) {
    tensorflow::tf$compat$v1$disable_eager_execution()
  }
  
  es <- keras3::callback_early_stopping(
    monitor = "val_loss",
    min_delta = 0,
    patience = epochs/5,
    verbose = 0,
    mode = "min")
  
  optimizer.to.use <- switch(optimizer,
                             "adadelta" = optimizer_adadelta(learning_rate = learning.rate),
                             "adagrad" = optimizer_adagrad(learning_rate = learning.rate),
                             "adam" = optimizer_adam(learning_rate = learning.rate),
                             "adamax" = optimizer_adamax(learning_rate = learning.rate),
                             "ftrl" = optimizer_ftrl(learning_rate = learning.rate),
                             "nadam" = optimizer_nadam(learning_rate = learning.rate),
                             "rmsprop" = optimizer_rmsprop(learning_rate = learning.rate),
                             "sgd" = optimizer_sgd(learning_rate = learning.rate),
                             stop("Please select a compatible optimizer function in the Keras R implementation."))
  
  if(verbose) {
    message("Converting to matrix....")
  }
  # Prepare the sequences matrix
  sequence.matrix <- switch(encoder.function,
                            "onehotEncoder" = onehotEncoder(input.sequences, 
                                                            sequence.dictionary = sequence.dictionary,
                                                            convert.to.matrix = TRUE),
                            "propertyEncoder" = propertyEncoder(input.sequences, 
                                                                method.to.use = aa.method.to.use,
                                                                convert.to.matrix = TRUE),
                            stop("Invalid encoder provided."))
  
  # Custom VAE Loss Layer
  vae_loss_layer <- function(original_dim) {
        layer_lambda(f = function(x) {
          x_decoded_mean <- x[[1]]
          x_input <- x[[2]]
          z_mean <- x[[3]]
          z_log_var <- x[[4]]
          xent_loss <- loss_binary_crossentropy(x_input, x_decoded_mean) * original_dim
          kl_loss <- -0.5 * tf$reduce_mean(1 + z_log_var - tf$square(z_mean) - tf$exp(z_log_var), axis = -1L)
          tf$reduce_mean(xent_loss + kl_loss)
        })
  }
  original_dim <- ncol(sequence.matrix)
  
  # Data splitting
  train_indices <- sample(seq_len(nrow(sequence.matrix)), n_train)
  x_train <- sequence.matrix[train_indices, ]
  x_test <- sequence.matrix[-train_indices, ]
  
  
  # Encoder
  encoder_input <- layer_input(shape = original_dim)
  h <- encoder_input
  for (dim in encoder.hidden.dim) {
       h <- layer_dense(h, units = dim, activation = activation.function)
  }
  z_mean <- layer_dense(h, units = latent.dim, name = "z_mean")
  z_log_var <- layer_dense(h, units = latent.dim, name = "z_log_var")
      
  # Sampling Layer
  z <- layer_lambda(f = function(args) {
        z_mean <- args[[1]]
        z_log_var <- args[[2]]
        batch <- tf$shape(z_mean)[1]
        dim <- tf$shape(z_mean)[2]
        epsilon <- tf$random$normal(shape = c(batch, dim), mean = 0., stddev = epsilon.std)
        z_mean + tf$exp(z_log_var / 2) * epsilon
      }, output_shape = c(latent.dim))(list(z_mean, z_log_var))
      
  # Decoder
  decoder_input <- layer_input(shape = latent.dim)
  d <- decoder_input
  if (is.null(decoder.hidden.dim)) {
    decoder.hidden.dim <- rev(encoder.hidden.dim)  # Default to mirroring the encoder layers
  }
  for (dim in decoder.hidden.dim) {
    d <- layer_dense(d, units = dim, activation = "relu")
  }
  decoder_output <- layer_dense(d, units = original_dim, activation = "sigmoid")
      
  # Encoder and Decoder Models
  encoder <- keras_model(encoder_input, z_mean)
  decoder <- keras_model(decoder_input, decoder_output)
      
  # VAE Model
  decoder_output <- decoder(z)
  vae <- keras_model(encoder_input, decoder_output)
      
  # Add custom loss layer
  loss_layer <- vae_loss_layer(original_dim)(list(decoder_output, encoder_input, z_mean, z_log_var))
  vae_with_loss <- keras_model(encoder_input, loss_layer)
      
  # Dummy loss function
  dummy_loss <- function(y_true, y_pred) {
    tf$reduce_mean(y_pred)
  }
      
  # Compile the model
  vae_with_loss %>% compile(optimizer = optimizer.to.use, loss = dummy_loss)
  
  if(verbose) {    
    message("Fitting Model....")
  }
  vae_with_loss %>% fit(
        x_train, x_train, 
        shuffle = TRUE,
        epochs = epochs,
        batch_size = batch.size,
        validation_data = list(x_test, x_test),
        verbose = 0,
        callbacks = es
  )
  if(verbose) {
    message("Generating New Sequences....")
  }
  encoded_sequences <- as.matrix(encoder(x_train))
  
  #Using the vectors/ranges of training sequences to form a new matrix
  lapply(seq_len(ncol(encoded_sequences)), function(x) {
    runif(number.of.sequences, min = min(encoded_sequences[,x]), max = max(encoded_sequences[,x]))
  }) -> z_sample
  
  z_sample <- do.call(cbind, z_sample)
  generated_matrix <- predict(decoder, z_sample)
        
  candidate.sequences <- sequenceDecoder(generated_matrix,
                                         encoder = encoder.function,
                                         aa.method.to.use = aa.method.to.use,
                                         call.threshold = call.threshold,
                                         sequence.dictionary = sequence.dictionary,
                                         padding.symbol = ".")
  return(candidate.sequences)
}
