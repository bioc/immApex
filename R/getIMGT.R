#' Get IMGT Sequences for Specific Loci
#'
#' Use this to access the ImMunoGeneTics (IMGT) sequences for a
#' specific species and gene loci via the
#' \href{https://github.com/BorchLab/immReferent}{immReferent} package,
#' which provides automatic local caching. More information on
#' IMGT can be found at \href{https://www.imgt.org/}{imgt.org}.
#'
#' @examples
#' \dontrun{
#' TRBV_aa <- getIMGT(species = "human",
#'                    chain = "TRB",
#'                    region = "v",
#'                    sequence.type = "aa")
#' }
#'
#' @param species One or two-word common designation of species.
#' @param chain Sequence chain to access, e.g., \strong{TRB} or \strong{IGH}.
#' @param region Gene loci to access - \strong{v}, \strong{d}, \strong{j},
#' or \strong{c}.
#' @param sequence.type Type of sequence - \strong{aa} (amino acid) or
#' \strong{nt} (nucleotide).
#' @param refresh Logical. If \code{TRUE}, forces a re-download of the
#' data even if cached locally. Default is \code{FALSE}.
#'
#' @export
#' @return A list containing \code{sequences} (a named list of allele
#' sequences) and \code{misc} (metadata including species, chain,
#' sequence.type, and region).
getIMGT <- function(species = "human",
                    chain = "TRB",
                    sequence.type = "aa",
                    region = "v",
                    refresh = FALSE) {

  validate_input(region, c("v", "d", "j", "c"), "region")
  validate_input(sequence.type, c("aa", "nt"), "sequence.type")

  species.update <- .mapSpecies(species)
  gene <- toupper(paste0(chain, region))
  type <- ifelse(tolower(sequence.type) == "aa", "PROT", "NUC")

  result <- tryCatch(
    immReferent::getIMGT(
      species = species.update,
      gene = gene,
      type = type,
      refresh = refresh
    ),
    error = function(e) {
      warning("Failed to retrieve IMGT data: ", e$message)
      return(NULL)
    }
  )

  if (is.null(result)) {
    return(NULL)
  }

  sequences <- .biostringsToList(result, sequence.type)

  list(
    sequences = sequences,
    misc = list(
      species = species, chain = chain,
      sequence.type = sequence.type, region = region
    )
  )
}

# Helper function to validate input parameters
validate_input <- function(value, valid_options, parameter_name) {
  if (tolower(value) %!in% valid_options) {
    stop(sprintf("Invalid %s. Choose one of: %s",
                 parameter_name, paste(valid_options, collapse = ", ")))
  }
}

# Convert Biostrings object to named list of character sequences
.biostringsToList <- function(biostrings_obj, sequence.type) {
  seq_chars <- as.character(biostrings_obj)
  allele_names <- names(biostrings_obj)

  # Extract allele names from FASTA headers (pipe-delimited)
  parsed_names <- vapply(allele_names, function(nm) {
    parts <- strsplit(nm, "\\|")[[1]]
    if (length(parts) >= 2) parts[2] else nm
  }, character(1), USE.NAMES = FALSE)

  fasta_list <- as.list(seq_chars)
  names(fasta_list) <- parsed_names
  fasta_list
}

# Map immApex species names to immReferent species names
.mapSpecies <- function(x) {
  species_map <- list(
    "human" = "human",
    "mouse" = "mouse",
    "rat" = "rat",
    "rabbit" = "rabbit",
    "rhesus monkey" = "rhesus_monkey",
    "pig" = "pig",
    "dog" = "dog",
    "cyno monkey" = "cyno_monkey"
  )

  x <- tolower(x)
  if (x %in% names(species_map)) {
    return(species_map[[x]])
  } else {
    stop(sprintf("Invalid species. Choose one of: %s",
                 paste(names(species_map), collapse = ", ")))
  }
}
