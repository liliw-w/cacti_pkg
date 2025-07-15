## File: R/helpers.R

#' Generate sliding-window peak groups
#'
#' Divides genomic peaks into non-overlapping windows to define groups of adjacent regulatory elements.
#'
#' @param file_pheno_meta Path to phenotype metadata with columns: phe_id, phe_chr, phe_start, phe_end.
#' @param window_size Window size specified as a string in kb (e.g., "50kb").
#'
#' @return A tibble with original columns plus `phe_group` indicating group.
#'
#' @importFrom data.table fread
#' @importFrom dplyr arrange group_by mutate ungroup
#' @importFrom stringr str_detect str_extract
#' @export
make_peak_group <- function(file_pheno_meta, window_size = "50kb") {
  # Parse window size
  if (str_detect(window_size, '^\\d+(\\.\\d*)?[kK][bB]$')) {
    ws <- as.numeric(str_extract(window_size, '^\\d+(\\.\\d*)?')) * 1000
  } else {
    stop("window_size must be in kb, e.g. '50kb'")
  }
  
  # Load metadata
  pheno_meta <- data.table::fread(file_pheno_meta)
  
  if (any(pheno_meta$phe_start >= pheno_meta$phe_end)) {
    stop("All phe_start must be < phe_end")
  }
  
  # For each chromosome, create breaks and assign groups
  pheno_meta <- pheno_meta %>%
    arrange(phe_chr, phe_start) %>%
    group_by(phe_chr) %>%
    mutate(
      phe_group = cut(
        phe_start,
        breaks = c(seq(0, max(phe_start) + 1, by = ws), Inf),
        labels = paste0(unique(phe_chr), seq_along(seq(0, max(phe_start) + 1, by = ws))),
        right = FALSE
      )
    ) %>%
    ungroup()
  return(pheno_meta)
}

#' Compute phenotype residuals after regressing out covariates
#'
#' Reads phenotype and covariate matrices, regresses each phenotype on covariates, and returns residuals.
#'
#' @param file_pheno Phenotype input: Tab-delimited file, samples in rows, phenotypes in columns.
#' @param file_cov Covariate input: samples in rows, covariates in columns.
#'
#' @return A matrix of residuals (samples x phenotypes).
#'
#' @importFrom utils read.table
#' @export
make_pheno_cov_residual <- function(file_pheno, file_cov) {
  pheno <- read.table(file_pheno, header = TRUE, sep = "\t", quote = "", row.names = 1)
  cov   <- read.table(file_cov,   header = TRUE, sep = "\t", quote = "", row.names = 1)
  
  # Ensure samples match
  pheno_mat <- t(as.matrix(pheno))
  cov_mat   <- t(as.matrix(cov))[rownames(pheno_mat), , drop = FALSE]
  
  # Regress out covariates for each phenotype
  residuals <- apply(pheno_mat, 2, function(y) lm(y ~ cov_mat)$residuals)
  
  return(residuals)
}

#' Add FDR per window
#'
#' Reads per-group SNP p-value files, computes ACAT combined p-value and q-values.
#'
#' @param file_p_peak_group Character vector of file paths. Each element is a path to a p-value file with columns: group, snp, p.
#'
#' @return A tibble with columns: group, n_var_in_cis, minp_nom, minp_acat, snp, q.
#'
#' @importFrom data.table fread
#' @importFrom dplyr group_by summarise ungroup mutate bind_rows
#' @importFrom purrr map_dfr
#' @importFrom ACAT ACAT
#' @importFrom qvalue qvalue
#' @export
add_fdr_cwindow <- function(file_p_peak_group) {
  df <- purrr::map_dfr(file_p_peak_group, function(x) {
    tmp_p <- data.table::fread(x) %>%
      group_by(group) %>%
      mutate(n_var_in_cis = n()) %>%
      ungroup() %>%
      mutate(p = ifelse(p == 1, 1 - 1/n_var_in_cis, p))
    
    tmp_p %>%
      group_by(group) %>%
      summarise(
        n_var_in_cis = first(n_var_in_cis),
        minp_nom = min(p),
        minp_acat = ACAT(Pvals = p),
        snp = snp[which.min(p)]
      ) %>%
      ungroup()
  })
  
  df$q <- qvalue::qvalue(df$minp_acat)$qvalues
  
  return(df)
}

