##############################################
########### Run PC-based test ###########
##############################################
# load packages -----
library(tidyverse)
library(qvalue)
library(ACAT)

dir_pco <- 'script/pco/'
## load pco
source("script/pco/ModifiedPCOMerged_acat.R")
source("script/pco/liu.R")
source("script/pco/liumod.R")
source("script/pco/davies.R")
dyn.load("script/pco/qfc.so")


cal_p <- function(
  chr, 
  file_qtl_cis_norm, 
  # cols - phe_id, phe_chr, var_id, nom_pval, slope
  
  file_peak_group = NULL, 
  # cols - phe_id, phe_group
  
  file_pheno_cov_residual = NULL,
  # sample by phenos, cols - phe_id
  
  n_pid = 2
  # number of peaks in the group considered for multivariate test
) {
  # read files -----
  qtl_cis_norm = data.table::fread(file_qtl_cis_norm, col.names = TRUE)
  
  if(is.null(file_peak_group)) 
    peak_group = make_peak_group()
  else
    peak_group = data.table::fread(file_peak_group, col.names = TRUE)
  
  if(is.null(file_pheno_cov_residual)) 
    pheno_cov_residual = make_pheno_cov_residual()
  else
    pheno_cov_residual = data.table::fread(file_pheno_cov_residual, col.names = TRUE)
  
  # organize data -----
  # add column of z
  qtl_cis_norm = filter(qtl_cis_norm, phe_chr == !!chr) %>%
    mutate(z = sqrt(qchisq(nom_pval, df = 1, lower.tail = FALSE)) * sign(slope))
  
  # select peak groups on the given chr
  peak_group = filter(peak_group, phe_id %in% unique(qtl_cis_norm$phe_id)) %>%
    group_by(phe_group) %>%
    mutate(size_phe_group = n()) %>%
    ungroup() %>%
    mutate(phe_group = str_glue("G_{phe_group}_n{size_phe_group}"))
  
  group_uniq = unique(peak_group$phe_group)
  n_group = length(group_uniq)
  
  
  # calculate p-values across peak groups -----
  p_peak_group_list = list()
  k = 0
  for(g in group_uniq){
    ## progress message
    k = k + 1
    if(k %% 100 == 0) cat(k, "-th group is running, out of", n_group, "group, on", chr, ". \n")
    
    
    # 1. Input Z matrix, peaks in columns, SNPs in rows
    tmp_phe_id = filter(peak_group, phe_group == !!g) %>% pull(phe_id)
    
    # single peak test
    if(length(tmp_phe_id) < n_pid) {
      p_peak_group_list[[g]] = filter(qtl_cis_norm, phe_id %in% !!tmp_phe_id) %>%
        select(var_id, nom_pval) %>%
        mutate("group" = !!g) %>%
        rename("snp" = 'var_id', 'p' = 'nom_pval')
      
      next
    }
    
    # multi peaks test
    z_mat = filter(qtl_cis_norm, phe_id %in% !!tmp_phe_id) %>%
      select(phe_id, var_id, z) %>%
      pivot_wider(names_from = phe_id, values_from = z) %>%
      filter(complete.cases(.)) %>%
      column_to_rownames(var = "var_id") %>%
      as.matrix()
    
    ## check if there is any overlapped snps for the peak group
    if(nrow(z_mat) == 0 | ncol(z_mat) < n_pid) next
    
    
    # 2. estimated Sigma of the group
    Sigma = cor(pheno_cov_residual[, colnames(z_mat)])
    
    
    # 3. Run test for a peak group
    p_peak_group_list[[g]] = enframe(
      ModifiedPCOMerged_acat(Z.mat = z_mat, Sigma = Sigma), 
      name = "snp", value = "p"
    ) %>%
      mutate("group" = !!group)
  }
  
  
  # return results -----
  return(
    list(
      'p_peak_group' = bind_rows(p_peak_group_list) %>%
        relocate(group, .before = everything()),
      'peak_group' = peak_group,
      'pheno_cov_residual' = pheno_cov_residual
    )
  )
}


make_peak_group <- function(
  file_pheno_meta, 
  # cols - phe_id, phe_chr, phe_start, phe_end
  
  window_size = '50kb'
) {
  # re-format window size
  if(str_detect(window_size, '^\\d+[.]?\\d*[k, K][b, B]$')){
    window_size = as.numeric(str_extract(window_size, '^\\d+[.]?\\d*')) * 1000
  }else stop("Specify window size in kb. \n")
  
  # read files -----
  pheno_meta = data.table::fread(file_pheno_meta, col.names = TRUE)
  
  # group adjacent peaks into a group by a sliding window -----
  if(any(pheno_meta$phe_start >= pheno_meta$phe_end)) stop("This pipeline requires feature start postion smaller than end postion.")
  
  
  # cluster peaks by dividing the genome using non-overlapping windows
  peak_group = arrange(pheno_meta, phe_chr, phe_start) %>%
    group_by(phe_chr) %>%
    mutate(
      "phe_group" = cut(
        phe_start, 
        breaks = c(seq(0, max(phe_start)+1, by = window_size), Inf), 
        labels = paste0(unique(phe_chr), seq_along(seq(0, max(phe_start)+1, by = window_size))),
        right = FALSE
      )
    ) %>%
    ungroup()
  
  return(peak_group)
}

make_pheno_cov_residual <- function(
  file_pheno,
  # pheno by sample, rownames are pheno, colnames are sample names
  
  file_cov
  # cov by sample, rownames are cov, colnames are sample names
) {
  # read files -----
  pheno = read.table(
    file_pheno, 
    header = TRUE, sep = "\t", quote = "", row.names = 1
  ) %>%
    rownames_to_column(var = "ID")
  
  cov = read.table(
    file_cov, 
    header = TRUE, sep = "\t", quote = "", row.names = 1
  ) %>%
    rownames_to_column(var = "ID")
  
  
  # make tibble a matrix for regression -----
  pheno_mat = select(pheno, -(ID)) %>%
    as.matrix() %>%
    t()
  colnames(pheno_mat) = pheno$ID
  
  cov_mat = select(cov, -(ID)) %>%
    as.matrix() %>%
    t()
  colnames(cov_mat) = cov$ID
  
  # match sample order
  cov_mat = cov_mat[match(rownames(pheno_mat), rownames(cov_mat)), , drop = FALSE]
  
  # regress out covariates -----
  pheno_cov_residual = apply(pheno_mat, 2, function(y) extract_residual(y, cov_mat))
  
  
  return(pheno_cov_residual)
}


extract_residual <- function(y, x){
  return(lm(y ~ x)$residuals)
}


add_fdr_cwindow <- function(
  file_p_peak_group
  # vector of files, each file cols - group, snp, p
) {
  p_fdr_added = lapply(
    X = file_p_peak_group,
    FUN = function(x){
      tmp_p = data.table::fread(x) %>%
        group_by(group) %>%
        mutate(n_var_in_cis = n()) %>%
        ungroup()
      
      # top p correction of each group -----
      ## change 1 p value for acat
      tmp_p$p[tmp_p$p == 1] = 1 - 1/tmp_p$n_var_in_cis[tmp_p$p == 1]
      
      tmp_p = group_by(tmp_p, group) %>%
        summarise(
          n_var_in_cis = n(),
          minp_nom = min(p),
          minp_acat = ACAT::ACAT(Pvals = p),
          snp = snp[which.min(p)]
        ) %>%
        ungroup()
    }
  ) %>%
    bind_rows()
  
  # q-values -----
  p_fdr_added$q = qvalue::qvalue(p_fdr_added$minp_acat)$qvalues
  
  return(p_fdr_added)
}

