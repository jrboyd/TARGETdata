.run_DESeq = function(mat, meta_dt, contrast_VAR = "ik_status"){
  des = DESeq2::DESeqDataSetFromMatrix(mat, meta_dt, design = formula(paste0("~", contrast_VAR)))
  des = DESeq2::DESeq(des)
}

#' run_DESeq
#'
#' @param mat expression count matrix, use load_RNA_counts() or load_miR_counts()
#' @param meta_dt meta data, suitable for colData() when running DESeq, use load_clinical_data()
#' @param cache_file File to use for caching results.  Will be loaded if it exists. Default is NULL and disables caching.
#' @param contrast_VAR Variable in meta_dt to use for simple contrast. Default is "ik_status".
#'
#' @return result of DESeq
#' @export
#'
#' @examples
#' mat = load_RNA_counts()
#'
#' mat = filter_expression_to_valid(mat)
#' meta_dt = make_meta_dt(mat)
#'
#' #filter to no more than 10 per group
#' meta_dt = meta_dt[ Phase %in% c("Phase_II_Discovery", "Phase_II_Validation") & sample_type_short == "TBM"]
#' meta_dt[, rnk := frank(runif(.N), ties.method = "random"), .(ik_status)]
#' meta_dt = meta_dt[rnk <= 10]
#'
#' mat = filter_expression_to_valid(mat, meta_dt)
#'
#' des = run_DESeq(mat, meta_dt)
run_DESeq = function(mat, meta_dt, cache_file = NULL, contrast_VAR = "ik_status"){
  if(is.character(meta_dt[[contrast_VAR]])){
    meta_dt[[contrast_VAR]] = factor(meta_dt[[contrast_VAR]])
  }
  if(is.null(cache_file)){
    des = .run_DESeq(mat, meta_dt, contrast_VAR)
  }else{
    if(file.exists(cache_file)){
      des = readRDS(cache_file)
    }else{
      des = .run_DESeq(mat, meta_dt, contrast_VAR)
      saveRDS(des, cache_file)
    }
  }
  des
}

#' write_DESeq_results
#'
#' @param des A "DESeqDataSet" object, as returned from run_DESeq.
#' @param analysis_name Required prefix for output files.
#' @param contrast_VAR Variable in meta_dt to use for simple contrast. Default is "ik_status".
#' @param min_log2FoldChange Minimum absolute value of log2FoldChange. Default is 1 (2 fold).
#' @param max_padj Maximum padj. Default is .05.
#' @param min_baseMean Minimum baseMean. Default is 1.
#'
#' @return returns list of file path to full and significant results.
#' @export
#'
#' @examples
#' #' mat = load_RNA_counts()
#'
#' mat = filter_expression_to_valid(mat)
#' meta_dt = make_meta_dt(mat)
#'
#' #filter to no more than 10 per group
#' meta_dt = meta_dt[ Phase %in% c("Phase_II_Discovery", "Phase_II_Validation") & sample_type_short == "TBM"]
#' meta_dt[, rnk := frank(runif(.N), ties.method = "random"), .(ik_status)]
#' meta_dt = meta_dt[rnk <= 10]
#'
#' mat = filter_expression_to_valid(mat, meta_dt)
#'
#' des = run_DESeq(mat, meta_dt)
#' write_DESeq_results(des, "example)
write_DESeq_results = function(des, analysis_name, contrast_VAR = "ik_status", min_log2FoldChange = 1, max_padj = .05, min_baseMean = 1){
  conds = as.character(unique(SummarizedExperiment::colData(des)[[contrast_VAR]]))
  output_files.full = character()
  output_files.sig = character()
  for(i in 1:(length(conds) - 1)){
    for(j in (i+1):length(conds)){
      a = conds[i]
      b = conds[j]
      print(paste(a, "vs", b))
      res <- DESeq2::results(des, contrast=c(contrast_VAR,a,b)) ## contrast specifies conditions to be tested

      full_file = res_file(paste0(analysis_name, "_DESeq2_from_", a, "_to_", b, "_Full.txt"))
      output_files.full = c(output_files.full, full_file)
      write.table(res, full_file, sep = "\t", quote = FALSE)

      res.sig = subset(res, baseMean >= min_baseMean & padj <= max_padj & abs(log2FoldChange) >= min_log2FoldChange)

      sig_file = res_file(paste0(analysis_name, "_from_", a, "_to_", b, "_Significant.txt"))
      output_files.sig = c(output_files.sig, sig_file)
      write.table(res.sig, sig_file, sep = "\t", quote = FALSE)
    }
  }
  invisible(list(
    full = output_files.full,
    significant = output_files.sig
  ))
}
