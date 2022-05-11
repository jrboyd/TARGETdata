#' run_DESeq
#'
#' @param mat
#' @param meta_dt
#' @param exp_cache_file
#' @param contrast_VAR
#'
#' @return
#' @export
#'
#' @examples
run_DESeq = function(mat, meta_dt, exp_cache_file, contrast_VAR = "ik_status"){
  if(file.exists(exp_cache_file)){
    exp_des = readRDS(exp_cache_file)
  }else{
    exp_des = DESeqDataSetFromMatrix(mat, meta_dt, design = paste0("~", contrast_VAR))
    exp_des = DESeq(exp_des)
    saveRDS(exp_des, exp_cache_file)
  }
  exp_des
}

#' write_DESeq_results
#'
#' @param des
#' @param analysis_name
#' @param contrast_VAR
#'
#' @return
#' @export
#'
#' @examples
write_DESeq_results = function(des, analysis_name, contrast_VAR = "ik_status"){
  conds = as.character(unique(colData(des)[[contrast_VAR]]))
  goi_todo.mirs = character()
  for(i in 1:(length(conds) - 1)){
    for(j in (i+1):length(conds)){
      a = conds[i]
      b = conds[j]
      print(paste(a, "vs", b))
      res <- DESeq2::results(des, contrast=c(contrast_VAR,a,b)) ## contrast specifies conditions to be tested

      write.table(res, res_file(paste0("miR_DESeq2_from_", a, "_to_", b, "_Full.txt")), sep = "\t", quote = FALSE)

      res.sig = subset(res, baseMean > 1 & padj < .05 & abs(log2FoldChange) > 1)
      write.table(res.sig, res_file(paste0("miR_DESeq2_from_", a, "_to_", b, "_Significant.txt")), sep = "\t", quote = FALSE)
      goi_todo.mirs = union(rownames(res.sig), goi_todo.mirs)
    }
  }
}
