#' sample_id_2_patient_id
#'
#' @param ids extracts patient_id from sample IDs
#'
#' @return TARGET patient ids
#' @export
#'
#' @examples
#' sample_ids = c(
#'   "TARGET-15-SJMPAL043508-04B-01R",
#'   "TARGET-15-SJMPAL043511-09B-01R",
#'   "TARGET-15-SJMPAL043512-03A-01R")
#' sample_id_2_patient_id(sample_ids)
sample_id_2_patient_id = function(ids){
  ids = sub(".+TARGET", "TARGET", ids)
  sapply(strsplit(ids, "-"), function(x){
    paste(x[1:3], collapse = "-")
  })
}

#' sample_id_2_sample_code
#'
#' @param ids
#'
#' @return
#' @export
#'
#' @examples
#' sample_ids = c(
#'   "TARGET-15-SJMPAL043508-04B-01R",
#'   "TARGET-15-SJMPAL043511-09B-01R",
#'   "TARGET-15-SJMPAL043512-03A-01R")
#' sample_id_2_sample_code(sample_ids)
sample_id_2_sample_code = function(ids){
  ids = sub(".+TARGET", "TARGET", ids)
  codes = sapply(strsplit(ids, "-"), function(x){
    x[4]
  })
  gsub("\\..+", "", gsub("[A-Z]", "", codes))
}

#' data.table_to_matrix
#'
#' create a matrix from a data.table using gene_id as rownmaes
#'
#' @param dt
#'
#' @return
#' @export
#'
#' @rdname conversion
#' @examples
#' mat = load_miR_counts()
#' dt = matrix_2_data.table(mat)
#' mat2 = data.table_2_matrix(dt)
data.table_2_matrix = function(dt){
  if(!is.null(dt$sample_id)){
    dt.wide = dcast(dt, gene_id~sample_id)
  }else{
    dt.wide = dt
  }

  mat = as.matrix(dt.wide[, setdiff(colnames(dt.wide), "gene_id"), with = FALSE])
  rownames(mat) = dt.wide$gene_id
  mat
}

#' matrix_2_data.table
#'
#' @param mat
#'
#' @return
#' @export
#'
#' @rdname conversion
matrix_2_data.table = function(mat){
  dt = data.table::as.data.table(reshape2::melt(mat))
  setnames(dt, c("gene_id", "sample_id", "value"))
  dt[]
}
