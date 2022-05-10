data_path = "/slipstream_old/home/joeboyd/R_workspace/TARGETdata_files"

#' load_sample_codes
#'
#' @return
#' @export
#' @importFrom data.table fread
#'
#' @examples
#' load_sample_codes()
load_sample_codes = function(){
  samp_dt = data.table::fread(file.path(data_path, "sample_codes.csv"),
                  colClasses = c("character"),
                  col.names = c("sample_code", "sample_type", "sample_type_short"))
  samp_dt[]
}

#' load_clinical_data
#'
#' @return
#' @export
#'
#' @examples
load_clinical_data = function(){
  clin_dt = data.table::fread(file.path(data_path, "clinical_merged.all.csv"))
  clin_dt[, Phase := sub("_2.+", "", sub(".+ClinicalData_", "", file))]
  clin_dt[]
}

#' load_miR_counts
#'
#' @return
#' @export
#'
#' @examples
load_miR_counts = function(){
  mir_cnt_dt = data.table::fread(file.path(data_path, "miRseq/miRseq_readCount.csv"))
  mir_cnt_dt[]
}

#' load_miR_RPM
#'
#' @return
#' @export
#'
#' @examples
load_miR_RPM = function(){
  mir_rpm_dt = data.table::fread(file.path(data_path, "miRseq/miRseq_RPM.csv"))
  mir_rpm_dt[]
}

#' load_RNA_counts
#'
#' @return
#' @export
#'
#' @examples
load_RNA_counts = function(){
  exp_cnt_dt = data.table::fread(file.path(data_path, "TARGET_expression.csv"))
  colnames(exp_cnt_dt)[-1] = sapply(strsplit(colnames(exp_cnt_dt)[-1], "\\."), function(x)x[3])
  exp_cnt_dt[]
}

#' load_RNA_RPM
#'
#' @return
#' @export
#'
#' @examples
load_RNA_RPM = function(){
  exp_cnt_dt = load_RNA_counts()
  exp_rpm_dt = melt(exp_cnt_dt, id.vars = "gene_id")
  exp_rpm_dt[, rpm := value / sum(value) * 1e6, .(variable)]
  exp_rpm_dt = dcast(exp_rpm_dt, gene_id~variable, value.var = "rpm")
}




