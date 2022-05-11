
#' make_meta_dt
#'
#' @param dt data.table or matrix from load_*.
#' @param clin_dt clinical info, as from load_clinical_data().
#' @param include_all if TRUE all clinical attributes are included.
#'
#' @return data.table of clinical data fro all sample_ids
#' @export
#' @rdname make_meta_dt
#' @examples
#' mat = load_miR_RPM()
#' make_meta_dt(mat)
#'
#' dt = matrix_2_data.table(mat)
#' make_meta_dt(dt)
make_meta_dt = function(dt, clin_dt = load_clinical_data(), include_all = FALSE){
  if(is.matrix(dt)){
    mat = dt
    meta_dt = data.table(sample_id = colnames(mat))
  }else{
    meta_dt = unique(dt[, .(sample_id)])
  }
  meta_dt[, patient_id := sample_id_2_patient_id(sample_id)]
  if(include_all){
    meta_dt = merge(meta_dt, clin_dt, by = "patient_id")
  }else{
    meta_dt = merge(meta_dt, clin_dt[, .(patient_id, Phase, Vital.Status, Overall.Survival.Time.in.Days)], by = "patient_id")
  }
  setnames(meta_dt, c("Vital.Status", "Overall.Survival.Time.in.Days"), c("vital_status", "days_to_last_follow_up"))
  meta_dt[, days_to_death := ifelse(vital_status == "Dead", days_to_last_follow_up, NA) ]

  meta_dt[, sample_code := sample_id_2_sample_code(sample_id)]
  meta_dt = merge(meta_dt, load_sample_codes(), by = "sample_code")

  meta_dt[]
}

#' make_meta_dt.full
#' @export
#' @rdname make_meta_dt
#'
#' @examples
#'
#' make_meta_dt.full(mat)
#' make_meta_dt.full(dt)
#'
make_meta_dt.full = function(dt, clin_dt = load_clinical_data()){
  make_meta_dt(dt, clin_dt, include_all = TRUE)
}
