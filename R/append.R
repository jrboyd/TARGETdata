
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
#' mat = filter_expression_to_valid(mat)
#' meta_dt = make_meta_dt(mat)
#'
#' dt = convert_matrix_2_data.table(mat)
#' make_meta_dt(dt)
make_meta_dt = function(dt,
                        clin_dt = load_clinical_data(),
                        splice_clust_dt = load_splice_clusters(),
                        include_all = FALSE,
                        extra_vars_included = character()){
  if(is.matrix(dt)){
    ids = colnames(dt)
  }else{
    ids = dt$sample_id
  }
  meta_dt = data.table(sample_id = ids)
  meta_dt[, patient_id := convert_sample_id_2_patient_id(sample_id)]
  if(include_all){
    meta_dt = merge(meta_dt, clin_dt, by = "patient_id")
  }else{
    meta_dt = merge(meta_dt,
                    clin_dt[, c("patient_id", "Phase", "Vital.Status", "Overall.Survival.Time.in.Days", extra_vars_included), with = FALSE],
                    by = "patient_id")
  }
  setnames(meta_dt, c("Vital.Status", "Overall.Survival.Time.in.Days"), c("vital_status", "days_to_last_follow_up"))
  meta_dt[, days_to_death := ifelse(vital_status == "Dead", days_to_last_follow_up, NA) ]

  meta_dt[, sample_code := convert_sample_id_2_sample_code(sample_id)]
  meta_dt = merge(meta_dt, load_sample_codes(), by = "sample_code")

  meta_dt = merge(meta_dt, splice_clust_dt, by = "sample_id")

  if(!setequal(ids, meta_dt$sample_id)){
    warning("Not all samples in expression data are present in clinical/meta data. Run filter_expression_to_valid.")
  }

  meta_dt[]
}

#' filter_expression_to_valid
#'
#' @param dt matrix or tidy data.table of expression data, as from load_RNA_counts() or convert_matrix_2_data.table()
#' @param clin_dt clinical data, as from load_clinical_data()
#'
#' @return Object of same class of dt with samples not covered by clin_dt removed
#' @export
#'
#' @examples
#' mat = load_RNA_counts()
#' clin_dt = load_clinical_data()
#' mat.valid = filter_expression_to_valid(mat, clin_dt)
filter_expression_to_valid = function(dt, clin_dt = load_clinical_data()){
  if(is.null(clin_dt$sample_id)){#if sample_id isn't present in clin_dt, use converted patient_id
    if(is.matrix(dt)){
      ids = colnames(dt)
    }else{
      ids = dt$sample_id
    }
    meta_dt = data.table(sample_id = ids)
    meta_dt[, patient_id := convert_sample_id_2_patient_id(sample_id)]

    k = meta_dt$patient_id %in% clin_dt$patient_id
    if(all(k)){
      return(dt)
    }else{
      message("Removing ", length(unique(meta_dt$patient_id[!k])), " of ", length(unique(meta_dt$patient_id)), " samples from expression due to missing from clinical.")
    }
    ids.valid = meta_dt$sample_id[k]

    if(is.matrix(dt)){
      return(dt[, ids.valid])
    }else{
      return(dt[ sample_id %in% ids.valid])
    }
  }else{#if sample_id is present in clin_dt, just use that
    if(is.matrix(dt)){
      n_original = ncol(dt)
      dt = dt[, clin_dt$sample_id]
      n_final = ncol(dt)
    }else{
      n_original = length(unique(dt$sample_id))
      dt[sample_id %in% clin_dt$sample_id]
      n_final = length(unique(dt$sample_id))
    }
    if(n_original != n_final){
      message("Removing ", n_original - n_final, " of ", n_original, " samples from expression due to missing from clinical.")
    }

    return(dt)
  }
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
