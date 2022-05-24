library(data.table)
library(IKdata)

#' load_wgs_ids
#'
#' @return list containing patient ids, sample ids, and data.table of files
#' @export
#' @rdname load_wgs
#'
#' @examples
#' wgs_desc = load_wgs_ids()
#' fusion_dt = load_wgs_fusion()
#' fusion_patients = unique(fusion_dt$patient_id)
#'
#' seqsetvis::ssvFeatureVenn(list(
#'   fusion = fusion_patients,
#'   all = wgs_desc$patient_ids
#' ))
#'
#' vcf_dt = wgs_desc$vcf_dt
#' vcf_dt.dbl = vcf_dt[group == "double_sample"]
#' vcf_dt.dbl[, ids := sub("_TARG", " TARG", ids)]
#' sample_ids.dbl = unique(unlist(vcf_dt.dbl[, tstrsplit(ids, " ")]))
#'
#' seqsetvis::ssvFeatureVenn(list(
#'   double = sample_ids.dbl,
#'   all = wgs_desc$sample_ids
#' ))
#'
load_wgs_ids = function(){
  all_vcf = dir(dir(dir(dir("/slipstream/home/dbgap/globus", full.names = TRUE), full.names = TRUE), full.names = TRUE), full.names = TRUE)
  all(file.exists(all_vcf))

  all_vcf.samples = all_vcf[grepl("TARGET", basename(all_vcf))]

  vcf_dt = data.table(file = all_vcf.samples)

  vcf_dt[, ids := sub("\\..+", "", basename(file))]
  vcf_dt[, code := TARGETdata::convert_sample_id_2_sample_code(ids)]
  vcf_dt[is.na(code), group := "patient" ]
  vcf_dt[is.na(group), group := ifelse(grepl("TARGET.+TARGET", ids), "double_sample", "single_sample")]

  vcf_dt[, suffix := sub(ids, "", basename(file)), .(file)]

  vcf_dt.need_fix = vcf_dt[grepl("-", suffix)]
  vcf_dt.ok = vcf_dt[!grepl("-", suffix)]

  while(nrow(vcf_dt.need_fix) > 0){
    vcf_dt.need_fix[, fix := tstrsplit(suffix, "\\.", keep = 2)]
    vcf_dt.need_fix[, ids := paste0(ids, ".", fix)]
    vcf_dt.need_fix[, suffix := sub(paste0(".", fix), "", suffix), .(file)]
    vcf_dt.need_fix$fix = NULL

    vcf_dt.fixed = vcf_dt.need_fix[!grepl("-", suffix)]
    vcf_dt.need_fix = vcf_dt.need_fix[grepl("-", suffix)]

    vcf_dt.ok = rbind(vcf_dt.ok, vcf_dt.fixed)
  }

  vcf_dt.ok[, .N, .(suffix)]
  all_single_ids = unique(unlist(strsplit(vcf_dt.ok[group == "single_sample", .N, .(ids)]$ids, "_")))
  all_single_ids.patients = unique(TARGETdata::convert_sample_id_2_patient_id(all_single_ids))
  list(patient_ids = all_single_ids.patients, sample_ids = all_single_ids, vcf_dt = vcf_dt.ok)
}

#' load_wgs_fusion
#'
#' @return data.table of fusion calls from WGS TARGET
#' @export
#' @rdname load_wgs
#'
load_wgs_fusion = function(){
  vcf_dirs = dir("/slipstream/home/dbgap/globus", full.names = TRUE)
  vcf_dirs = file.path(vcf_dirs, c("structural/BCCA"))

  fusion_files = dir(vcf_dirs, pattern = "fusions.somatic.large.summary.tsv$", full.names = TRUE)
  patient_ids.fusion = sapply(strsplit(basename(fusion_files), "\\."), function(x)x[1])
  names(fusion_files) = patient_ids.fusion

  f = fusion_files[1]
  dtl = lapply(fusion_files, function(f){
    # message(f)
    dt = fread(f, sep = "\t")
    if(ncol(dt) > 1){
      dt = rbindlist(lapply(as.list(dt), function(x)data.table(x)))
      dt = dt[x != ""]
    }
    setnames(dt, "content")
    dt = dt[, tstrsplit(content, "[_|]")]
    setnames(dt, c("type", "posA", "posB", "geneA", "geneB", "dunno"))
    dt
  })
  fusion_dt = rbindlist(dtl, idcol = "patient_id")

  fusion_dt[]
}



