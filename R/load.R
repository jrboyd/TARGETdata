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
#' load_clinical_data()
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
#' load_miR_counts()
load_miR_counts = function(){
  mir_cnt_dt = data.table::fread(file.path(data_path, "miRseq/miRseq_readCount.csv"))
  mir_cnt_dt[, miRNA_ID := paste(miRNA_ID, `cross-mapped`, sep = "-")]
  mir_cnt_dt$`cross-mapped` = NULL
  data.table::setnames(mir_cnt_dt, "miRNA_ID", "gene_id")
  convert_data.table_2_matrix(mir_cnt_dt[])
}

#' load_miR_RPM
#'
#' @return
#' @export
#'
#' @examples
#' load_miR_RPM()
load_miR_RPM = function(){
  mir_rpm_dt = data.table::fread(file.path(data_path, "miRseq/miRseq_RPM.csv"))
  mir_rpm_dt[, miRNA_ID := paste(miRNA_ID, `cross-mapped`, sep = "-")]
  mir_rpm_dt$`cross-mapped` = NULL
  data.table::setnames(mir_rpm_dt, "miRNA_ID", "gene_id")
  convert_data.table_2_matrix(mir_rpm_dt[])
}

#' load_RNA_counts
#'
#' @return
#' @export
#'
#' @examples
#' load_RNA_counts()
load_RNA_counts = function(){
  exp_cnt_dt = data.table::fread(file.path(data_path, "TARGET_expression.csv"))
  colnames(exp_cnt_dt)[-1] = sapply(strsplit(colnames(exp_cnt_dt)[-1], "\\."), function(x)x[3])
  convert_data.table_2_matrix(exp_cnt_dt[])
}

#' load_RNA_RPM
#'
#' @return
#' @export
#'
#' @examples
#' load_RNA_RPM()
load_RNA_RPM = function(){
  exp_cnt_dt = load_RNA_counts()
  exp_rpm_dt = data.table::as.data.table(reshape2::melt(exp_cnt_dt, id.vars = "gene_id"))
  data.table::setnames(exp_rpm_dt, c("gene_id", "sample_id", "value"))
  exp_rpm_dt[, rpm := value / sum(value) * 1e6, list(sample_id)]
  exp_rpm_dt = dcast(exp_rpm_dt, gene_id~sample_id, value.var = "rpm")
  convert_data.table_2_matrix(exp_rpm_dt[])
}

#' load_ref_gr
#'
#' @return
#' @export
#'
#' @examples
#' load_ref_gr
load_ref_gr = function(){
  ref_gr = rtracklayer::import.gff(file.path(data_path, "gencode.v36.annotation.gtf"), feature.type = "gene")
  names(ref_gr) = ref_gr$gene_id
  ref_gr
}

#' load_splice_clusters
#'
#' @return
#' @export
#'
#' @examples
#' load_splice_clusters()
load_splice_clusters = function(){
  splice_clust_dt = data.table::fread(file.path(data_path, "TARGET_IKZF1_splicing_clusters.csv"))
  splice_clust_dt[, sample_id := tstrsplit(sample, "\\.", keep = 3)]
  splice_clust_dt[, ik_status := "normal_splicing"]
  splice_clust_dt[cluster_id == 3, ik_status := "IK6"]
  splice_clust_dt[cluster_id %in% c(1, 2), ik_status := "early_termination"]
  splice_clust_dt[, list(sample_id, cluster_id, ik_status)]
}


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





