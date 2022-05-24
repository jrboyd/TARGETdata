#### Setup ####
library(TARGETdata)
library(ggalluvial)

#### Array Data ####
library(data.table)
wxs_clin_dt = fread("../TARGETdata_files/WXS_clinical.txt")
ids_wxs = unique(wxs_clin_dt$case_submitter_id)

arr_clin_dt = fread("../TARGETdata_files/GenotypeArray_clinical.txt")
ids_arr_genotype =  unique(arr_clin_dt$case_submitter_id)

arr_xlsx = openxlsx::read.xlsx("../TARGETdata_files/TARGET_ALL_CopyNumberArray_Phase2_20160812.xlsx", sheet = 2)
arr_xlsx[1:10, 1:10]
ids_arr_copynumber = unique(arr_xlsx$Source.Name)


#Objective, describe breakdown of TARGET data regarding sequencing type (RNAseq,
#miRseq, WGS).


#### Seq Data ####

if(!exists("rpm_mat.full")){
  message("loading RNA RPM")
  rpm_mat.full = load_RNA_RPM()
}
if(!exists("clin_dt")){
  message("loading clin_dt")
  clin_dt = load_clinical_data()
}
if(!exists("rpm_mat")){
  rpm_mat = filter_expression_to_valid(rpm_mat.full, clin_dt)
}
if(!exists("meta_dt")){
  meta_dt = make_meta_dt(rpm_mat, clin_dt, extra_vars_included = c("Cell.of.Origin"))
  setnames(meta_dt, "Cell.of.Origin", "cell_of_origin")
}

if(!exists("rpm_mat.mir")){
  message("loading miR RPM")
  rpm_mat.mir = load_miR_RPM()
}

active_filters = c("Phase", "sample_type", "cell_of_origin", "ik_status")

type_dt = unique(meta_dt[, .(sample_code, sample_type_short, sample_type)])
type_dt[, caption := paste0(sample_type_short, paste(rep(" ", 5 - nchar(sample_type_short)), collapse = ""), "= ", sample_type), .(sample_type_short)]
nchar(type_dt$caption)
type_dt[, caption := paste0(caption, paste0(rep(" ", 57 - nchar(caption)), collapse = "")), .(caption)]
type_dt

if(!exists("wgs_data")){
  message("loading WGS")
  wgs_data = load_wgs_ids()
}

convert_sample_id_2_code_id = function(ids){
  ids = sub(".+TARGET", "TARGET", ids)
  out = sapply(strsplit(ids, "-"), function(x){
    paste(x[1:4], collapse = "-")
  })
  out = sub("\\..+", "", out)
  out = sub("[A-Z]$", "", out)
  out

}

if(TRUE){#sample codes between RNA/MIR and WGS are not compatible beyond sample code
  ids_wgs = convert_sample_id_2_code_id(wgs_data$sample_ids)
  ids_mir = convert_sample_id_2_code_id(colnames(rpm_mat.mir))
  ids_rna = convert_sample_id_2_code_id(colnames(rpm_mat.full))
}else{
  ids_wgs = wgs_data$sample_ids
  ids_mir = colnames(rpm_mat.mir)
  ids_rna = colnames(rpm_mat.full)

  intersect(wgs_data$sample_ids, colnames(rpm_mat.full))
  intersect(wgs_data$sample_ids, colnames(rpm_mat.mir))
}


#### Compare ####
cov_dt = rbind(
  data.table(seq = "WXS", patient_id = ids_wxs),
  data.table(seq = "ArrayGenotype", patient_id = ids_arr_genotype),
  data.table(seq = "ArrayCopyNumber", patient_id = ids_arr_copynumber),
  data.table(seq = "WGS", patient_id = convert_sample_id_2_patient_id(ids_wgs)),
  data.table(seq = "miR", patient_id = convert_sample_id_2_patient_id(ids_mir)),
  data.table(seq = "RNA", patient_id = convert_sample_id_2_patient_id(ids_rna))
)
cov_dt = unique(cov_dt)

table(table(cov_dt$patient_id))

cov_dt$present = "present"
cov_dt = dcast(cov_dt, patient_id~seq, value.var = "present", fill = "-")

full_clin_dt = clin_dt[, .(patient_id, Cell.of.Origin, Phase)]
full_clin_dt = merge(full_clin_dt, cov_dt, by = "patient_id", all.x = TRUE)

full_clin_dt[is.na(RNA), RNA := "-"]
full_clin_dt[is.na(WXS), WXS := "-"]
full_clin_dt[is.na(ArrayGenotype), ArrayGenotype := "-"]
full_clin_dt[is.na(ArrayCopyNumber), ArrayCopyNumber := "-"]
full_clin_dt[is.na(WGS), WGS := "-"]
full_clin_dt[is.na(miR), miR := "-"]

full_clin_dt.alluv = full_clin_dt[, .N, .(Cell.of.Origin, Phase, ArrayCopyNumber, ArrayGenotype, RNA, WGS, WXS, miR)]

theme_set(theme(
  panel.background = element_blank(),
  panel.grid = element_blank(),
  legend.position = "bottom",
  legend.background = element_blank(),
  legend.box.background = element_blank(),
  legend.key = element_blank(),
  plot.title = element_text(size = 10),
  plot.subtitle = element_text(size = 8),
  legend.text = element_text(size = 8),
  axis.text = element_text(size = 8, angle = 45, hjust = 1, vjust = 1),
  plot.caption = element_text(family = "mono", size = 8)
))

ax1 = "Phase"
ax2 = "Cell.of.Origin"
ax3 = "ArrayCopyNumber"
ax4 = "WXS"
ax5 = "RNA"
ax6 = "ArrayGenotype"
ax7 = "miR"
ax8 = "WGS"

full_clin_dt.alluv

sel = "Phase"
ggplot(full_clin_dt.alluv, aes_string(axis1 = ax1, axis2 = ax2, axis3 = ax3, axis4 = ax4, axis5 = ax5, axis6 = ax6, axis7 = ax7, axis8 = ax8, y = "N")) +
  ggalluvial::geom_alluvium(aes_string(fill = sel), width = 1/5) +
  ggalluvial::geom_stratum(width = 1/8, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum)), size = 1.8) +
  scale_x_discrete(limits = c(ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8), expand = c(.15, .15)) +
  # scale_fill_manual(values = c("-" = "gray", "present" = "red")) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  guides(fill = guide_legend(nrow = 3)) +
  labs(title = "Phase vs Cell of Origin")

sel = ax3
plots = lapply(c(ax3, ax4, ax5, ax6, ax7, ax8), function(sel){
  ggplot(full_clin_dt.alluv, aes_string(axis1 = ax1, axis2 = ax2, axis3 = ax3, axis4 = ax4, axis5 = ax5, axis6 = ax6, axis7 = ax7, axis8 = ax8, y = "N")) +
    ggalluvial::geom_alluvium(aes_string(fill = sel), width = 1/5) +
    ggalluvial::geom_stratum(width = 1/8, fill = "black", color = "grey") +
    geom_label(stat = "stratum", aes(label = after_stat(stratum)), size = 1.8) +
    scale_x_discrete(limits = c(ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8), expand = c(.15, .15)) +
    scale_fill_manual(values = c("-" = "gray", "present" = "red")) +
    guides(fill = "none") +
    # scale_fill_brewer(type = "qual", palette = "Set1") +
    labs(title = sel)
})

plots = lapply(c(ax3, ax4, ax5, ax6, ax7, ax8), function(sel){
  ggplot(full_clin_dt.alluv[get(sel) == "present"], aes_string(axis1 = ax1, axis2 = ax2, axis3 = ax3, axis4 = ax4, axis5 = ax5, axis6 = ax6, axis7 = ax7, axis8 = ax8, y = "N")) +
    ggalluvial::geom_alluvium(aes_string(fill = sel), width = 1/5) +
    ggalluvial::geom_stratum(width = 1/8, fill = "black", color = "grey") +
    geom_label(stat = "stratum", aes(label = after_stat(stratum)), size = 1.8) +
    scale_x_discrete(limits = c(ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8), expand = c(.15, .15)) +
    scale_fill_manual(values = c("-" = "gray", "present" = "red")) +
    guides(fill = "none") +
    # scale_fill_brewer(type = "qual", palette = "Set1") +
    labs(title = sel)
})

phase_cols = seqsetvis::safeBrew(full_clin_dt.alluv$Phase)

plots = lapply(c(ax3, ax4, ax5, ax6, ax7, ax8), function(sel){
  ggplot(full_clin_dt.alluv[get(sel) == "present"], aes_string(axis1 = ax1, axis2 = ax2, axis3 = ax3, axis4 = ax4, axis5 = ax5, axis6 = ax6, axis7 = ax7, axis8 = ax8, y = "N")) +
    ggalluvial::geom_alluvium(aes(fill = Phase), width = 1/5) +
    ggalluvial::geom_stratum(width = 1/8, fill = "black", color = "white", size = 1.5) +
    geom_label(stat = "stratum", aes(label = after_stat(stratum)), size = 1.8) +
    scale_x_discrete(limits = c(ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8), expand = c(.15, .15)) +
    scale_fill_manual(values = phase_cols) +
    guides(fill = "none") +
    # scale_fill_brewer(type = "qual", palette = "Set1") +
    labs(title = sel)
})
gc()
pg_plots = cowplot::plot_grid(plotlist = plots)
ggsave("array_data_alluvial.png", pg_plots, width = 14, height = 7)



full_clin_dt.alluv$id = seq_len(nrow(full_clin_dt.alluv))
full_clin_dt.alluv[, RNA_and_WGS := {
  if(RNA == "present" & WGS == "present"){
    "RNA and WGS"
  }else if(RNA == "present"){
    "RNA"
  }else if(WGS == "present"){
    "WGS"
  }else{
    "none"
  }
}, .(id)]
gc()
p_rna_and_wgs_alluvial = ggplot(full_clin_dt.alluv,#[get(sel) == "present"],
                                aes_string(axis1 = ax5, axis2 = ax8, axis3 = ax3, axis4 = ax4, axis5 = ax6, axis6 = ax7, y = "N")) +
  ggalluvial::geom_alluvium(aes_string(fill = "RNA_and_WGS"), width = 1/5) +
  ggalluvial::geom_stratum(width = 1/8, fill = "black", color = "white", size = 1.5) +
  geom_label(stat = "stratum", aes(label = after_stat(stratum)), size = 1.8) +
  scale_x_discrete(limits = c(ax5, ax8, ax3, ax4, ax6, ax7), expand = c(.15, .15)) +
  scale_fill_manual(values = c("RNA and WGS" = "purple", "RNA" = "blue", "WGS" = "red", "none" = "gray"))
p_rna_and_wgs_alluvial
ggsave("RNA_and_WGS.png", p_rna_and_wgs_alluvial, width = 4.7, height = 4.4)

full_clin_dt.alluv[RNA_and_WGS == "RNA and WGS"]
full_clin_dt.alluv[grepl("Phase_II_", Phase)]


axes = aes_string(axis1 = ax5, axis2 = ax3, axis3 = ax4, axis4 = ax6, axis6 = ax7, y = "N")
lims = c(ax5, ax3, ax4, ax6, ax7)

plots = lapply(c(ax3, ax4, ax6), function(ax){
  ggplot(full_clin_dt.alluv[WGS == "-" & miR == "present"],#[get(sel) == "present"],
         axes) +
    ggalluvial::geom_alluvium(aes_string(fill = ax), width = 1/5) +
    ggalluvial::geom_stratum(width = 1/8, fill = "black", color = "white", size = 1.5) +
    geom_label(stat = "stratum", aes(label = after_stat(stratum)), size = 1.8) +
    scale_x_discrete(limits = lims, expand = c(.15, .15)) +
    scale_fill_manual(values = c("-"= "gray", "present" = "green"))
})
pg_new_data =
  cowplot::plot_grid(ggplot() + labs(title = "miR samples without WGS") + theme_void(),
                     ncol = 1, rel_heights = c(1, 10),
                     cowplot::plot_grid(plotlist = plots, nrow = 1)
  )
ggsave("unexplored_data.miR_without_WGS.png", pg_new_data, width = 11, height = 4)

ggplot(full_clin_dt.alluv,#[get(sel) == "present"],
       axes) +
  ggalluvial::geom_alluvium(aes_string(fill = ax3), width = 1/5) +
  ggalluvial::geom_stratum(width = 1/8, fill = "black", color = "white", size = 1.5) +
  geom_label(stat = "stratum", aes(label = after_stat(stratum)), size = 1.8) +
  scale_x_discrete(limits = lims, expand = c(.15, .15)) +
  scale_fill_manual(values = c("-"= "gray", "present" = "green"))

ggplot(full_clin_dt.alluv,#[get(sel) == "present"],
       axes) +
  ggalluvial::geom_alluvium(aes_string(fill = ax4), width = 1/5) +
  ggalluvial::geom_stratum(width = 1/8, fill = "black", color = "white", size = 1.5) +
  geom_label(stat = "stratum", aes(label = after_stat(stratum)), size = 1.8) +
  scale_x_discrete(limits = lims, expand = c(.15, .15)) +
  scale_fill_manual(values = c("-"= "gray", "present" = "green"))

ggplot(full_clin_dt.alluv,#[get(sel) == "present"],
       axes) +
  ggalluvial::geom_alluvium(aes_string(fill = ax6), width = 1/5) +
  ggalluvial::geom_stratum(width = 1/8, fill = "black", color = "white", size = 1.5) +
  geom_label(stat = "stratum", aes(label = after_stat(stratum)), size = 1.8) +
  scale_x_discrete(limits = lims, expand = c(.15, .15)) +
  scale_fill_manual(values = c("-"= "gray", "present" = "green"))

# scale_fill_manual(values = phase_cols) +
# guides(fill = "none") #+
# scale_fill_brewer(type = "qual", palette = "Set1") +
# labs(title = sel)
