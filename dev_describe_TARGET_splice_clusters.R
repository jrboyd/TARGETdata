
#Objective, describe breakdown of TARGET data regarding sequencing type (RNAseq,
#miRseq, WGS).

#### Setup ####
library(TARGETdata)
library(ggalluvial)

#### Load Data ####

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

sample_codes_dt = load_sample_codes()

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

#### WGS data has distinct codes ####

ids_wgs.normals = ids_wgs[grepl("(14$)|(15$)|(10$)", ids_wgs)]
ids_wgs = ids_wgs[!grepl("(14$)|(15$)|(10$)", ids_wgs)]
seqsetvis::ssvFeatureVenn(list(
  WGS = ids_wgs, miR = ids_mir, RNA = ids_rna
))


intersect(ids_rna, ids_wgs)
setdiff(ids_wgs, ids_rna)

#### Test Plots ####

as.data.frame(UCBAdmissions)
ggplot(as.data.frame(UCBAdmissions),
       aes(y = Freq, axis1 = Gender, axis2 = Dept)) +
  ggalluvial::geom_alluvium(aes(fill = Admit), width = 1/12) +
  ggalluvial::geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Gender", "Dept"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  ggtitle("UC Berkeley admissions and rejections, by sex and department")


ggplot(meta_dt, aes(axis1 = Phase, axis2 = ik_status)) +
  ggalluvial::geom_alluvium(aes(fill = Phase), width = 1/5) +
  ggalluvial::geom_stratum(width = 1/8, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Phase", "ik_status"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1")


ggplot(meta_dt[grepl("Phase_II_", Phase)], aes(axis1 = Phase, axis2 = ik_status)) +
  ggalluvial::geom_alluvium(aes(fill = Phase), width = 1/5) +
  ggalluvial::geom_stratum(width = 1/8, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Phase", "ik_status"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1")

ax1 = "Phase"
ax2 = "sample_type_short"
ax3 = "ik_status"

ggplot(meta_dt[grepl("Phase_II_", Phase)], aes_string(axis1 = ax1, axis2 = ax2, axis3 = ax3)) +
  ggalluvial::geom_alluvium(aes(fill = Phase), width = 1/5) +
  ggalluvial::geom_stratum(width = 1/8, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c(ax1, ax2, ax3), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  labs(caption = paste(type_dt$caption, collapse = "\n")) +
  theme(plot.caption = element_text(family = "mono"))



#### sample_id coverage ####

cov_dt = rbind(
  data.table(seq = "WGS", sample_id = ids_wgs),
  data.table(seq = "miR", sample_id = ids_mir),
  data.table(seq = "RNA", sample_id = ids_rna)
)
cov_dt = unique(cov_dt)
cov_dt$present = "present"
cov_dt = dcast(cov_dt, sample_id~seq, value.var = "present", fill = "-")

ax1 = "RNA"
ax2 = "miR"
ax3 = "WGS"

universal_ids = cov_dt[RNA != "-" & WGS != "-" & miR != "-"]$sample_id
cov_dt[, universal := sample_id %in% universal_ids]
cov_dt[, presence := sum(RNA == "present") + sum(WGS == "present") + sum(miR == "present"), .(sample_id)]
cov_dt$presence = factor(cov_dt$presence)
cov_dt[, RNA_and_miR := RNA == "present" & miR == "present"]
cov_dt[, RNA_and_WGS := RNA == "present" & WGS == "present"]
cov_dt[, RNA_and_miR_and_WGS := RNA == "present" & WGS == "present"  & miR == "present"]

splice_dt = load_splice_clusters()
splice_dt[, patient_id := convert_sample_id_2_patient_id(sample_id)]


meta_dt[, sample_id := convert_sample_id_2_code_id(sample_id)]
meta_dt = merge(meta_dt, cov_dt[, .(sample_id, RNA, WGS, miR)], by = 'sample_id')

todo = unique(meta_dt$ik_status)
todo

phase_colors = seqsetvis::safeBrew(meta_dt$Phase)

plots = lapply(todo, function(sel){
  plot_meta_alluvial(meta_dt[ik_status == sel & sample_type_short == "TBM"],
                     shown_filters = c("ik_status", "Phase", "RNA", "WGS", "miR"),
                     fill_var = "Phase") +
    labs(title = sel) +
    guides(fill = "none") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    scale_fill_manual(values = phase_colors)

})
# cowplot::plot_grid(plotlist = plots, ncol = 1)

p_ref = plot_meta_alluvial(meta_dt[ik_status == "IK6"],
                           shown_filters = c("ik_status", "Phase", "cell_of_origin", "sample_type_short", "RNA", "WGS", "miR"),
                           fill_var = "Phase")+
  scale_fill_manual(values = phase_colors) +
  guides(fill = guide_legend(ncol = 1))
leg = cowplot::get_legend(p_ref)
pg_alluvial = cowplot::plot_grid(nrow = 1, rel_widths = c(4, 1),
                   cowplot::plot_grid(plotlist = plots, ncol = 1),
                   leg
)
pg_alluvial

p_phase_and_seq = plot_meta_alluvial(meta_dt, shown_filters = c("Phase", "RNA", "WGS", "miR")) +
  labs(title = "Phase and Sequencing Done")
ggsave("phase_and_seq_alluvial.pdf", p_phase_and_seq, width = 10, height = 5.6 )

sample_codes_dt[sample_type_short == "TBM"]$sample_type


plot_upset_by_splice = function(meta_dt, plot_title, vars = c("ik_status", "Phase", "RNA", "WGS", "miR"), ...){
  # meta_dt.by_status = split(meta_dt[sample_type_short == "TBM"], meta_dt[sample_type_short == "TBM"]$ik_status)
  meta_dt.by_status = split(meta_dt, meta_dt$ik_status)
  memb_plots = lapply(names(meta_dt.by_status), function(nam){
    by_var = lapply(vars, function(var){
      xs = split(meta_dt.by_status[[nam]]$sample_id, meta_dt.by_status[[nam]][[var]])
      xs = xs[names(xs) != "-"]
      names(xs)[names(xs) == "present"] = var
      xs
    })
    memb = by_var[[1]]
    for(i in seq(2, length(by_var))){
      memb = c(memb, by_var[[i]])
    }
    seqsetvis::ssvFeatureUpset(memb, nsets = 10) + labs(title = nam)
  })
  cowplot::plot_grid(
    ggplot() + theme_void() + labs(title = plot_title),
    ncol = 1,
    rel_heights = c(1, 10),
    cowplot::plot_grid(plotlist = memb_plots, ...)
  )
}
pg_all = plot_upset_by_splice(meta_dt,
                              "All sample types",
                              vars = c("ik_status", "Phase", "WGS", "miR"),
                              nrow = 1)
pg_tbm = plot_upset_by_splice(meta_dt[sample_type_short == "TBM"],
                              "Primary Bone Marrow",
                              vars = c("ik_status", "Phase", "WGS", "miR"),
                              nrow = 1)
gc()
cowplot::plot_grid(nrow = 2,
  pg_all,
  pg_tbm
)
ggsave("splice_clusters_upset.simple.pdf", width = 12.3, height = 8.9)

fwrite(meta_dt, "sample_meta_dt.csv")
pg_all = plot_upset_by_splice(meta_dt, "All sample types",
                              nrow = 1)
pg_tbm = plot_upset_by_splice(meta_dt[sample_type_short == "TBM"], "Primary Bone Marrow",
                              nrow = 1)
gc()
cowplot::plot_grid(nrow = 2,
                   pg_all,
                   pg_tbm
)
ggsave("splice_clusters_upset.details.pdf", width = 12.3, height = 8.9)

lapply(by_var, unlist)
unlist(by_var)

tmp = meta_dt.by_status$IK6
tmp


plot_meta_alluvial(meta_dt[ik_status == "IK6"],
                   shown_filters = c("ik_status", "Phase", "cell_of_origin", "sample_type_short", "RNA", "WGS", "miR"),
                   fill_var = "Phase")

plot_meta_alluvial(meta_dt[ik_status == "early_termination"],
                   shown_filters = c("Phase", "cell_of_origin", "sample_type_short", "ik_status", "RNA", "WGS", "miR"),
                   fill_var = "Phase")

plot_meta_alluvial(meta_dt[ik_status == "normal_splicing"],
                   shown_filters = c("Phase", "cell_of_origin", "sample_type_short", "ik_status", "RNA", "WGS", "miR"),
                   fill_var = "Phase")
