
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

todo = c("RNA_and_miR", "RNA_and_WGS", "RNA_and_miR_and_WGS")
sel = todo[2]

plot_sample_alluvial_coverage = function(sel){
  ax1 = "RNA"
  ax2 = "miR"
  ax3 = "WGS"
  ggplot(cov_dt, aes_string(axis1 = ax1, axis2 = ax2, axis3 = ax3)) +
    ggalluvial::geom_alluvium(aes_string(fill = sel), width = 1/5) +
    ggalluvial::geom_stratum(width = 1/8, fill = "black", color = "grey") +
    geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_x_discrete(limits = c(ax1, ax2, ax3), expand = c(.15, .15)) +
    scale_fill_manual(values = c("FALSE" = "gray", "TRUE" = "forestgreen")) +
    # scale_fill_brewer(type = "qual", palette = "Set1") +
    guides(fill = "none") +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(size = 10),
      plot.subtitle = element_text(size = 8),
      plot.caption = element_text(family = "mono", size = 8)
    ) +
    labs(title = sel, subtitle = "TARGET sequencing at sample resolution", caption = paste(sum(cov_dt[[sel]]), "shared"))

}
plots = lapply(c("RNA_and_miR", "RNA_and_WGS", "RNA_and_miR_and_WGS"), plot_sample_alluvial_coverage)

gc()
pg = cowplot::plot_grid(plotlist = plots, nrow = 1)
ggsave("describe_sample_sequencing.pdf", pg, width = 9, height = 3.7)

#### patient_id coverage ####

cov_dt.pat = rbind(
  data.table(seq = "WGS", patient_id = unique(convert_sample_id_2_patient_id(ids_wgs))),
  data.table(seq = "miR", patient_id = unique(convert_sample_id_2_patient_id(ids_mir))),
  data.table(seq = "RNA", patient_id = unique(convert_sample_id_2_patient_id(ids_rna)))
)
cov_dt.pat$present = "present"
cov_dt.pat = dcast(cov_dt.pat, patient_id~seq, value.var = "present", fill = "-")

ax1 = "RNA"
ax2 = "miR"
ax3 = "WGS"

universal_ids = cov_dt.pat[RNA != "-" & WGS != "-" & miR != "-"]$patient_id
cov_dt.pat[, universal := patient_id %in% universal_ids]
cov_dt.pat[, presence := sum(RNA == "present") + sum(WGS == "present") + sum(miR == "present"), .(patient_id)]
cov_dt.pat$presence = factor(cov_dt.pat$presence)
cov_dt.pat[, RNA_and_miR := RNA == "present" & miR == "present"]
cov_dt.pat[, RNA_and_WGS := RNA == "present" & WGS == "present"]
cov_dt.pat[, RNA_and_miR_and_WGS := RNA == "present" & WGS == "present"  & miR == "present"]

todo = c("RNA_and_miR", "RNA_and_WGS", "RNA_and_miR_and_WGS")
sel = todo[2]
plot_patient_alluvial_coverage = function(sel){
  ax1 = "RNA"
  ax2 = "miR"
  ax3 = "WGS"
  ggplot(cov_dt.pat, aes_string(axis1 = ax1, axis2 = ax2, axis3 = ax3)) +
    ggalluvial::geom_alluvium(aes_string(fill = sel), width = 1/5) +
    ggalluvial::geom_stratum(width = 1/8, fill = "black", color = "grey") +
    geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_x_discrete(limits = c(ax1, ax2, ax3), expand = c(.15, .15)) +
    scale_fill_manual(values = c("FALSE" = "gray", "TRUE" = "forestgreen")) +
    # scale_fill_brewer(type = "qual", palette = "Set1") +
    guides(fill = "none") +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(size = 10),
      plot.subtitle = element_text(size = 8),
      plot.caption = element_text(family = "mono", size = 8)
    ) +
    labs(title = sel, subtitle = "TARGET sequencing at patient resolution", caption = paste(sum(cov_dt.pat[[sel]]), "shared"))
}
plots = lapply(c("RNA_and_miR", "RNA_and_WGS", "RNA_and_miR_and_WGS"), plot_patient_alluvial_coverage)

gc()
pg = cowplot::plot_grid(plotlist = plots, nrow = 1)
ggsave("describe_patient_sequencing.pdf", pg, width = 9, height = 3.7)

#### Full clinical ####
full_clin_dt = clin_dt[, .(patient_id, Cell.of.Origin, Phase)]
full_clin_dt = merge(full_clin_dt, cov_dt.pat[, .(patient_id, RNA, WGS, miR)], by = "patient_id", all.x = TRUE, )
full_clin_dt[is.na(RNA), RNA := "-"]
full_clin_dt[is.na(WGS), WGS := "-"]
full_clin_dt[is.na(miR), miR := "-"]
ax1 = "Phase"
ax2 = "Cell.of.Origin"
ax3 = "RNA"
ax4 = "miR"
ax5 = "WGS"
sel = "Phase"

full_clin_dt.alluv = full_clin_dt[, .N, c(ax1, ax2, ax3, ax4, ax5)]

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
  axis.text = element_text(size = 8),
  plot.caption = element_text(family = "mono", size = 8)
))

p_phase_and_origin = ggplot(full_clin_dt.alluv, aes_string(axis1 = ax1, axis2 = ax2, y = "N")) +
  ggalluvial::geom_alluvium(aes_string(fill = sel), width = 1/5) +
  ggalluvial::geom_stratum(width = 1/8, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum)), size = 1.8) +
  scale_x_discrete(limits = c(ax1, ax2), expand = c(.15, .15)) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  guides(fill = guide_legend(nrow = 3)) +
  labs(title = "Phase vs Cell of Origin")
p_phase_and_origin
ggsave("full_phase_and_origin.png", p_phase_and_origin, width = 5, height = 3.6, dpi = "print")

full_clin_dt$presence = NULL
full_clin_dt[, presence := sum(RNA == "present") + sum(WGS == "present") + sum(miR == "present"), .(patient_id)]
# full_clin_dt$presence = factor(full_clin_dt$presence)

full_clin_dt[, RNA_and_WGS := ifelse(RNA == "present" & WGS == "present", "present", "-")]
full_clin_dt[, RNA_WGS_and_miR := ifelse(RNA == "present" & WGS == "present" & miR == "present", "present", "-")]
full_clin_dt[, RNA_and_miR := ifelse(RNA == "present" & miR == "present", "present", "-")]
full_clin_dt[, WGS_and_miR := ifelse(WGS == "present" & miR == "present", "present", "-")]

full_clin_dt.alluv_presence = full_clin_dt[, .N, c(ax1, ax2, ax3, ax4, ax5, "presence", "RNA_and_WGS", "RNA_WGS_and_miR", "RNA_and_miR", "WGS_and_miR")]
todo_presence = c("RNA", "miR", "WGS", "RNA_and_WGS", "RNA_and_miR", "WGS_and_miR", "RNA_WGS_and_miR")
stopifnot(todo_presence %in% colnames(full_clin_dt.alluv_presence))
sel = todo_presence[1]
seq_plots = lapply(todo_presence, function(sel){
  ggplot(full_clin_dt.alluv_presence, aes_string(axis1 = ax1, axis2 = ax2, axis3 = ax3, axis4 = ax4, axis5 = ax5, y = "N")) +
    ggalluvial::geom_alluvium(aes_string(fill = sel), width = 1/5) +
    ggalluvial::geom_stratum(width = 1/8, fill = "black", color = "grey") +
    geom_label(stat = "stratum", aes(label = after_stat(stratum)), size = 1.8) +
    scale_x_discrete(limits = c(ax1, ax2, ax3, ax4, ax5), expand = c(.15, .15)) +
    scale_fill_manual(values = c("-" = "gray", "present" = "red")) +
    guides(fill = guide_legend(nrow = 3)) +
    labs(title = sel) +
    facet_wrap(~Phase, scales = "free_y")
})
pg_seq_plots = cowplot::plot_grid(plotlist = seq_plots)
ggsave("full_seq_by_phase.png", pg_seq_plots, width = 32, height = 16, dpi = 60)

full_clin_dt
ggplot(full_clin_dt.alluv_presence, aes(x = Phase, y = "N"))

tmp = copy(full_clin_dt)
tmp$presence = NULL
tmp = melt(tmp, id.vars = c("Phase", "Cell.of.Origin", "patient_id"), value.name = "N")

tmp[, .N, .(Phase, variable, N)][Phase == "Dicentric"]
tmp[, .N, .(Phase, variable, N)][variable == "RNA" & Phase == "Phase_III"]
tmp[, .N, .(Phase, variable, N)][variable == "RNA_WGS_and_miR" & Phase == "Phase_III"]
tmp[, .N, .(Phase, variable, N)][variable == "RNA" & Phase == "Dicentric"]
tmp$variable = factor(tmp$variable, levels = c("RNA", "miR", "WGS", "RNA_and_WGS", "RNA_and_miR", "WGS_and_miR", "RNA_WGS_and_miR"))
p_seq_barplots = ggplot(tmp, aes(x = Phase, fill = N)) +
  geom_bar() +
  facet_wrap(~variable, ncol = 4) +
  scale_fill_manual(values = c("-" = "gray", "present" = "red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(fill = 'sequencing done')
ggsave("full_seq_by_phase_barplot.png", p_seq_barplots, width = 8.6, height = 5.8)

# p_phase_and_origin_sequencing = ggplot(full_clin_dt.alluv_presence, aes_string(axis1 = ax1, axis2 = ax2, axis3 = ax3, axis4 = ax4, axis5 = ax5, y = "N")) +
#   ggalluvial::geom_alluvium(aes_string(fill = "presence"), width = 1/5) +
#   ggalluvial::geom_stratum(width = 1/8, fill = "black", color = "grey") +
#   geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
#   scale_x_discrete(limits = c(ax1, ax2, ax3, ax4, ax5), expand = c(.15, .15)) +
#   # scale_fill_brewer(type = "qual", palette = "Set1") +
#   guides(fill = guide_legend(nrow = 3)) +
#   labs(title = "Phase vs Cell of Origin") +
#   facet_wrap(~Phase, scales = "free_y")
# p_phase_and_origin_sequencing
# ggsave("full_phase_and_origin.png", p_phase_and_origin, width = 8, height = 8)



ggplot(full_clin_dt, aes_string(axis1 = ax1, axis2 = ax2, axis3 = ax3, axis4 = ax4, axis5 = ax5)) +
  ggalluvial::geom_alluvium(aes_string(fill = "presence"), width = 1/5) +
  ggalluvial::geom_stratum(width = 1/8, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c(ax1, ax2, ax3, ax4, ax5), expand = c(.15, .15)) +
  scale_fill_brewer(type = "qual", palette = "Dark2") +
  # guides(fill = "none") +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(size = 10),
    plot.subtitle = element_text(size = 8),
    plot.caption = element_text(family = "mono", size = 8)
  ) +
  labs(title = "Phase vs Cell of Origin") +
  facet_wrap(~Phase, scales = "free_y")


sel = "Cell.of.Origin"
ggplot(full_clin_dt, aes_string(axis1 = ax1, axis2 = ax2, axis3 = ax3, axis4 = ax4, axis5 = ax5)) +
  ggalluvial::geom_alluvium(aes_string(fill = sel), width = 1/5) +
  ggalluvial::geom_stratum(width = 1/8, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c(ax1, ax2, ax3, ax4, ax5), expand = c(.15, .15)) +
  # scale_fill_manual(values = c("FALSE" = "gray", "TRUE" = "forestgreen")) +
  scale_fill_brewer(type = "qual", palette = "Dark2") +
  guides(fill = "none") +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(size = 10),
    plot.subtitle = element_text(size = 8),
    plot.caption = element_text(family = "mono", size = 8)
  ) +
  labs(title = sel, subtitle = "TARGET sequencing at patient resolution", caption = paste(sum(cov_dt.pat[[sel]]), "shared"))


plot_full_clinical_alluvial = function(sel){
  ax1 = "RNA"
  ax2 = "miR"
  ax3 = "WGS"
  ggplot(cov_dt.pat, aes_string(axis1 = ax1, axis2 = ax2, axis3 = ax3)) +
    ggalluvial::geom_alluvium(aes_string(fill = sel), width = 1/5) +
    ggalluvial::geom_stratum(width = 1/8, fill = "black", color = "grey") +
    geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_x_discrete(limits = c(ax1, ax2, ax3), expand = c(.15, .15)) +
    scale_fill_manual(values = c("FALSE" = "gray", "TRUE" = "forestgreen")) +
    # scale_fill_brewer(type = "qual", palette = "Set1") +
    guides(fill = "none") +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(size = 10),
      plot.subtitle = element_text(size = 8),
      plot.caption = element_text(family = "mono", size = 8)
    ) +
    labs(title = sel, subtitle = "TARGET sequencing at patient resolution", caption = paste(sum(cov_dt.pat[[sel]]), "shared"))
}


#### Cross sequencing groups in detail ####
cov_dt
cov_dt[, sample_code := convert_sample_id_2_sample_code(sample_id)]
merge(type_dt, cov_dt, by = "sample_code")


