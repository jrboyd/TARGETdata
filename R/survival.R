#' run_survival_scan.phaseII
#'
#' @param goi_todo Genes to analyze.  Should be rownames of rpm_mat.
#' @param rpm_mat RPM normalized gene expression.
#' @param meta_dt result of load_clinical_data()
#' @param analysis_name Prefix for analysis.
#' @param n_todo Should just be length of goi_todo
#' @param n_cores number of cores to use to parallelize processing
#' @param title_FUN function to run on gene_id to generate a title.
#'
#' @return data.table of pvalues for Discovery and Validation
#' @export
#'
#' @examples
#' mat_rpm = load_RNA_RPM()
#' clin_dt = load_clinical_data()
#' meta_dt = make_meta_dt(mat_rpm, clin_dt)
#' meta_dt = meta_dt[Phase %in% c("Phase_II_Discovery", "Phase_II_Validation")]
#' meta_dt[, rnk := frank(runif(.N), ties.method = "first"), .(Phase)]
#' meta_dt = meta_dt[rnk <= 10]
#' mat_rpm = filter_expression_to_valid(mat_rpm, meta_dt)
#' surv_res = run_survival_scan.phaseII(goi_todo = rownames(mat_rpm)[1:5],
#'     rpm_mat = mat_rpm,
#'     meta_dt = meta_dt,
#'     analysis_name = "test",
#'     n_cores = 1)
run_survival_scan.phaseII = function(goi_todo, rpm_mat, analysis_name, meta_dt = load_clinical_data(), n_todo = length(goi_todo), n_cores = getOption("mc.cores", 1), title_FUN = function(goi){goi}){
  phase_tab = table(meta_dt$Phase)
  if(is.na(phase_tab["Phase_II_Validation"]) | is.na(phase_tab["Phase_II_Discovery"])){
    stop("Phase_II_Validation and Phase_II_Discovery must be present in Phase of meta_dt")
  }
  if(n_todo > length(goi_todo)){
    n_todo = length(goi_todo)
  }
  surv_cache = res_file(paste0(analysis_name, ".survival_scan.", n_todo, ".csv"))

  if(file.exists(surv_cache)){
    message("loading cached results")
    surv_res_dt = fread(surv_cache)
  }else{
    already_done = dir(res_file(""), pattern = paste0(analysis_name, ".survival_scan..+.csv$"), full.names = TRUE)

    set.seed(0)
    goi_todo.sel =  sample(goi_todo, n_todo)

    if(length(already_done) > 0){
      surv_res_dt.done = unique(rbindlist(lapply(already_done, fread)))
      surv_res_dt.done[gene_id %in% goi_todo.sel]

      message("found results for ", nrow(surv_res_dt.done))
      goi_todo.sel = setdiff(goi_todo.sel, surv_res_dt.done$gene_id)
    }else{
      surv_res_dt.done = NULL
    }

    message("running survival for ", length(goi_todo.sel), " genes")
    surv_res_dt = pbmcapply::pbmclapply(goi_todo.sel, mc.cores = n_cores, function(goi){
      goi.name = title_FUN(goi)

      res.disc = run_survival.goi_zscore(rpm_mat = rpm_mat,
                                         survival_info = meta_dt[Phase == "Phase_II_Discovery"],
                                         goi = goi, goi.name = goi.name)

      res.vali = run_survival.goi_zscore(rpm_mat = rpm_mat,
                                         survival_info = meta_dt[Phase == "Phase_II_Validation"],
                                         goi = goi, goi.name = goi.name)
      data.table(
        gene_id = goi,
        gene_name = goi.name,
        discovery = res.disc$result$pval,
        validation = res.vali$result$pval
      )
    })
    k = sapply(surv_res_dt, function(x)class(x)[1]) == "data.table"
    surv_res_dt = rbindlist(surv_res_dt[k])

    if(!is.null(surv_res_dt.done)){
      surv_res_dt = rbind(
        surv_res_dt,
        surv_res_dt.done
      )
    }

    fwrite(surv_res_dt, surv_cache)
  }
  surv_res_dt
}

#' plot_survival_discovery_vs_validation
#'
#' @param surv_res_dt
#'
#' @return
#' @export
#'
#' @examples
plot_survival_discovery_vs_validation = function(surv_res_dt, n_labels = 10){
  to_label = list(
    surv_res_dt[discovery < 10^-3.8]$gene_id,
    surv_res_dt[validation < 10^-3.8]$gene_id,
    surv_res_dt[validation < 10^-3 & discovery < 10^-3]$gene_id
  )
  to_label = unlist(lapply(to_label, function(x){
    if(length(x) > n_labels){
      x = sample(x, n_labels)
    }
    x
  }))

  if(length(to_label) < n_labels){
    to_label = unique(c(
      to_label,
      sample(surv_res_dt$gene_id, min(nrow(surv_res_dt), n_labels - length(to_label)))
    ))
  }

  subtitle = "genes from DE lncRNAs between discovery and validation"
  p_surv.exp = ggplot(surv_res_dt, aes(x = -log10(validation), y = -log10(discovery))) +
    geom_point(size = .3) +
    ggrepel::geom_text_repel(data = surv_res_dt[gene_id %in% to_label], aes(label = gene_name), color = "gray40") +
    labs(
      title = "discovery vs validation significance for lncRNAs",
      subtitle = subtitle,
      x = "-log10(validation p value)",
      y = "-log10(discovery p value)") +
    theme(panel.background = element_blank(), panel.grid = element_blank(), axis.line = element_line())
  p_surv.exp
}

#' plot_survival_goi.phaseII
#'
#' @param goi
#'
#' @return
#' @export
#'
#' @examples
plot_survival_goi.phaseII = function(goi, rpm_mat, meta_dt, title_FUN = function(goi){goi}){
  goi.name = title_FUN(goi)
  res.disc = run_survival.goi_zscore(rpm_mat = rpm_mat,
                                     survival_info = meta_dt[Phase == "Phase_II_Discovery"],
                                     goi = goi, goi.name = goi.name)

  res.vali = run_survival.goi_zscore(rpm_mat = rpm_mat,
                                     survival_info = meta_dt[Phase == "Phase_II_Validation"],
                                     goi = goi, goi.name = goi.name)

  p_disc = res.disc$plots
  p_vali = res.vali$plots

  rpm_dt = as.data.table(reshape2::melt(rpm_mat[goi,]), keep.rownames = "sample_id")
  rpm_dt = merge(rpm_dt, meta_dt, by = "sample_id")
  p_expression = ggplot(rpm_dt, aes(x = Phase, y = log2(value+.01))) +
    geom_violin() +
    labs(title = goi.name)
  pg = cowplot::plot_grid(nrow = 1,
                          cowplot::plot_grid(p_disc, p_vali, ncol = 1),
                          p_expression
  )
  list(
    plot = pg,
    plot_parts = list(
      survival_discovery = res.disc,
      survival_validation = res.vali,
      expression = p_expression
    ),
    pvals = data.frame(goi = goi, discovery = res.disc$result$pval, validation = res.vali$result$pval)
  )
}




#' calc_zscore
#'
#' @param x
#' @param apply_log
#' @param pseudo
#'
#' @return z score values
#' @export
#'
#' @examples
#' calc_zscore(runif(100))
calc_zscore = function(x, apply_log = TRUE, pseudo = .01){
  if(apply_log){
    x = log10(x + pseudo)
  }
  z = (x - mean(x)) / sd(x)
  z[is.nan(z)] = 0
  z
}

#' run_survival.goi_zscore
#'
#' @param rpm_mat
#' @param survival_info
#' @param goi
#' @param goi.name
#' @param pseudo
#' @param breaks
#' @param bin_colors
#' @param return_data_only if TRUE, omit all plots and only return data
#'
#' @return
#' @export
#'
#' @examples
run_survival.goi_zscore = function(rpm_mat,
                                   survival_info,
                                   goi = "LEF1",
                                   goi.name = goi,
                                   pseudo = .01,
                                   breaks = c(-1,1),
                                   bin_colors = viridisLite::viridis(length(breaks)+1),
                                   return_data_only = FALSE,
                                   apply_zscore = TRUE){
  breaks = sort(breaks)

  common = intersect(colnames(rpm_mat), survival_info$sample_id)
  message("Using ", length(common), " of ", ncol(rpm_mat), " matrix samples")
  message("Using ", length(common), " of ", nrow(survival_info), " survival samples")

  rpm_mat = rpm_mat[, common, drop = FALSE]

  exp_sub<- data.frame(sample_id = colnames(rpm_mat), counts=rpm_mat[goi,])
  if(apply_zscore){
    exp_sub$z = calc_zscore(exp_sub$counts, apply_log = TRUE, pseudo = pseudo)
  }else{
    exp_sub$z = exp_sub$counts
  }
  breaks.full = c(-Inf, breaks, Inf)
  rect_dt = rbindlist(lapply(seq_along(breaks.full)[-1], function(i){
    data.table(xmin = breaks.full[i-1], xmax = breaks.full[i], group = as.character(i-1))
  }))

  if(!return_data_only){
    p_exp_hist = ggplot(exp_sub) +
      geom_rect(data = rect_dt,
                aes(xmin = xmin,
                    xmax = xmax,
                    ymin = 0,
                    ymax = Inf,
                    fill = group),
                alpha = .3) +
      geom_histogram(aes(x = z), bins = 30) +
      geom_vline(xintercept = breaks, color = "red", size = 1.5) +
      labs(title = paste("Expression distribution for", goi.name), x= "z-score of log10 expression") +
      scale_fill_manual(values = bin_colors) +
      scale_y_continuous(expand = expansion(0, 0)) +
      theme(panel.background = element_blank(), panel.grid = element_blank())
  }else{
    p_exp_hist = ggplot() + labs(title = "return_data_only = TRUE")
  }

  exp_sub$group = cut(exp_sub$z, breaks = c(-Inf, breaks, Inf))

  #merge exp_sub to patient data
  exp_sub.surv = merge(exp_sub, survival_info, by = "sample_id")
  # run_survival
  exp_sub.surv = exp_sub.surv[order(exp_sub.surv$group),]

  #TCGA biolink does not use factor order for some reason, this forces colors to be alphabetical to match TCGA biolinks
  names(bin_colors) = levels(exp_sub.surv$group)
  group_present = table(exp_sub.surv$group) > 0
  bin_colors = bin_colors[group_present]
  bin_colors = bin_colors[order(names(bin_colors))]
  names(bin_colors) = NULL

  if(sum(group_present) < 2){
    res.goi = list(
      plot = ggplot()+theme_void()+labs(title = 'not enough expression distribution to test'),
      pval = 1
    )
  }else{
    res.goi<- run_survival(exp_sub.surv, "group", legend_title = paste(goi.name), color = bin_colors)
  }
  if(!return_data_only){
    pg = cowplot::plot_grid(p_exp_hist, res.goi$plot)
  }else{
    pg = p_exp_hist
  }
  list(plots = pg, result = res.goi, expression_plot = p_exp_hist, expression_data = exp_sub)
}

#' run_survival
#'
#' @param surv_dt
#' @param group_var
#' @param legend_title
#' @param color
#'
#' @return
#' @export
#'
#' @examples
run_survival = function(surv_dt, group_var, legend_title = group_var, color = NULL){
  res = TCGAbiolinks::TCGAanalyze_survival(as.data.frame(surv_dt),
                                           group_var, legend = legend_title,
                                           height = 10,
                                           width=10,
                                           filename = NULL,
                                           color = color )
  n_layers = length(res$plot$layers)
  pval = res$plot$layers[[n_layers]]$aes_params$label
  if(is.null(pval)){
    pval = res$plot$layers[[n_layers]]$geom_params$label
  }
  pval = as.numeric(sub("p [=<>] ", "", pval))
  res$pval = pval
  res
}
