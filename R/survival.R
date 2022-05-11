#' run_survival_scan.phaseII
#'
#' @param goi_todo
#' @param n_todo
#' @param analysis_name
#'
#' @return
#' @export
#'
#' @examples
run_survival_scan.phaseII = function(goi_todo, analysis_name, n_todo = length(goi_todo)){

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
      surv_res_dt.done = lapply(already_done, fread) %>% rbindlist %>% unique
      surv_res_dt.done[gene_id %in% goi_todo.sel]

      message("found results for ", nrow(surv_res_dt.done))
      goi_todo.sel = setdiff(goi_todo.sel, surv_res_dt.done$gene_id)
    }else{
      surv_res_dt.done = NULL
    }

    message("running for ", length(goi_todo.sel))
    surv_res_dt = pbmcapply::pbmclapply(goi_todo.sel, function(goi){
      goi.name = ref_gr[goi]$gene_name

      res.disc = surv_fun.z(rpm_mat = exp_rpm_mat,
                            survival_info = meta_dt[Phase == "Phase_II_Discovery"],
                            goi = goi, goi.name = goi.name)

      res.vali = surv_fun.z(rpm_mat = exp_rpm_mat,
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
    surv_res_dt = surv_res_dt[k] %>% rbindlist

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
plot_survival_discovery_vs_validation = function(surv_res_dt){
  to_label = list(
    surv_res_dt[discovery < 10^-3.8]$gene_id,
    surv_res_dt[validation < 10^-3.8]$gene_id,
    surv_res_dt[validation < 10^-3 & discovery < 10^-3]$gene_id
  )
  to_label = lapply(to_label, function(x){
    if(length(x) > 10){
      x = sample(x, 10)
    }
    x
  }) %>% unlist

  subtitle = ifelse(n_todo < length(goi_todo.lncs),
                    paste(n_todo, "random genes from DE lncRNAs between discovery and validation"),
                    "genes from DE lncRNAs between discovery and validation")
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

#' plot_survival_goi
#'
#' @param goi
#'
#' @return
#' @export
#'
#' @examples
plot_survival_goi = function(goi){
  goi = subset(ref_gr, gene_name == "LINC00958")$gene_id
  goi.name = ref_gr[goi]$gene_name
  res.disc = surv_fun.z(rpm_mat = exp_rpm_mat,
                        survival_info = surv_dt.exp[Phase == "Phase_II_Discovery"],
                        goi = goi, goi.name = goi.name)

  res.vali = surv_fun.z(rpm_mat = exp_rpm_mat,
                        survival_info = surv_dt.exp[Phase == "Phase_II_Validation"],
                        goi = goi, goi.name = goi.name)

  p_disc = res.disc$plots
  p_vali = res.vali$plots

  rpm_dt = as.data.table(reshape2::melt(exp_rpm_mat[goi,]), keep.rownames = "sample_id")
  rpm_dt = merge(rpm_dt, surv_dt.exp, by = "sample_id")
  p_expression = ggplot(rpm_dt, aes(x = Phase, y = log2(value+.01))) +
    geom_violin() +
    labs(title = goi.name)
  pg = cowplot::plot_grid(
    cowplot::plot_grid(p_disc, p_vali, nrow = 1),
    p_expression
  )
  list(
    plot = pg,
    plot_parts = list(
      survival_discovery = res.disc,
      survival_validation = res.vali,
      expression = p_expression
    )
  )
}
