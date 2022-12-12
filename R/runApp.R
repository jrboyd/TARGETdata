
capitalize <- function(x) {
  x <- strsplit(x, " ")
  for (i in seq(along = x)) {
    substr(x[[i]], 1, 1) <- toupper(substr(x[[i]], 1, 1))
  }
  sapply(x, function(z) paste(z, collapse = " "))
}

parse_breaks = function(brks){
  out = suppressWarnings({
    as.numeric(strsplit(brks, ",")[[1]])
  })
  if(any(is.na(out)) | length(out) == 0)
    stop("invalid breaks input, must be comma delimited numbers")
  unique(sort(out))
}

#' plot_meta_alluvial
#'
#' @param meta_dt
#' @param shown_filters
#' @param fill_var
#'
#' @return
#' @export
#' @import ggalluvial
#'
#' @examples
#' clin_dt = load_clinical_data()
#' plot_meta_alluvial(clin_dt, shown_filters = c("Phase", "Cell.of.Origin", "Vital.Status"))
plot_meta_alluvial = function(meta_dt, shown_filters, fill_var = unlist(shown_filters[1])){
  shown_filters = as.list(shown_filters)
  max_axis = 8
  if(length(shown_filters) > max_axis){
    stop("too many filters for max axis")
  }
  len = length(shown_filters)
  for(i in seq_len(max_axis - length(shown_filters))){
    shown_filters = c(shown_filters, list(NULL))
  }
  meta_dt = meta_dt[, .N, c(unlist(shown_filters[!sapply(shown_filters, is.null)]))]

  ggplot(meta_dt, aes_string(axis1 = shown_filters[[1]],
                             axis2 = shown_filters[[2]],
                             axis3 = shown_filters[[3]],
                             axis4 = shown_filters[[4]],
                             axis5 = shown_filters[[5]],
                             axis6 = shown_filters[[6]],
                             axis7 = shown_filters[[7]],
                             axis8 = shown_filters[[8]],
                             y = "N"
  )) +
    geom_alluvium(aes_string(fill = fill_var), width = 1/5) +
    geom_stratum(width = 1/8, fill = "black", color = "grey") +
    geom_label(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
    scale_x_discrete(limits = unlist(shown_filters), expand = c(.15, .15)) +
    scale_fill_brewer(type = "qual", palette = "Set1") +
    guides(fill = guide_legend(nrow = 3)) +
    # labs(title = "Phase vs Cell of Origin") +
    theme(panel.background = element_blank(), panel.grid = element_blank())
}

plot_meta_upset = function(meta_dt, vars, plot_title = "", ...){
  by_var = lapply(vars, function(var){
    xs = split(meta_dt$sample_id, meta_dt[[var]])
    xs = xs[names(xs) != "-"]
    names(xs)[names(xs) == "present"] = var
    xs
  })
  memb = by_var[[1]]
  for(i in seq(2, length(by_var))){
    memb = c(memb, by_var[[i]])
  }
  seqsetvis::ssvFeatureUpset(memb, nsets = 10) + labs(title = plot_title)
}


#' TARGETdata.runApp.preload
#'
#' @param force
#'
#' @return
#' @export
#'
#' @examples
TARGETdata.runApp.preload = function(force = FALSE, cache_path = NULL){
  for(req_lib in c("data.table",
                   "ggalluvial",
                   "ggplot2",
                   "shiny")){
    if(!req_lib %in% (.packages())){
      library(req_lib, character.only = TRUE)
    }
  }

  .get_cache_path = function(obj_name){
    if(is.null(cache_path)){
      "/dev/null/not_a_file"
    }else{
      file.path(cache_path, paste0(obj_name, ".Rds"))
    }
  }
  .saveRDS = function(obj, file){
    if(!is.null(cache_path)){
      saveRDS(obj, file = file)
    }
  }

  .cached_load = function(obj_name, msg, FUN){
    if(!exists(obj_name) | force){
      c_file = .get_cache_path(obj_name)
      message(msg)
      if(file.exists(c_file)){
        obj = readRDS(c_file)
      }else{
        obj = FUN()
        .saveRDS(obj, file = c_file)
      }
      assign(obj_name, obj, envir = .GlobalEnv)
    }
    invisible(c_file)
  }

  .cached_load("mir_count_mat.full", "loading miR counts", load_miR_counts)
  .cached_load("mir_rpm_mat.full", "loading miR RPM", load_miR_RPM)
  .cached_load("clin_dt", "loading clin_dt", function(){
    clin_dt = load_clinical_data()
    data.table::setnames(clin_dt, "Cell.of.Origin", "cell_of_origin")
    clin_dt
  })
  .cached_load("ref_gr", "loading ref_gr", load_ref_gr)

  conv_dt = data.table::data.table(gene_id = ref_gr$gene_id, gene_name = ref_gr$gene_name)

  .cached_load("mrna_count_mat.full", "loading RNA counts", function(){
    mrna_count_mat.full = load_RNA_counts()
    if(any(rownames(mrna_count_mat.full) %in% names(ref_gr))){
      tmp = convert_matrix_2_data.table(mrna_count_mat.full)
      tmp = merge(tmp, conv_dt, by = "gene_id")
      tmp = tmp[, .(value = mean(value)), .(sample_id, gene_id = gene_name)]
      mrna_count_mat.full = convert_data.table_2_matrix(tmp)
    }
    mrna_count_mat.full
  })

  .cached_load("mrna_rpm_mat.full", "loading RNA RPM", function(){
    mrna_rpm_mat.full = load_RNA_RPM()
    if(any(rownames(mrna_rpm_mat.full) %in% names(ref_gr))){
      tmp = convert_matrix_2_data.table(mrna_rpm_mat.full)
      tmp = merge(tmp, conv_dt, by = "gene_id")
      tmp = tmp[, .(value = mean(value)), .(sample_id, gene_id = gene_name)]
      mrna_rpm_mat.full = convert_data.table_2_matrix(tmp)
    }
    mrna_rpm_mat.full
  })
  "done"
}


#' TARGETdata.runApp
#'
#' @param et TARGETdata object
#'
#' @return
#' @export
#' @import shiny shinycssloaders tippy
#' @rawNamespace import(shinyjs, except = runExample)
#' @examples
#' TARGETdata.runApp()
TARGETdata.runApp = function(){
  # Define UI for application that draws a histogram
  #### Setup ####
  TARGETdata.runApp.preload()
  #filter for clinical

  mrna_count_mat = filter_expression_to_valid(mrna_count_mat.full, clin_dt)
  mrna_rpm_mat = filter_expression_to_valid(mrna_rpm_mat.full, clin_dt)
  mir_count_mat = filter_expression_to_valid(mir_count_mat.full, clin_dt)
  mir_rpm_mat = filter_expression_to_valid(mir_rpm_mat.full, clin_dt)



  active_filters = c("Phase", "sample_type", "cell_of_origin", "ik_status")

  manual_selections = list(
    Phase = c("Phase_II_Validation", "Phase_II_Discovery"),
    sample_type = "Primary Blood Derived Cancer - Bone Marrow"
  )

  all_filters = active_filters
  filter = all_filters[2]
  make_filter_ui.single = function(filter, meta_dt){
    f_tab = table(meta_dt[[filter]])
    f_cnts = as.numeric(f_tab)
    f_names = names(f_tab)

    o = order(f_cnts, decreasing = TRUE)
    f_cnts = f_cnts[o]
    f_names = f_names[o]

    if(filter %in% names(manual_selections)){
      f_names.sel = manual_selections[[filter]]
    }else{
      f_names.sel = f_names
    }

    f_labels = paste0(f_names, " (", f_cnts, ")")
    ui_label = gsub("[_\\.]", " ", filter)
    ui_label = capitalize(ui_label)
    checkboxGroupInput(paste0("checkbox_", filter),
                       label = ui_label,
                       choiceNames = f_labels,
                       selected = f_names.sel,
                       choiceValues = f_names)
  }
  make_filter_ui = function(all_filters, meta_dt){
    lapply(all_filters, make_filter_ui.single, meta_dt = meta_dt)
  }

  #### UI ####
  ui <- fluidPage(
    theme = "bootstrap.css",
    # Application title
    shinyjs::useShinyjs(),
    titlePanel("TARGET data"),
    tabsetPanel(
      #### Panel: Selection  ####
      tabPanel("Sample Selection",
               # Sidebar with a slider input for number of bins
               sidebarLayout(
                 sidebarPanel(
                   actionButton("btn_showNote", "Show Note"),
                   tags$span(),
                   radioButtons("selDataType", "Seq Type", choices = c("RNAseq", "miRseq"), selected = "RNAseq"),
                   radioButtons("selMetaPlotType", "Plot Type", choices = c("Alluvial", "UpSet"), selected = "Alluvial"),
                   uiOutput("reactMetaFilter")
                 ),
                 mainPanel(
                   tags$span(),
                   shinycssloaders::withSpinner(plotOutput("plot_alluvial", width = "600px", height = "500px"))
                 )
               )
      ),
      #### Panel: Gene Query  ####
      tabPanel("Gene Query",
               sidebarLayout(
                 sidebarPanel(
                   selectizeInput(inputId = "txtGene",
                                  label = "Select Genes",
                                  choices = NULL,
                                  multiple = FALSE),
                   textInput(inputId = "txtBreaks",
                             label = "Expression Bin Breaks",
                             value = "-1,1"),
                   tags$span(id = "txtBreaksHelp", "breaks must be numeric and comma separated", style="color:red;"),
                   selectizeInput(inputId = "txtFacetVar",
                                  label = "Facet by",
                                  choices = c("none", active_filters),
                                  selected = "Phase",
                                  multiple = FALSE),
                   tippy::tippy_this(
                     "txtGene-container",
                     tooltip = "Delete and start typing gene names to see matching genes."
                   )
                 ),

                 # Show a plot of the generated distribution
                 mainPanel(
                   uiOutput("uiSurvival"),
                   # shinycssloaders::withSpinner(plotOutput("plotSurvivalGOI",
                   #                                         width = "550px", height = "310px")),
                   shinycssloaders::withSpinner(
                     plotOutput("plotExpressionGOI",
                                width = "550px", height = "310px")
                   )
                 )
               )
      ),
      tabPanel("Group Query",
               sidebarLayout(
                 sidebarPanel(
                   selectizeInput(inputId = "txtGroupVar", label = "Group by", choices = active_filters, selected = "ik_status", multiple = FALSE)
                 ),

                 # Show a plot of the generated distribution
                 mainPanel(
                   shinycssloaders::withSpinner(plotOutput("plotSurvivalGroup", width = "550px", height = "500"))
                 )
               )
      )#,
      #### Panel: Group Query  ####
      # tabPanel("Run DE",
      #          tabsetPanel(
      #            tabPanel("simple",
      #                     actionButton("btn_runDEfast", "Run DE"),
      #                     shinycssloaders::withSpinner(plotOutput("plot_volcano", width = "400px", height = "400px")),
      #                     shinycssloaders::withSpinner(plotOutput("plot_boxes", width = "400px", height = "400px"))
      #
      #            ),
      #            tabPanel("DESeq2",
      #                     actionButton("btn_runDE", "Run DE"),
      #                     shinycssloaders::withSpinner(DT::dataTableOutput("dt_DE_res"))
      #            )
      #          )
      # )
    )
  )

  #### Server ####
  # Define server logic required to draw a histogram
  server <- function(input, output, session) {

    rpm_mat = reactiveVal(mrna_rpm_mat)
    count_mat = reactiveVal(mrna_count_mat)
    meta_dt = reactiveVal(NULL)
    meta_dt.sel = reactiveVal(NULL)

    observeEvent({
      input$btn_showNote
    }, {
      showNotification("button pushed")
      print("button pushed")
      print(input$checkbox_Phase)
      print(input$checkbox_sample_type)
      print(input$checkbox_cell_of_origin)
      print(input$checkbox_ik_status)
    })


    observeEvent({
      #can't get this to work without hardcoding
      input$checkbox_Phase
      input$checkbox_sample_type
      input$checkbox_cell_of_origin
      input$checkbox_ik_status
      meta_dt()
    }, {
      stopifnot(active_filters %in% colnames(meta_dt()))
      meta_dt.filtered = meta_dt()[
        Phase %in% input$checkbox_Phase &
          sample_type %in% input$checkbox_sample_type &
          cell_of_origin %in% input$checkbox_cell_of_origin &
          ik_status %in% input$checkbox_ik_status]
      meta_dt.sel(meta_dt.filtered)
    })

    observeEvent({
      input$selDataType
    }, {
      #radioButtons("selDataType", "Seq Type", choices = c("RNAseq", "miRseq"), selected = "RNAseq"),
      req(input$selDataType)
      switch (input$selDataType,
              RNAseq = {
                rpm_mat(mrna_rpm_mat)
                count_mat(mrna_count_mat)
              },
              miRseq = {
                rpm_mat(mir_rpm_mat)
                count_mat(mir_count_mat)
              }
      )
    })


    observeEvent({
      rpm_mat()
    }, {
      suppressWarnings({
        meta_dt(make_meta_dt(rpm_mat(), clin_dt, extra_vars_included = c("cell_of_origin")))
        rpm_mat(filter_expression_to_valid(rpm_mat(), meta_dt()))
      })
    })

    output$reactMetaFilter = renderUI({
      req(meta_dt())
      make_filter_ui(active_filters, meta_dt())
    })

    output$dt_DE_res = DT::renderDataTable({
      DT::datatable()
    })

    output$plot_boxes = renderPlot({
      ggplot()
    })

    output$plot_volcano = renderPlot({
      ggplot()
    })
    # dpi = 150
    # output$plot_alluvial = renderPlot(width = dpi*4, height = dpi*6, res = dpi, {
    #### Alluvial Plot ####
    output$plot_alluvial = renderPlot({
      #radioButtons("selMetaPlotType", "Plot Type", choices = c("Alluvial", "UpSet"), selected = "Alluvial"),
      plot_meta_dt = req(meta_dt.sel())
      switch (input$selMetaPlotType,
              Alluvial = {
                plot_meta_alluvial(plot_meta_dt, active_filters) + theme(legend.position = "bottom")
              },
              UpSet = {
                plot_dt = plot_meta_dt[, c("sample_id", active_filters), with = FALSE]
                membs = lapply(active_filters, function(fil){
                  split(plot_dt$sample_id, plot_dt[[fil]])
                })
                plot_meta_upset(plot_meta_dt, "temp", vars = active_filters)
              },
              {
                ggplot() + labs(title = "bad plot type")
              }
      )

    })

    observeEvent({
      rpm_mat()
    }, {
      req(rpm_mat())
      message("+update txtGene")
      default = "CD34"
      all_genes = rownames(rpm_mat())
      # browser()
      updateSelectizeInput(session, 'txtGene', choices = all_genes, selected = ifelse(default %in% all_genes, default, all_genes[1]), server = TRUE)
      message("-update txtGene")
    })

    # uiOutput("uiSurvival"),
    # # shinycssloaders::withSpinner(plotOutput("plotSurvivalGOI",
    # #                                         width = "550px", height = "310px")),
    output$uiSurvival = renderUI({
      req(meta_dt.sel())
      req(input$txtFacetVar)
      if(input$txtFacetVar == "none"){
        nplot = 1
      }else{
        nplot = length(unique(meta_dt.sel()[[input$txtFacetVar]]))
      }
      showNotification(paste(nplot, " plots"))
      shinycssloaders::withSpinner(plotOutput("plotSurvivalGOI",
                                              width = ifelse(nplot == 1, "600px", "850"),
                                              height = paste0(300*nplot, "px")))
    })

    # withSpinner(plotOutput("plotSurvivalGOI", width = "550px", height = "310px")),
    # withSpinner(plotOutput("plotExpressionGOI", width = "550px", height = "310px")

    exp_breaks = reactiveVal(c(-1, 1))

    observeEvent({
      input$txtBreaks
    }, {
      brks = input$txtBreaks
      brk_res = tryCatch(
        expr = {
          brk_vals = parse_breaks(brks)
          message(paste(brk_vals, collapse = ", "))
          shinyjs::hide("txtBreaksHelp")
          exp_breaks(brk_vals)
        },
        error = function(e){
          shinyjs::show("txtBreaksHelp")
        }
      )
    })

    output$plotSurvivalGOI = renderPlot({
      req(meta_dt.sel())
      req(rpm_mat())
      req(input$txtGene)
      req(input$txtFacetVar)
      req(exp_breaks())
      message("plotSurvivalGOI")
      # res = run_survival(meta_dt.sel(), input$txtFacetVar)
      # res$plot

      goi = input$txtGene
      plot_rpm_mat = rpm_mat()
      plot_meta_dt = meta_dt.sel()
      facet_var = input$txtFacetVar
      brks = exp_breaks()

      cleanup_survival_plots = function(res, goi, xlim = NULL){
        p_surv = res$result$plot +
          guides(color = guide_legend(ncol = 1)) +
          # theme(legend.position = "bottom") +
          theme(legend.position = "right") +
          labs(title = paste0("Survival by ", goi, " expression bin"),
               color = paste(goi, "z-score"),
               fill = paste(goi, "z-score"))
        p_exp = res$expression_plot + guides(fill = "none")
        if(!is.null(xlim)){
          p_exp = p_exp + coord_cartesian(xlim = xlim)
        }
        pg = cowplot::plot_grid(p_exp, p_surv, rel_widths = c(1, 1.4))
        pg
      }

      if(facet_var == "none"){
        res = run_survival.goi_zscore(plot_rpm_mat, plot_meta_dt, goi, breaks = brks)
        pg = cleanup_survival_plots(res, goi)
        # p_surv = res$result$plot +
        #   guides(color = guide_legend(ncol = 1)) +
        #   theme(legend.position = "bottom") +
        #   labs(title = paste0("Survival by ", goi, " expression bin"),
        #        color = paste(goi, "z-score"),
        #        fill = paste(goi, "z-score"))
        # p_exp = res$expression_plot + guides(fill = "none")
        # pg = cowplot::plot_grid(p_exp, p_surv, rel_widths = c(1, 1.4))
      }else{
        todo = split(plot_meta_dt, plot_meta_dt[[facet_var]])
        z_mat = calc_zscore(plot_rpm_mat[goi, , drop = FALSE])
        plot_parts = lapply(names(todo), function(nam){
          meta_dt.split = todo[[nam]]
          res = run_survival.goi_zscore(z_mat, meta_dt.split, goi, apply_zscore = FALSE, breaks = brks)
          # res$plots
          pg_i = cleanup_survival_plots(res, goi, xlim = range(z_mat))
          pg_i
        })
        plot_labels = lapply(names(todo), function(nam){
          ggplot() +
            coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
            theme_void() +
            annotate("text", x= .5, y = .5, label = nam, hjust = .5, vjust = .5)
        })
        pg = cowplot::plot_grid(
          nrow = 1,
          rel_widths = c(1, 4),
          cowplot::plot_grid(plotlist = plot_labels, ncol = 1),
          cowplot::plot_grid(plotlist = plot_parts, ncol = 1)
        )
      }
      pg
    })

    output$plotSurvivalGroup = renderPlot({
      req(meta_dt.sel())
      req(input$txtGroupVar)
      message("plotSurvivalGOI")
      res = run_survival(meta_dt.sel(), input$txtGroupVar)
      res$plot
    })

    output$plotExpressionGOI = renderPlot({
      req(input$txtGene)
      req(input$txtFacetVar)
      req(rpm_mat())
      message("plotExpressionGOI")

      goi = input$txtGene
      facet_var = input$txtFacetVar
      plot_mat = rpm_mat()

      cnt_dt = data.table::data.table(sample_id = names(plot_mat[goi, ]), value = plot_mat[goi, ])
      cnt_dt = merge(cnt_dt, meta_dt.sel(), by = 'sample_id')
      cnt_dt[, lg_value := log10(value+.01)]
      # browser()
      if(facet_var == "none"){
        ggplot() +
          theme_void() +
          labs(title = "Nothing to plot when Facet Variable is none.")
      }else{
        ggplot(cnt_dt, aes_string(x = facet_var, y = "lg_value")) +
          geom_boxplot() +
          geom_jitter(width = .05) +
          scale_y_continuous(labels = function(x){
            ys = 10^(round(x-.01))
            # k = ys >= 1
            # ys[k] = round(ys[k])
            ys
          }, breaks = seq(-2, ceiling(max(cnt_dt$lg_value)))+.01) +
          # scale_color_manual(values = c("WT" = "black", "dF4" = "red")) +
          cowplot::theme_cowplot() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
          labs(title = paste(goi, "expression"), x = "group", y = "log10 RPM")
      }
    })
  }

  shiny::runApp(list(ui = ui, server = server))
}
