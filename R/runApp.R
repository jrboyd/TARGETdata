
capitalize <- function(x) {
  x <- strsplit(x, " ")
  for (i in seq(along = x)) {
    substr(x[[i]], 1, 1) <- toupper(substr(x[[i]], 1, 1))
  }
  sapply(x, function(z) paste(z, collapse = " "))
}

#' Title
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
plot_meta_alluvial = function(meta_dt, shown_filters, fill_var = shown_filters[1]){
  shown_filters = as.list(shown_filters)
  max_axis = 5
  if(length(shown_filters) > max_axis){
    stop("too many filters for max axis")
  }
  if(length(shown_filters) < max_axis){
    shown_filters = c(shown_filters, list(rep(NULL, max_axis - length(shown_filters ))))
  }
  ggplot(meta_dt, aes_string(axis1 = shown_filters[[1]],
                             axis2 = shown_filters[[2]],
                             axis3 = shown_filters[[3]],
                             axis4 = shown_filters[[4]],
                             axis5 = shown_filters[[5]]
  )) +
    ggalluvial::geom_alluvium(aes_string(fill = fill_var), width = 1/5) +
    ggalluvial::geom_stratum(width = 1/8, fill = "black", color = "grey") +
    geom_label(stat = "stratum", aes(label = after_stat(stratum)), size = 1.8) +
    scale_x_discrete(limits = unlist(shown_filters), expand = c(.15, .15)) +
    scale_fill_brewer(type = "qual", palette = "Set1") +
    guides(fill = guide_legend(nrow = 3)) +
    labs(title = "Phase vs Cell of Origin")
}


#' TARGETdata.runApp
#'
#' @param et TARGETdata object
#'
#' @return
#' @export
#' @import shiny shinycssloaders
#' @rawNamespace import(shinyjs, except = runExample)
#' @examples
#' TARGETdata.runApp()
TARGETdata.runApp = function(){
  # Define UI for application that draws a histogram
  #### Setup ####
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

  active_filters = c("Phase", "sample_type", "cell_of_origin", "ik_status")
  stopifnot(active_filters %in% colnames(meta_dt))

  all_filters = active_filters
  filter = all_filters[3]
  make_filter_ui.single = function(filter, meta_dt){
    f_tab = table(meta_dt[[filter]])
    f_cnts = as.numeric(f_tab)
    f_names = names(f_tab)

    o = order(f_cnts, decreasing = TRUE)
    f_cnts = f_cnts[o]
    f_names = f_names[o]

    f_labels = paste0(f_names, " (", f_cnts, ")")
    ui_label = gsub("[_\\.]", " ", filter)
    ui_label = capitalize(ui_label)
    checkboxGroupInput(paste0("checkbox_", filter),
                       label = ui_label,
                       choiceNames = f_labels,
                       selected = f_names,
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
      tabPanel("Main",
               # Sidebar with a slider input for number of bins
               sidebarLayout(
                 sidebarPanel(
                   actionButton("btn_showNote", "Show Note"),
                   tags$span(),
                   make_filter_ui(active_filters, meta_dt)
                 ),
                 mainPanel(
                   tags$span()
                 )
               )
      ),
      tabPanel("Add Gene Set",
               tags$span()
      ),
      tabPanel("DE",
               tabsetPanel(
                 tabPanel("simple",
                          actionButton("btn_runDEfast", "Run DE"),
                          shinycssloaders::withSpinner(plotOutput("plot_volcano", width = "400px", height = "400px")),
                          shinycssloaders::withSpinner(plotOutput("plot_boxes", width = "400px", height = "400px"))

                 ),
                 tabPanel("DESeq2",
                          actionButton("btn_runDE", "Run DE"),
                          shinycssloaders::withSpinner(DT::dataTableOutput("dt_DE_res"))
                 )
               )
      )
    )
  )

  #### Server ####
  # Define server logic required to draw a histogram
  server <- function(input, output, session) {

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
    }, {
      meta_dt.filtered = meta_dt[
        Phase %in% input$checkbox_Phase &
          sample_type %in% input$checkbox_sample_type &
          cell_of_origin %in% input$checkbox_cell_of_origin &
          ik_status %in% input$checkbox_ik_status]
      meta_dt.sel(meta_dt.filtered)
      showNotification(nrow(meta_dt.sel()))
    })

    meta_dt.sel = reactiveVal(meta_dt)

    output$dt_DE_res = DT::renderDataTable({
      DT::datatable()
    })

    output$plot_boxes = renderPlot({
      ggplot()
    })

    output$plot_volcano = renderPlot({
      ggplot()
    })
  }

  shiny::runApp(list(ui = ui, server = server))
}
