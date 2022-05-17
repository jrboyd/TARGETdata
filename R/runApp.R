
#' TARGETdata.runApp
#'
#' @param et TARGETdata object
#'
#' @return
#' @export
#' @import shiny shinycssloaders
#' @rawNamespace import(shinyjs, except = runExample)
#' @examples
#' ex_data = system.file("extdata/test_TARGETdata", package = "TARGETdata", mustWork = TRUE)
#' et = TARGETdata.load(ex_data)
#' TARGETdata.runApp(et)
TARGETdata.runApp = function(et){
  dataset_names = c(A = "fix", B = "this")
  gene_lists = list(
    All = rownames(et$norm_counts)
  )

  FACET_VAR = list(NONE = "none", SAMPLE_TYPE = "sample type", PAM50 = "PAM50")
  clean_list = function(l){
    tmp = unlist(l)
    names(tmp) = NULL
    tmp
  }

  # Define UI for application that draws a histogram
  ui <- fluidPage(
    theme = "bootstrap.css",
    # Application title
    useShinyjs(),
    titlePanel("TARGET data"),
    tabsetPanel(
      tabPanel("Main",
               # Sidebar with a slider input for number of bins
               sidebarLayout(
                 sidebarPanel(
                   tags$span()
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
               ui_point_selection(),
               tabsetPanel(
                 tabPanel("simple",
                          actionButton("btn_runDEfast", "Run DE"),
                          withSpinner(plotOutput("plot_volcano", width = "400px", height = "400px")),
                          withSpinner(plotOutput("plot_boxes", width = "400px", height = "400px"))

                 ),
                 tabPanel("DESeq2",
                          actionButton("btn_runDE", "Run DE"),
                          withSpinner(DT::dataTableOutput("dt_DE_res"))
                 )
               )
      )
    )
  )

  # Define server logic required to draw a histogram
  server <- function(input, output, session) {
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
