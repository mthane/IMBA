
mod_aggregationTable_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      DTOutput(ns("aggregation_table"))
    ),
    fluidRow(
      
      downloadButton(ns('downloadAggregation'), 'Download')
    )
  )
}


mod_aggregationTable_server <- function(id,data) {
  moduleServer(
    id,
    function(input, output, session) {
      
      output$downloadAggregation <- downloadHandler(
        filename = function() {
          "aggregationData.csv"
        },
        content = function(file) {
          write.csv(data(), file, row.names = FALSE)
        }
      )
      
      output$aggregation_table <- renderDT({
        data()
      },options=list(scrollX=T))
    })}
