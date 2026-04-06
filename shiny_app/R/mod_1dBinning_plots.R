

mod_1dbinning_plots_ui <- function(id) {
  ns <- NS(id)
  tagList(
    wellPanel(
      
      
      
      plotlyOutput(ns("plot_line_binned"))
      
    ),
    wellPanel(
      wellPanel(
        checkboxInput(ns("use_upload"), "Use upload",F),
        fileInput(ns("binningFile"), "Choose CSV File",
                  multiple = FALSE,
                  accept = c("text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv"))
      ),
      DTOutput(ns("binned_dt")),
      
      downloadButton(ns('downloadBinning'), 'Download')
    )
  )
  
}

mod_1dbinning_plots_server <- function(id, upload, binnedData) {
  moduleServer(id,
               function(input, output, session) {
                 
                 xvariable <- reactive({
                   global_vars()$BINNING_VARTABLE %>% 
                     filter(variable==input$line_binned_variable_x)
                 })
                 
                 yvariable <- reactive({
                   global_vars()$BINNING_VARTABLE %>% 
                     filter(variable==input$line_binned_variable_y)
                 })
                 
                 output$description <- renderPrint({
                   text <- yvariable()$description%>%
                     str_replace("[1]","")%>%
                     str_replace('"',"")
                   cat(text)
                 })
                 
                 output$plot_line_binned <- renderPlotly({
                   binning_data_plot(
                     binnedData()$data,
                     binnedData()$x_var,
                     input$line_binned_variable_y,
                     as.numeric(upload()$radius),
                     global_vars()$BINNING_VARTABLE,
                     input$line_binned_npoints,
                     FALSE,
                     input$binning_mode,
                     input$use_polar
                   )
                   
                 })
                 
                 
                 
                 output$binned_dt <- renderDT({
                   binnedData()$data
                 },
                 options=list(scrollX=T))
                 
                 
                 output$downloadBinning <- downloadHandler(
                   filename = function() {
                     paste("binnedData", ".csv", sep = "")
                   },
                   content = function(file) {
                     write.csv(binnedData()$data, file, row.names = FALSE)
                   }
                 )
                 
               })
}
