#' aggregation UI Function
#'
#' @description Shiny module for creating grouped data based on the time series data.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_aggregation_ui <- function(id){
  ns <- NS(id)
  tagList(
    uiOutput(ns("uploadAggregation")),
    checkboxInput(ns("use_upload"), "Use upload",F),
    checkboxInput(ns("rm_outliers"), "Remove Outliers",F)

  )
}


#' aggregation Server Functions
#'
#' @noRd 
mod_aggregation_server <- function(id,upload,threshold){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    output$uploadAggregation <- renderUI({
      

   
      if(input$use_upload){
        wellPanel(
          fileInput(ns("aggregationFile"), "Choose CSV File",
                    multiple = FALSE,
                    accept = c("text/csv",
                               "text/comma-separated-values,text/plain",
                               ".csv"))
        )
        
      }else{
        
        if(is.null(upload()$data)){
          tagList(
            span("There is no analysis data uploaded.
                 Upload a valid .Rds file or upload a summary file by selecting the checkbox 'Use upload'! ",
                 style="color:red")
          )
        }else{
        wellPanel(
          fluidRow(
            column(4,
                   
                   radioButtons(
                     ns("grouping_mode"),
                     label = "Grouping",
                     choices =
                       c("Dish",
                         "Individual")
                   )
            ),
            column(8,
                   mod_filtering_ui(ns("aggregationFilter"))
                   
            )
          ),
          fluidRow(
            actionButton(ns("summarise"), "Update aggregation")
          ),
          
          pre(id = ns("consoleOutput"))
        )
        
      }
      }

    })
    
    filters <- mod_filtering_server("aggregationFilter",upload)
    
    data <- eventReactive(input$summarise, {
      ns <- session$ns
      message("garbage collection")
      gc()
      
      if (input$grouping_mode == "Dish") {
        grouped_by = "trial"
      }
      else
      {
        grouped_by = "id"
      }
      tic()
      
      message = paste("summarise data... ", "grouped by", grouped_by , sep = " ")
      shinyjs::html("consoleOutput", message)
      fdata <- upload()$data%>%as.data.table()
      filters <- filters()
      data <- create_summarized_analysis_dt(
        fdata,
        grouped_by,
        filters$time,
        filters$odor_distance,
        filters$hc_size,
        filters$bearing_angle,
        filters$spinepoint_y_6_conv, #new after publication
        filters$spinepoint_x_6_conv, #new after publication
        filters$abs_y_angle, #new after publication
        upload()$radius,
        upload()$video_length,
        upload()$frame_rate,
        filters$direction,
        threshold()
      )

      it <- toc()
      exectime <- it$toc - it$tic
      message = paste("summarising data: ",round(exectime) , "seconds elapsed",sep=" ")
      shinyjs::html("consoleOutput", message)
      data%>%as.data.frame()
    })
    
    
    
    aggregationFile <- reactive({
      req(input$aggregationFile)
      fread(input$aggregationFile$datapath)
    })
    
    aggregationData <- reactive({
      remove_outliers <- function(x, na.rm = TRUE, ...) {
        qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
        H <- 1.5 * IQR(x, na.rm = na.rm)
        y <- x
        y[x < (qnt[1] - H)] <- NA
        y[x > (qnt[2] + H)] <- NA
        y
      }
      if(input$use_upload){
        df <- aggregationFile()
      }
      else{
        df <- data()
      }
      if(input$rm_outliers){
        df <- df%>%mutate_if(is.numeric,remove_outliers)
      }

      df
    })
    aggregationData
  })
}

## To be copied in the UI
# mod_aggregation_ui("aggregation_ui_1")

## To be copied in the server
# mod_aggregation_server("aggregation_ui_1")
