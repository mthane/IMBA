#' filtering UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_filtering_ui <- function(id){
  ns <- NS(id)
  tagList(
    checkboxInput(ns("showFilter"),"Show filters",T),
    uiOutput(ns("filtering"))
  )
}

#' filtering Server Functions
#'
#' @noRd 
mod_filtering_server <- function(id,upload){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    
    output$filtering <- renderUI({
      
      ns <- session$ns
      if(input$showFilter & !is.null(upload()$video_length)){
        wellPanel(
          fluidRow(
            column(2,
                   radioButtons(ns("direction"),
                                      "Run direction",
                                      choices = c("None"= "none",
                                                  "Forwards" = "forwards",
                                                  "Backwards" = "backwards",
                                                  "Both" = "both"
                                                  )
                                      )
                   ),
            column(5,
                   sliderInput(
                     ns("time"),
                     "Time [s]:",
                     -180,
                     upload()$video_length,
                     c(0, upload()$video_length),
                     step = 5),
                   sliderInput(
                     ns("hc_size"),
                     "Absolute HC Angle [°]:",
                     0,
                     360,
                     c(20, 360),
                     step = 10)
            ),
            column(5,
                   sliderInput(
                     ns("bearing_angle"),
                     "Absolute bearing angle [°]:",
                     0,
                     180,
                     c(0, 180),
                     step = 10
                   ),
                   sliderInput(
                     ns("odor_distance"),
                     "Distance to odor:",
                     0,
                     as.numeric(upload()$radius*2),
                     c(0, as.numeric(upload()$radius*2)),
                     step = 1
                   )
            )
          )
        )
      }
    })
    
    reactive(
      list(time= input$time,
           hc_size= input$hc_size,
           bearing_angle= input$bearing_angle,
           odor_distance= input$odor_distance,
           direction = input$direction
      )
    )
  })
}

## To be copied in the UI
# mod_filtering_ui("filtering_ui_1")

## To be copied in the server
# mod_filtering_server("filtering_ui_1")
