
#' name_of_module1 UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
#' @import shinyWidgets
#' @importFrom DT DTOutput 
#' 
#' 
#' 
#' 
#' 
mod_analysis_ui <- function(id, label = "Base directory") {
  ns <- NS(id)

    fluidRow(
      column(12,
             wellPanel(
               h3("Analyze datasets"),
               tags$div(
                 "
               You have to select the radius of the petri dish. This can be done by 
               using the preset option which would be the choice between 42.5mm and 69mm
               or using a manual radius.
               On top of that the length of the video has to be given in seconds which will 
               be used for filtering out tracks that are shorter then half of the video and the frame rate in FPS.
               ",
                 tags$hr()
               ),
               fluidRow(
                 column(2,
                        textInput(
                          ns('base_dir'),
                          label = "Base directory",
                          value = "C:"),
                        shinyDirButton(
                          ns("exp_dir"),
                          "Choose directory",
                          "Choose directory"
                        )
                 ),
                 column(2, 
                        textInput(ns('save_as'),
                                  label = "Folder name:",
                                  value = paste("experiment",Sys.Date(),sep="-")),
                        
                 ),
                 column(6,
                        uiOutput(ns("output_dir"))
                 )
               ),
               wellPanel(
                 fluidRow(
                   column(4,
                          h4("Video Options"),
                          radioButtons(ns("radiusAB"),
                                       "Select radius",
                                       choices = c("Preset"="button",
                                                   "Manual"="manual")),
                          uiOutput(ns("radius_choice")),
                          numericInput(ns("frame_rate"),label = "Frames per second", value = 16),
                          numericInput(ns("video_length"),label = "Video length",value = 180),
                          checkboxInput(ns("drop_contours"), "Drop contours", F)
                   ),
                   column(4,
                          h4("Head cast Options"),
                          numericInput(ns("threshold_hva_speed"),"Threshold for the head vector angular speed", 35),
                          numericInput(ns("threshold_tva_speed"),"Threshold for the tail vector angular speed", 45),
                          sliderInput(ns("hc_size"),label = "HC size",0,360,value=c(0,360))
                   )#,
                   #column(4,
                   #       h4("Step Options"),
                   #       numericInput(ns("threshold_f_step"),"Tail velocity threshold for forward step detection", 0.6),
                   #       numericInput(ns("threshold_b_step"),"Tail velocity threshold for backward step detection", -0.4)
                   #)
                 
                 )
               ),
               
               actionButton(ns("process_per_frame"),
                            "Analyze!")
             ),
             
        tags$hr(),        
        pre(id = ns("consoleOutput"),
            style = "overflow-y:scroll;
          max-height: 200px;
          #display:flex;
          flex-direction:column-reverse;"
        )
      
    )
  )
}



#' Upload server
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#' @import shiny
#' @import shinyFiles
#' @importFrom DT renderDT
#' @importFrom shinyjs html
# Define the server logic for a module
mod_analysis_server <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      ns <- session$ns
      
      ##### CHOOSING FOLDER #####
      
      updateDirChoose <- function(folder) {
        shinyFiles::shinyDirChoose(input = input,
                                   "exp_dir",
                                   roots = c(chosenFolder = folder))
      }
      
      observe({
        updateDirChoose(input$base_dir)
      })
      
      output$output_dir <- renderUI({
        req(input$exp_dir)
        ns <- session$ns
        path <- input$exp_dir[["path"]]
        
        textInput(ns("outputDir"),"Output directory... ",paste(path[seq(1,length(path)-1)], collapse = "/"))
      })
      
      # choose radius
      output$radius_choice <- renderUI({
        ns <- session$ns
        
        if(input$radiusAB=="manual"){
          numericInput(ns("radiusA"),
                       "Radius",
                       value = 42.5)
        }else{
          radioButtons(
            ns("radiusB"),
            "Radius",
            choices = c(42.5, 69),
            selected = 42.5
          )
        }
      })
      
      radius <- reactive({
        if(input$radiusAB=="manual"){
          input$radiusA
        }else{
          input$radiusB
        }
      })
      output_dir <- reactive(input$outputDir)
      
      
      metadata_dirs <- reactive({
        req(input$exp_dir)        
        message = "fetch metadata ..."
        shinyjs::html("consoleOutput", message)

        directory = paste(input$exp_dir[["path"]],
                          collapse = "/")
        directory = paste0(input$base_dir,
                           directory)
        dirs <- dir(directory, 
                    "metadata.txt",
                    recursive = TRUE)
        data.frame(dirs)
      })
      

      ##### REACTIVE VALUES #####
      
      # create data according to processing mode
      
      data_reactive <- reactiveValues(data = NULL)
      
      observeEvent(input$process_per_frame, {
        track_dir = paste0(input$base_dir, paste(input$exp_dir[["path"]], collapse = "/"))
        
        analysis_dir <- paste(output_dir(),
                              paste(input$save_as,"analysis",sep="-"), sep = "/")
        analysis_dir <- paste0(input$base_dir,analysis_dir)
       
        dir.create(analysis_dir)
        withCallingHandlers(
          data_reactive$data <-
            process_data_dt(
              input_dir = track_dir,
              output_dir = analysis_dir,
              save_as = input$save_as,
              radius = radius(),
              drop_contours = input$drop_contours,
              frame_rate = input$frame_rate,
              video_length = input$video_length,
              threshold_head_vector_angular_speed = input$threshold_hva_speed,
              threshold_tail_vector_angular_speed = input$threshold_tva_speed,
              min_hc_size = input$hc_size[1],
              max_hc_size = input$hc_size[2]#,
              #threshold_tail_speed_foward = input$threshold_f_step,
              #threshold_tail_speed_backward = input$threshold_b_step
            )
          
          ,
          # can use "warning" instead/on top of "message" to catch warnings too
          message = function(m) {
            shinyjs::html("consoleOutput", m$message, F)
          }
        )
      })
      

    }
  )
}