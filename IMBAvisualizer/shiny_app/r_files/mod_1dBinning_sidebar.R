
mod_1dbinning_sidebar_ui <- function(id) {
  ns <- NS(id)
  fluidRow(
    
    pickerInput(ns("line_binned_variable_y"),
                "Choose variable y",
                choices = global_vars()$BINNING_VARIABLES),
    sliderInput(ns("line_binned_npoints"),
                "Convolution interval ",
                1,
                50,
                11),
    checkboxInput(ns("use_polar"),"Use polar coordinates (For bearing angle)",value = T),
    wellPanel(textOutput(ns("description"))),
    radioButtons(
      ns("xvar"),
      "Choose attribute x",
      choices = c(
        "Time" =
          "frame",
        "Bearing angle" =
          "bearing_angle",
        "Distance to odor" =
          "distance_to_odor",
        "Time before entering odour half" = #new after publication
          "time_to_trans_to",
        "Time before leaving odour half" = #new after publication
          "time_to_trans_from",
        "Y coordinate" = #new after publication
          "spinepoint_y_6_conv",
        "X coordinate" = #new after publication
          "spinepoint_x_6_conv",
        "Direction relative to Y" =
          "y_angle"
      )),
    
    uiOutput(ns("binwidth")),
    mod_filtering_ui(ns("binningFilter")),
    radioButtons(
      ns("binning_mode"),
      "Binning Mode",
      choices = c("Population" = "all",
                  "Dish" = "trial",
                  "Individual" = "id")
    ),      
    actionButton(ns("binning_button"), "Update binning")
    
  )
  
}




mod_1dbinning_sidebar_server <- function(id, upload) {
  moduleServer(id,
               function(input, output, session) {
                 ns <- session$ns
                 
                 filters = mod_filtering_server("binningFilter",upload)
                 
                 output$binwidth <- renderUI({
                   
                   if(input$xvar=="frame"){
                     limits = c(1,upload()$video_length/4,5) #upload()$video_length
                   }
                   
                   if(input$xvar=="bearing_angle"){
                     limits = c(5,30,10)
                   }
                   
                   if(input$xvar=="distance_to_odor"){
                     limits = c(2,20,4)
                   }
                   
                   if(input$xvar=="time_to_trans_to"){ #new after publication
                     limits = c(1,upload()$video_length/4,5)
                   }
                   
                   if(input$xvar=="time_to_trans_from"){ #new after publication
                     limits = c(1,upload()$video_length/4,5)
                   }
                   
                   if(input$xvar=="spinepoint_y_6_conv"){ #new after publication
                     limits = c(2,20,4)
                   }
                   
                   if(input$xvar=="spinepoint_x_6_conv"){ #new after publication
                     limits = c(2,20,4)
                   }
                   
                   if(input$xvar=="y_angle"){ #new after publication
                     limits = c(5,30,10)
                   }
                   
                   sliderInput(
                     ns("width"),
                     "Bin width",
                     limits[1],
                     limits[2],
                     limits[3],
                     limits[4],
                     limits[5]#,
                     #limits[6],
                     #limits[7],
                     #limits[8],
                     #limits[9]
                   )
                   
                 })
                 
                 binnedData <- eventReactive(input$binning_button, {
                   x_var <- input$xvar
                   
                   
                   withCallingHandlers({
                     data <- upload()$data
                     data$group_condition <- paste(data$group,data$condition,sep="-")
                     groups = unique(data$group_condition)
                     dfs = list()
                     filters <- filters()
                     i = 1
                       suppressWarnings({
                         
                         if(input$binning_mode =="all"){
                           df <-
                             binning_data(
                               df = data,
                               variable =  x_var,
                               width = input$width,
                               frame_interval = filters$time,
                               distance_to_odor_interval = filters$odor_distance,
                               Abs_HC_Angle_interval =  filters$hc_size,
                               Abs_bearing_angle = filters$bearing_angle,
                               spinepoint_y_6_interval = filters$spinepoint_y_6_conv,
                               spinepoint_x_6_interval = filters$spinepoint_x_6_conv,
                               abs_y_angle_interval = filters$abs_y_angle,
                               radius = upload()$radius,
                               frame_rate=upload()$frame_rate,
                               filters$direction
                             )
                           
                         }else{
                           df <-
                             binning_data_grouped(
                               df = data,
                               variable =  x_var,
                               width = input$width,
                               frame_interval = filters$time,
                               distance_to_odor_interval = filters$odor_distance,
                               Abs_HC_Angle_interval =  filters$hc_size,
                               Abs_bearing_angle = filters$bearing_angle,
                               spinepoint_y_6_interval = filters$spinepoint_y_6_conv,
                               spinepoint_x_6_interval = filters$spinepoint_x_6_conv,
                               abs_y_angle_interval = filters$abs_y_angle,
                               radius = upload()$radius,
                               frame_rate=upload()$frame_rate,
                               grouping=input$binning_mode,
                               filters$direction
                             )
                           
                         }
                         
                       })
                     
                     list(data = df,
                          x_var=x_var)
                   },
                   
                   # can use "warning" instead/on top of "message" to catch warnings too
                   message = function(m) {
                     shinyjs::html("console_binning", m$message, T)
                   })
                   
                 })
                 
                 binningFile <- reactive({
                   req(input$binningFile)
                   fread(input$binningFile$datapath)
                 })
                 
                 data <- reactive({
                   if(input$use_upload){
                     df <- binningFile()
                     
                   }else{
                     df <- binnedData()$data
                   }
                   
                   
                   list(data = df,
                        x_var=binnedData()$x_var)
                 })
                 
                 data
                 
               })
}

