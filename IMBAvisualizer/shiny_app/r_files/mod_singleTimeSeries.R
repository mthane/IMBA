#' name_of_module1 UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom plotly plotlyOutput
#'
mod_singleTimeSeries_ui <- function(id, label = "Track selection") {
  ns <- NS(id)
  tagList(wellPanel(
    h4("Time Series Visualization"),
    tagList(
      
      tabsetPanel(
        tabPanel(
                 "Trajectories",
                 fluidRow(
                   pickerInput(
                     inputId = ns("TS_variable"),
                     label = "Variable",
                     choices = global_vars()$SINGLE_VARIABLES,
                     options = list(`live-search` = TRUE),
                     selected = c("tail_vel_forward"),
                     multiple = F
                   ),
                   sliderInput(
                     ns('bins_single'),
                     'Number of bins',
                     value = 10,
                     min = 10,
                     max = 100
                   ),
                   checkboxInput(ns("density"),"Use Density Plot"),
                   
                 ),
                 fluidRow(
                   column(6,
                          plotOutput(ns("track_on_dish_plot")),
                          
                          downloadButton(ns("download_track_on_dish"))
                          ),
                   column(6,
                          plotOutput(ns(
                            "plot_trackOnDish_colorCode_TS"
                          )),
                          
                          downloadButton(ns("download_track_on_dish_colorCode_TS"))
                          )
                          
                   ),
                        
                 fluidRow(
                   column(6,
                          plotOutput(ns(
                            "plot_trackOnDish_colorCode_ID"
                          )),
                          
                          downloadButton(ns("download_track_on_dish_colorCode_ID"))
                          ),
                   
                   
                   column(6,
                          
                          plotlyOutput(ns(
                            "plot_single_histogram_TS"
                          ))
                          )
                 )
                 

               ),
        tabPanel(
          "Time Series Plot",
               fluidRow(
                 column(
                   
                   6,
                   pickerInput(
                     inputId = ns("TS_variables"),
                     label = "Variable",
                     choices = global_vars()$SINGLE_VARIABLES,
                     options = list(`live-search` = TRUE),
                     selected = c("tail_vel_forward"),
                     multiple = T
                   ),
                   checkboxInput(ns("scale"), label = "Scale time series", value = F),
                   checkboxInput(ns("use_limits"), label = "Use axis limits", value = F)
                 ),
                 column(
                   6,
                   sliderInput(
                     ns("slider_time"),
                     "Time:",-60,
                     900,
                     c(0, 180),
                     step = 1,
                     width = 1500
                   ),
                   sliderInput(ns("ylimits"), "Y-axis Limits", -10, 10, value = c(-5, 5),step = 0.1),
                 )
               ),
               
               fluidRow(
                 uiOutput(ns("plot_timeseries"))
               )
               )
        
      )


      
    )
  ))
}


#' name_of_module1 UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom plotly plotlyOutput renderPlotly
#'
mod_singleTimeSeries_server <- function(id, track, upload) {
  moduleServer(id,
               function(input, output, session) {
                 output$plot_TS_variables_timeseries <- renderPlotly({
                   req(track())
                   plot_timeseries(
                     track(),
                     input$TS_variables,
                     global_vars()$SINGLE_VARTABLE,
                     input$slider_time,
                     input$scale,
                     input$use_limits,
                     input$ylimits,
                     upload()$frame_rate
                   )
                 })
                 
                 output$plot_timeseries <- renderUI({
                   ns <- session$ns
                   tagList(plotlyOutput(ns('plot_TS_variables_timeseries')),
                           linebreaks(15 * length(unique(track(
                           )$id))))
                   
                 })
                 
                 track_on_dish <- reactive({
                   req(track())
                   data <- track()
                   odor_x_r = unlist(unique(data$odor_a_x_rot))
                   odor_y_r = unlist(unique(data$odor_a_y_rot))
                   radius = as.numeric(upload()$radius)
                   data$HC_LR = as.factor(data$HCs_right -data$HCs_left)
                   track_end = tail(data %>%select(spinepoint_x_6_conv,
                                                   spinepoint_y_6_conv,
                   )%>%na.omit(), 1)
                   name = unique(data$name)[1]
                   id = unique(data$id)[1]
                   
                   plot <- ggplot(data) +
                     xlim(-radius, radius) +
                     ylim(-radius, radius) +
                     annotate(
                       "path",
                       color = 'red',
                       alpha = 0.7,
                       x = odor_x_r + 3 * cos(seq(0, 2 * pi, length.out = 100)),
                       y = odor_y_r + 3 * sin(seq(0, 2 * pi, length.out = 100))
                     ) +
                     annotate("path",
                              x = 0 + radius * cos(seq(0, 2 * pi, length.out = 100)),
                              y = 0 + radius * sin(seq(0, 2 * pi, length.out = 100))) +
                     geom_point(
                       aes(x = spinepoint_x_6_conv,
                           y = spinepoint_y_6_conv,
                           
                           text = frame/16
                       ),
                       alpha = 0.1,
                       #size=0.6
                     )+
                     geom_point(
                       aes(x = spinepoint_x_6_conv,
                           y = spinepoint_y_6_conv,
                           colour= HC_LR
                       ),
                       alpha = 0.3,
                       #size=1,
                       data=data%>%filter(HC_LR!=0)
                     )+
                     scale_colour_manual(values =c("green",
                                                   "red"),
                                         labels =c("left","right"),
                                         breaks = c("-1", "1")
                     )+
                     labs(title=paste("Track:", name, "ID:",id))+
                     geom_point(
                       data = track_end,
                       mapping = aes(x = spinepoint_x_6_conv,
                                     y = spinepoint_y_6_conv),
                       #size =2,
                       color = "blue", alpha = 0.5)+
                     theme_classic()+
                     coord_fixed()
                   plot
                   
                 })
                 output$download_track_on_dish = downloadHandler(
                   filename = 'track_on_dish.png',
                   content = function(file) {
                     device <- function(..., width, height) {
                       grDevices::png(..., width = width, height = height,
                                      res = 300, units = "in")
                     }
                     ggsave(file, plot = track_on_dish(), device = device)
                   })
                 
                 output$track_on_dish_plot <- renderPlot({
                    track_on_dish()
                 })
                 
                 output$plot_single_histogram_TS <- renderPlotly({
                   req(track())
                   plot_single_histogram(track(),
                                         input$TS_variable,
                                         input$slider_time,
                                         input$bins_single,
                                         use_density = input$density
                                         )
                 })
                 
                 track_on_dish_colorCode_TS <- reactive({
                   req(track())
                   data <- track() %>%
                     filter(
                       frame >=  input$slider_time[1] * upload()$frame_rate &
                         frame <=  input$slider_time[2] * upload()$frame_rate
                     )
                   odor_x_r = unlist(unique(data$odor_a_x_rot))
                   odor_y_r = unlist(unique(data$odor_a_y_rot))
                   radius = as.numeric(upload()$radius)
                   track_end = tail(data %>% select(spinepoint_x_6_conv,
                                                    spinepoint_y_6_conv,) %>% na.omit(),
                                    1)
                   name = unique(data$name)[1]
                   id = unique(data$id)[1]
                   
                   plot <- ggplot(data) +
                     xlim(-radius, radius) +
                     ylim(-radius, radius) +
                     annotate(
                       "path",
                       color = 'red',
                       alpha = 0.7,
                       x = odor_x_r + 3 * cos(seq(0, 2 * pi, length.out = 100)),
                       y = odor_y_r + 3 * sin(seq(0, 2 * pi, length.out = 100))
                     ) +
                     annotate("path",
                              x = 0 + radius * cos(seq(0, 2 * pi, length.out = 100)),
                              y = 0 + radius * sin(seq(0, 2 * pi, length.out = 100))) +
                     geom_point(
                       aes_string(
                         x = "spinepoint_x_6_conv",
                         y = "spinepoint_y_6_conv",
                         color = input$TS_variable
                       ),
                       alpha = 0.1#,
                       #size=0.6
                     ) +
                     labs(title = paste("Track:", name, "ID:", id)) +
                     geom_point(
                       data = track_end,
                       mapping = aes(x = spinepoint_x_6_conv,
                                     y = spinepoint_y_6_conv),
                       #size =2,
                       color = "blue",
                       alpha = 0.5
                     ) +
                     scale_color_viridis() +
                     theme_classic() +
                     coord_fixed()
                   plot
                 })
                 
                 output$download_track_on_dish_colorCode_TS = downloadHandler(
                   filename = 'track_on_dish_color_TS.png',
                   content = function(file) {
                     device <- function(..., width, height) {
                       grDevices::png(..., width = width, height = height,
                                      res = 300, units = "in")
                     }
                     ggsave(file, plot = track_on_dish_colorCode_TS(), device = device)
                   })
                 
                 output$plot_trackOnDish_colorCode_TS <- renderPlot({
                    track_on_dish_colorCode_TS()
                 })
                 
                 track_on_dish_colorCode_ID <- reactive({
                   req(track())
                   data <- track() %>%
                     filter(
                       frame >=  input$slider_time[1] * upload()$frame_rate &
                         frame <=  input$slider_time[2] * upload()$frame_rate
                     )%>%
                     mutate(id = factor(id))
                   odor_x_r = unlist(unique(data$odor_a_x_rot))
                   odor_y_r = unlist(unique(data$odor_a_y_rot))
                   radius = as.numeric(upload()$radius)
                   track_end = tail(data %>% select(spinepoint_x_6_conv,
                                                    spinepoint_y_6_conv,) %>% na.omit(),
                                    1)
                   name = unique(data$name)[1]
                   id = unique(data$id)[1]
                   
                   plot <- ggplot(data) +
                     xlim(-radius, radius) +
                     ylim(-radius, radius) +
                     annotate(
                       "path",
                       color = 'red',
                       alpha = 0.7,
                       x = odor_x_r + 3 * cos(seq(0, 2 * pi, length.out = 100)),
                       y = odor_y_r + 3 * sin(seq(0, 2 * pi, length.out = 100))
                     ) +
                     annotate("path",
                              x = 0 + radius * cos(seq(0, 2 * pi, length.out = 100)),
                              y = 0 + radius * sin(seq(0, 2 * pi, length.out = 100))) +
                     geom_point(
                       aes_string(
                         x = "spinepoint_x_6_conv",
                         y = "spinepoint_y_6_conv",
                         color = "id"
                       ),
                       alpha = 0.5#,
                       #size=0.6
                     ) +
                     labs(title = paste("Track:", name, "ID:", id)) +
                     geom_point(
                       data = track_end,
                       mapping = aes(x = spinepoint_x_6_conv,
                                     y = spinepoint_y_6_conv),
                       #size =2,
                       color = "blue",
                       alpha = 0.5
                     ) +
                     theme_classic() +
                     coord_fixed()
                   plot
                 })
                 output$plot_trackOnDish_colorCode_ID <- renderPlot({
                   track_on_dish_colorCode_ID()
                 })
                 
                 output$download_track_on_dish_colorCode_ID = downloadHandler(
                   filename = 'track_on_dish_color_ID.png',
                   content = function(file) {
                     device <- function(..., width, height) {
                       grDevices::png(..., width = width, height = height,
                                      res = 300, units = "in")
                     }
                     ggsave(file, plot = track_on_dish_colorCode_ID(), device = device)
                   })
               })
}
