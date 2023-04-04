

mod_2dbinning_ui <- function(id) {
  ns <- NS(id)
  
  
  fluidRow(
    column(4,
           wellPanel(
             selectInput(ns("heatmap_variable_z"), "Variable",
                         global_vars()$BINNING_VARIABLES),
             
             uiOutput(ns("group_selection")),
             sliderInput(ns("heatmap_nbins"),
                         
                         
                         "Bins",
                         min= 20,
                         max = 300,
                         value = 50),
             
             sliderInput(ns("heatmap_bandwidth"),
                         
                         
                         "Bandwidth",
                         min=1,
                         max=50,
                         value =5),
             uiOutput(ns("colorbarLimits"))
           )
           
           
           
    ),
    
    column(4,
           wellPanel(
             tabsetPanel(
               tabPanel("Dish",
                        actionButton(ns("actionHeatmapDish"),"Plot"),
                        plotlyOutput(ns("plot_heatmap"))
               ),
               
               tabPanel("Bearing-Distance",
                        actionButton(ns("actionHeatmap_bearing_distance"),"Plot"),
                        plotlyOutput(ns("plot_heatmap_bearing_distance"))
               ),                
               tabPanel("Bearing-Time",
                        actionButton(ns("actionHeatmap_bearing_time"),"Plot"),
                        plotlyOutput(ns("plot_heatmap_bearing_time"))
               ),
               tabPanel("Time-Distance",
                        actionButton(ns("actionHeatmap_time_distance"),"Plot"),
                        plotlyOutput(ns("plot_heatmap_time_distance"))
               )
             )
           )
           
           
    ))
  
  
}



mod_2dbinning_server <- function(id, upload) {
  moduleServer(id,
               function(input, output, session) {
                 
                 output$group_selection <- renderUI({
                   ns <- session$ns
                   selectInput(ns("group"),
                               "Experimental Condition",
                               choices=upload$groups())
                   
                 })
                 heatmap_variable <- reactive({
                   global_vars()$HEATMAP_VARTABLE%>%
                     filter(variable == input$heatmap_variable_z)
                 })
                 
                 output$colorbarLimits <- renderUI({
                   ns <- session$ns
                   print(input$heatmap_variable_z)
                   print(heatmap_variable())
                   sliderInput(
                     ns("heatmap_limits"),
                     heatmap_variable()$name,
                     heatmap_variable()$min,
                     heatmap_variable()$max,
                     c(heatmap_variable()$min, heatmap_variable()$max)
                   )
                   
                 })
                 
                 limits <- reactive(input$heatmap_limits)
                 
                 ##### Heatmap x and y #####
                 
                 heatmapDish_plot <- eventReactive(
                   input$actionHeatmapDish,{
                     create2dHeatmap(upload$frameData()%>%filter(group_condition==input$group),
                                     x = "spinepoint_x_6_conv",
                                     y = "spinepoint_y_6_conv",
                                     
                                     z = input$heatmap_variable_z,
                                     nbins = input$heatmap_nbins,
                                     filter_radius = input$heatmap_bandwidth,
                                     
                                     
                                     xlabel= "X [mm]",
                                     ylabel= "y [mm]",
                                     ratio_fixed = T,
                                     limits=limits(),
                                     radius=as.numeric(upload$radius())
                     )
                     
                   })
                 output$plot_heatmap <- renderPlotly(heatmapDish_plot())
                 
                 heatmap_bearing_distance_plot <- eventReactive(
                   input$actionHeatmap_bearing_distance,{
                     create2dHeatmap(upload$frameData()%>%filter(group_condition==input$group),
                                     x = "bearing_angle",
                                     y = "distance_to_odor",
                                     
                                     z = input$heatmap_variable_z,
                                     nbins = input$heatmap_nbins,
                                     filter_radius = input$heatmap_bandwidth,
                                     
                                     xlabel= "Bearing Angle [°]",
                                     ylabel= "Distance to Odor [mm]",
                                     ratio_fixed = F,
                                     limits=limits()
                     )
                   })
                 output$plot_heatmap_bearing_distance <- renderPlotly(heatmap_bearing_distance_plot())
                 
                 
                 heatmap_bearing_time_plot <- eventReactive(
                   input$actionHeatmap_bearing_time,{
                     create2dHeatmap(upload$frameData()%>%filter(group_condition==input$group),
                                     x = "bearing_angle",
                                     y = "time",
                                     z = input$heatmap_variable_z,
                                     nbins = input$heatmap_nbins,
                                     filter_radius = input$heatmap_bandwidth,
                                     xlabel= "Bearing Angle [°]",
                                     ylabel= "Time [s]",
                                     ratio_fixed = F,
                                     limits=limits()
                     )
                   })
                 output$plot_heatmap_bearing_time <- renderPlotly(heatmap_bearing_time_plot())
                 
                 heatmap_time_distance_plot <- eventReactive(
                   input$actionHeatmap_time_distance,{
                     create2dHeatmap(upload$frameData()%>%filter(group_condition==input$group),
                                     x = "time",
                                     y = "distance_to_odor",
                                     z = input$heatmap_variable_z,
                                     nbins = input$heatmap_nbins,
                                     filter_radius = input$heatmap_bandwidth,
                                     xlabel= "Time [s]",
                                     ylabel= "Distance [mm]",
                                     ratio_fixed = F,
                                     limits=limits()
                     )
                   })
                 output$plot_heatmap_time_distance <- renderPlotly(heatmap_time_distance_plot())
                 
                 
                 
               }
  )
}