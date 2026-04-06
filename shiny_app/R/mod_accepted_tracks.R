
mod_accepted_tracks_ui <- function(id) {
  ns <- NS(id)
  tagList(
    wellPanel(
      fluidRow(numericInput(ns("threshold"), label = "Threshold",value = 90,min = 0, max = 240)),
      fluidRow(
        column(4,
               plotOutput(ns("plot_track_lengths")),
               downloadButton(ns('download_track_lengths'))    
              ),
        column(4,
               
               plotOutput(ns("plot_accepted_tracks")),
               downloadButton(ns('download_accepted_tracks'))
        ),
        column(4,
               plotOutput(ns("plot_non_accepted_tracks")),
               downloadButton(ns('download_non_accepted_tracks'))
        )
      )
    )

    
  )
}
  
mod_accepted_tracks_server <- function(id,upload) {
  moduleServer(
    id,
    function(input, output, session) {
      
      threshold <- reactive(input$threshold) 
      
      
      plotTrackLengths = reactive({
        threshold <- input$threshold
      data <- upload()$data%>%
        as.data.frame()%>%
        select(id)%>%
        group_by(id)%>%
        summarise(n= n()/upload()$frame_rate)%>%
        mutate(acc = as.factor(n>=threshold))
      
      ggplot(data)+
        geom_histogram(aes(n,fill=acc),bins=50,alpha=0.5)+
        scale_fill_manual(
          name = "Condition",
          labels = c("Non-accepted","Accepted"),
          values = c("#bababa","#2cf28f")
        )+
        geom_vline(xintercept=threshold, linetype="dashed",size=1) +
        labs(
          title = "Distribution of track lengths",
          x = "Track length",
          y = "N"
        )+
        theme_classic()
      })
      
      output$plot_track_lengths <- renderPlot({
        plotTrackLengths()
      })
      
      output$download_track_lengths = downloadHandler(
        filename = 'plot_track_lengths.png',
        content = function(file) {
          ggplot2::ggsave(file, plot = plotTrackLengths())
        })
      
      
      plotAcceptedTracks <- reactive({
        radius <- as.numeric(upload()$radius)
        
        data_acc <-upload()$data%>%
          as.data.frame()%>%
          group_by(id)%>%
          filter(n() > upload()$frame_rate*input$threshold)%>%
          ungroup()%>%
          #sampling every 0.5 sec
          filter(row_number() %% 4 == 1)
        
        perc_acc <- round(length(unique(data_acc$id))/length(unique(upload()$data$id)),4)*100
        perc_acc_frames <- round(nrow(data_acc)*4/nrow(upload()$data)*100,2)
        
        ggplot(data_acc)+
          geom_point(aes(spinepoint_x_6_conv,spinepoint_y_6_conv,group=id),size=0.5,alpha=0.1,color="#2cf28f")+
          xlim(-radius, radius) +
          ylim(-radius, radius) +
          annotate("path",
                   x = 0 + radius * cos(seq(0, 2 * pi, length.out = 100)),
                   y = 0 + radius * sin(seq(0, 2 * pi, length.out = 100)))+
          theme_classic()+
          coord_fixed()+
          labs(
            title = paste("Accepted tracks ","(",perc_acc,"%) ","Accepted frames ","(",perc_acc_frames,"%)",sep=""),
            caption = paste0("Tracks with a length > ",input$threshold),
            x = "x [mm]",
            y = "y [mm]"
          )
      })
      
      output$plot_accepted_tracks <- renderPlot({
        plotAcceptedTracks()
      })
      
      output$download_accepted_tracks = downloadHandler(
        filename = 'plot_accepted_tracks.png',
        content = function(file) {
          ggplot2::ggsave(file, plot = plotAcceptedTracks())
        })
      
      
      
      plotNonAcceptedTracks <- reactive({
        radius <- as.numeric(upload()$radius)
        data_nacc <-   upload()$data%>%
          as.data.frame()%>%
          group_by(id)%>%
          filter(n() <= upload()$frame_rate*input$threshold)%>%
          ungroup()%>%
          #sampling every 0.5 sec
          filter(row_number() %% 4 == 1)
        
        perc_nacc <- round(length(unique(data_nacc$id))/length(unique(upload()$data$id)),4)*100
        perc_nacc_frames <- round(nrow(data_nacc)*4/nrow(upload()$data)*100,2)
        
        ggplot(data_nacc)+
          geom_point(aes(spinepoint_x_6_conv,spinepoint_y_6_conv,group=id),size=0.5,alpha=0.1,color= "#bababa")+        
          xlim(-radius, radius) +
          ylim(-radius, radius) +
          annotate("path",
                   x = 0 + radius * cos(seq(0, 2 * pi, length.out = 100)),
                   y = 0 + radius * sin(seq(0, 2 * pi, length.out = 100)))+
          theme_classic()+
          coord_fixed()+
          labs(
            title = paste("Non-accepted tracks ","(",perc_nacc,"%) ","Non-accepted frames ","(",perc_nacc_frames,"%)",sep=""),
            caption = paste0("Tracks with a length <= ",input$threshold),
            x = "x [mm]",
            y = "y [mm]"
          )
        
      })
      
      output$plot_non_accepted_tracks <- renderPlot({
        plotNonAcceptedTracks()
      })
      
      output$download_non_accepted_tracks = downloadHandler(
        filename = 'plot_non_accepted_tracks.png',
        content = function(file) {
          ggplot2::ggsave(file, plot = plotNonAcceptedTracks())
        })
      
      
      threshold
    })
    }