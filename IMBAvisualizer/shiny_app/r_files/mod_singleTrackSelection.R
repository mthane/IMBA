
mod_singleTrackSelection_ui <- function(id, label = "Track selection"){
  ns <- NS(id)
  tagList(
    wellPanel(
      h4("Track Selection"),
      fluidRow(
               uiOutput(ns("empty_data")),
               DTOutput(ns('tracks_tbl'))
        )
        
      )
    )
}

mod_singleTrackSelection_server <- function(id,upload,threshold) {
  moduleServer(
    id,
    function(input, output, session) {
      ns <- session$ns
      ##### SINGLE MODE: select track ####
      output$empty_data <- renderUI({
        ns <- session$ns
        span(textOutput(ns('empty_data_text')),style="color:red")

      })
      observe(print(upload()))
      output$empty_data_text <- renderText({

        if(is.null(upload()$data)){
        "No data has been uploaded! Go to 'Upload analysis' and select a valid .Rds file!"
        }
        else{
          ""
        }
      })

      frameDataIndividual <- reactive({
        print(upload()$data)
        req(upload()$data)
          upload()$data%>%
            group_by(id) %>%
            filter(n() > upload()$frame_rate*threshold())
      })

      tracks_table <- reactive({
        req(frameDataIndividual())
        frameDataIndividual() %>%
            group_by(id, name, trial, condition, group) %>%
            summarise()
      })


      selected_tracks <- reactive({
        req(tracks_table())
        df <- tracks_table()
        df$id[input$tracks_tbl_rows_selected]%>%unlist()
      })

      track <- reactive({
        req(frameDataIndividual())
        frameDataIndividual()%>%
          filter(id%in%selected_tracks())
      })


      output$tracks_tbl <- DT::renderDT({
        req(tracks_table())
        tracks_table()
      },options=list(scrollX=T),
      selection = list(mode = 'multiple',
                       selected = 1))


      
      
      return(track)

    }
  )}
