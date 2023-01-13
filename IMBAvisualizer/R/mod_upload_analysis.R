
mod_upload_analysis_ui <- function(id) {
  ns <- NS(id)
  tagList(         
    fluidRow(
    column(2,
           textInput(
             ns('base_dir_upload'),
             label = "Base directory",
             value = "C:"),
           shinyFilesButton(
             ns("analyzed_dir"),
             "Choose analyzed data",
             "Choose file",
             multiple = F
           )
    ),



    ),
    wellPanel(
      h4("Metadata:"),
      verbatimTextOutput(ns("metadata_text"))
    ),
    wellPanel(
      h4("Data:"),
      DTOutput(ns('analyzed_table'))
    )
    )
  
  }

mod_upload_analysis_server <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      ns <- session$ns
      
      updateFileChooseA <- function(folder) {
        shinyFiles::shinyFileChoose(input = input,
                                    "analyzed_dir",
                                    roots = c(chosenFolder = folder))
      }
      
      observe({
        updateFileChooseA(input$base_dir_upload)
      })
      
      upload_rds <- reactiveValues(
          upload = list(timestamp = NULL,
          output_dir =NULL,
          save_as = NULL,
          radius =NULL,
          drop_contours = NULL,
          frame_rate = NULL,
          data = NULL,
          video_length= NULL,
          threshold_head_vector_angular_speed =NULL,
          threshold_tail_vector_angular_speed =NULL,
          min_hc_size = NULL,
          max_hc_size = NULL)
      )
      
      observeEvent(input$analyzed_dir,{
        print("trigger upload")
        if(length(input$analyzed_dir) <= 1){
          print("none selected... ")
          return({})
        }
        withCallingHandlers({
          message("read analyzed data...")
          directory = paste(input$analyzed_dir$files$`0`, collapse = "/")
          directory = paste0(input$base_dir_upload, directory)
          upload_rds$upload <- readRDS(directory)
          upload_rds$upload$data$group_condition <- paste(upload_rds$upload$data$group,
                                                          upload_rds$upload$data$condition,sep="-")
          message("finished reading!")
          #upload_rds$upload
        },
        message = function(m) {
          shinyjs::html("console", m$message, FALSE)
        }
        )
      })
      
      
      
      
    
      output$metadata_text <- renderText({
        paste(
          "Groups:",
          paste("Group:",unique(upload_rds$upload$data$group_condition),collapse="\n"),
          paste("Time analyzed:",upload_rds$upload$timestamp,sep=" "),
          paste("Folder name:",upload_rds$upload$save_as,sep=" "),
          paste("Output location:",upload_rds$upload$output_dir,sep=" "),
          paste("Radius:", upload_rds$upload$radius,sep=" "),
          paste("Frame rate:",upload_rds$upload$frame_rate ,sep=" "),
          paste("Video length:",upload_rds$upload$video_length ,sep=" "),
          paste("Drop contours:",upload_rds$upload$drop_contours,sep=" "),
          paste("Threshold head_vector_angular_speed:",
                upload_rds$upload$threshold_head_vector_angular_speed,sep=" "),
          paste("Threshold tail_vector_angular_speed:",
                upload_rds$upload$threshold_tail_vector_angular_speed,sep=" "),
          paste("Minimum head cast size:",upload_rds$upload$min_hc_size,sep=" "),
          paste("Maximum head cast size:",upload_rds$upload$max_hc_size,sep=" "),
          sep = "\n"
        )
        
      })
      
      output$analyzed_table <- renderDT({
        req(upload_rds$upload$data)
        message = "render analyzed table ..."
        shinyjs::html("consoleOutput", message)
        upload_rds$upload$data
      },
      options=list(scrollX=T)
      )
      
      upload <- reactive(upload_rds$upload)
      
      ##### RETURN #####
        return(
          upload
        )
    })
}