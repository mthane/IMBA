
mod_parcoord_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidRow(
      column(4,wellPanel(
        
        sliderInput(ns("spline"),"Spline Factor",1.,10.,1),
        radioButtons(ns("order"),"Ordering",c("allClass","skewness")),
        checkboxInput(ns("costum_order"),"Use costum order",F),
        uiOutput(ns("ranked_list"))
      )),
      column(8,
             
             wellPanel(
               plotlyOutput(ns("parcoord")),
               linebreaks(30)
             )
      )
      
    )
  )
}

mod_parcoord_server <- function(id,data,variables) {
  moduleServer(
    id,
    function(input, output, session){
      df <- reactive({data()%>%
          select(c(variables(),"group_condition"))})
      
      # M<- reactive({
      #   print(df())
      #   data <- df()%>%
      #     select(c(variables()))%>%
      #     na.omit()
      #   ### changing labels to variable names
      #   varNames <- summaryNames %>%
      #     filter(variable %in% colnames(data))
      #   varNames <- varNames[match(colnames(data), varNames$variable), ]
      #   colnames(data) <- varNames$name
      #   
      #   data
      # })
      
      
      output$ranked_list <- renderUI({
        ns <- session$ns
        rank_list(
          text = "Drag the items in any desired order",
          labels = variables(),
          input_id = ns("rank_vars")
        )
      })
      
      output$parcoord <- renderPlotly({
        df = df()
        if(input$costum_order){
          
          ordering = rev(match(input$rank_vars, colnames(df)))
          
        }else{
          ordering = input$order
        }
        p<- ggparcoord(df,columns = seq(1,ncol(df)-1),
                       groupColumn = "group_condition",
                       order=ordering,
                       splineFactor=input$spline,
                       
                       missing = "mean",
                       alphaLines = 0.2)+ 
          theme_classic()+
          coord_flip()
        #theme(axis.text.x = element_text(angle = 90))
        ggplotly(p,height=900)
        
      })
      
    }
    
  )
}