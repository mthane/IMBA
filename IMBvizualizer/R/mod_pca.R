
mod_pca_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidRow(
      column(4,wellPanel(
        numericInput(ns("ncomp"),"Number of components",value = 10,max = 100)
        
      )),
      column(8,
             
             wellPanel(
               plotlyOutput(ns("scree_plot")),
               plotlyOutput(ns("loading_plot")),
               linebreaks(30)
             )
      )
      
    )
  )
}

mod_pca_server <- function(id,data,variables) {
  moduleServer(
    id,
    function(input, output, session){
      df <- reactive({data()%>%
          select(c(variables(),"group_condition"))%>%
          na.omit()})
      
      M<- reactive({
        
        data <- df()%>%
          select(c(variables()))%>%
          na.omit()
        ### changing labels to variable names
        varNames <- global_vars()$SUMMARY_VARTABLE %>%
          filter(variable %in% colnames(data))
        varNames <- varNames[match(colnames(data), varNames$variable), ]
        colnames(data) <- varNames$name
        
        data
      })
      
      pca_res <- reactive(prcomp(M(), scale. = TRUE))
      
      output$loading_plot <- renderPlotly({
        print(pca_res())
        
        p<- autoplot(pca_res(), data = df(), colour = 'group_condition',
                     loadings = TRUE, loadings.colour = 'blue',
                     loadings.label = TRUE, loadings.label.size = 3, frame.type = 'norm')+
          theme_classic()
        ggplotly(p,height=900)
        
      })
      output$scree_plot <- renderPlotly({
        p <- fviz_eig(pca_res(),ncp = input$ncomp)
        ggplotly(p)
      })
      
    }
    
  )
}