
mod_scatterMulti_ui <- function(id) {
  ns <- NS(id)
  

      wellPanel(
        selectInput(
          inputId = ns("scatter_variable"),
          label = "Variables: ",
          choices = global_vars()$SUMMARY_VARIABLES,
          multiple=T
        ),
        plotlyOutput(ns("plot_scatter"),
                     height="800px",
                     width="800px")

  )
  
  
}


mod_scatterMulti_server <- function(id,data,colorValues) {
  moduleServer(
    id,
    function(input, output, session) {
      
      output$plot_scatter <- renderPlotly({
        
        req(data())
        row = global_vars()$SUMMARY_VARTABLE%>%
          filter(variable%in%input$scatter_variable)
        
        names = row$name
        data <- data() %>%
          select(c(row$variable, "group_condition"))%>%
          as.data.frame()
        N = ncol(data)-1
        #print(count(data,group_condition))
        if(N>1 ){
          N = ncol(data)-1
          data$group_condition <- factor(data$group_condition)
          p <- ggpairs(data, mapping = aes(color=stringr::str_wrap(group_condition,10),
                                           alpha = 0.7),
                       columns = 1:N,
                       columnLabels =names
          )+
            theme_classic()
          
          
          #changing the colors
          if(!is.null(unlist(colorValues()$groupColorValues))){
            for(i in 1:p$nrow) {
              for(j in 1:p$ncol){
                p[i,j] <- p[i,j] + 
                  scale_fill_manual(values=colorValues()$groupColorValues) +
                  scale_color_manual(values=colorValues()$groupColorValues)
              }
            }
          }
          ggplotly(p)%>% 
            highlight("plotly_selected")%>%
            config(
              toImageButtonOptions = list(
                format = "svg",
                filename = paste0("Boxplot.svg"),
                width = 600,
                height = 700
              )
            )
          
        }
        
        
      })
    })
  
}




