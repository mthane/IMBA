
mod_uTests_ui <- function(id){
  ns <- NS(id)
  tagList(
      wellPanel(
           uiOutput(ns("groupSelection")),
           plotlyOutput(ns("uTests"))

      )

    
  )
  
}

mod_uTests_server <- function(id,data,variables) {
  moduleServer(
    id,
    function(input, output, session) {
      
      output$groupSelection <- renderUI({
        ns <- session$ns
        groups <- unique(data()$group_condition)
        tagList(
          selectInput(ns("group1"),label = "Group 1", choices = groups),
          selectInput(ns("group2"),label = "Group 2", choices = groups,selected=groups[2])
          
        )
      })
      
      output$uTests <- renderPlotly({
        # helper function needed for tranformation
        reverselog_trans <- function(base = exp(1)) {
          trans <- function(x) -log(x, base)
          inv <- function(x) base^(-x)
          trans_new(paste0("reverselog", format(base)), trans, inv, 
                    log_breaks(base = base), 
                    domain = c(1e-100, Inf))
        }
        
        df <- data()%>% 
          filter(group_condition %in% c(input$group1,input$group2))%>%
          uTests(variables())
        #get real names
        varNames <- global_vars()$SUMMARY_VARTABLE %>%
          filter(variable %in% variables())
        varNames <- varNames[match(df$variable, varNames$variable), ]
        print(varNames)
        df$variable <- unlist(varNames$name)
        
        p <- df%>%
          ggplot(aes(x=pvalue,y=reorder(variable,-pvalue),fill=pvalue))+
          geom_bar(stat="identity", position=position_dodge())+
          geom_vline(xintercept=0.05,color="red")+
          theme_classic()+
          scale_fill_viridis(direction=-1)+
          ylab("Variable")+
          xlab("P-Value")+
          theme(axis.title.y = element_text(element_text(vjust = -0.15)))
        #scale_x_continuous(trans=reverselog_trans(10))
        ggplotly(p)%>%
          config(
            toImageButtonOptions = list(
              format = "svg",
              filename = paste0("uTests"),
              width = 600,
              height = 700
            )
          )
      })
      
    }
  )
}


