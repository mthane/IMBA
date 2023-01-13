

mod_correlationAnalysis_ui <- function(id){
  ns <- NS(id)
               wellPanel(
                 radioButtons(ns("correlationType"),"Correlation type", c("spearman", "pearson")),
                 checkboxInput(ns("absolute"),"Absolute Correlation",value = FALSE),

                 tabsetPanel(
                   tabPanel("Correlation Matrix",
                            plotlyOutput(ns("cmatrix")),
                            
                            linebreaks(30)
                            
                   ),
                   tabPanel("Dendogram",
                            
                            radioButtons(ns("hc_method"),"Clustering Type", c("ward.D",
                                                                              "ward.D2", 
                                                                              "average", 
                                                                              "single",
                                                                              "complete")),
                            radioButtons(ns("dend_type"),"Dendogram Type", c("rectangle", "circular", "phylogenic")),
                            plotOutput(ns("correlogram")),
                            linebreaks(30),
                            sliderInput(ns("K"),"K",2,30,5),
                            plotlyOutput(ns("sil")),
                            linebreaks(40),
                            plotlyOutput(ns("clusterVis")),
                            linebreaks(30)
                            
                   ),
                   tabPanel("Network",
                            sliderInput(ns("minCor"),
                                        label = "Minimum Correlation",min = 0,max=1,value = 0.3),
                            
                            plotOutput(ns("corrnetwork")#,
                            #linebreaks(30)
                            )
                            
                   )
                   
                 )
               )
  
}

mod_correlationAnalysis_server <- function(id,data,variables) {
  moduleServer(
    id,
    function(input, output, session){
      labels <- reactive({
        global_vars()$SUMMARY_VARTABLE%>%
          filter(variable %in% variables())%>%
          select(name)%>%
          unlist()
        
      })
      M<- reactive({
        data <-data()%>%
          select(variables())%>%
          na.omit()%>%
          as.matrix()
        ### changing labels to variable names
        varNames <- global_vars()$SUMMARY_VARTABLE %>%
          filter(variable %in% colnames(data))
        varNames <- varNames[match(colnames(data), varNames$variable), ]
        colnames(data) <- varNames$name
        
        data
      })
      
      output$cmatrix <- renderPlotly({
        
        pM <- cor.mtest(M(), conf.level = .95)
        P <- pM$p
        cm <- cor(M(),method=input$correlationType )
        if(input$absolute){
          cm <- abs(cm)
        }
        p<-ggcorrplot(cm,
                      hc.order = TRUE, outline.color = "white")+
          theme(axis.text.x = element_text(size=7))+
          theme(axis.text.y = element_text(size=7))
        
        ggplotly(p)%>%
          layout(height=1000,width=1000)
      })
      cm <- reactive({
        cm <- cor(M(),method=input$correlationType)
        if(input$absolute){
          cm <- abs(cm)
        }
        cm
      })
      res.hc <- reactive({
        
        eclust(cm(),
               "hclust",
               k = input$K,
               hc_method = input$hc_method,
               graph = FALSE)
        
      })
      
      
      output$correlogram <- renderPlot({
        
        p<-fviz_dend(res.hc() ,type = input$dend_type,repel=T)
        if(input$dend_type=="rectangle"){
          
          p <- p+ theme(axis.text.x = element_text(angle = 90))
        }
        p
      },height=900,width=900)
      
      output$sil <- renderPlotly({
        p <- fviz_silhouette(res.hc())+
          theme_classic()+
          coord_flip()
        #theme(axis.text.x = element_text(angle = 90))
        ggplotly(p,height=1200)
      })
      
      output$clusterVis <- renderPlotly({
        p <- fviz_cluster(res.hc(),data=M(),ggtheme=theme_classic())
        ggplotly(p,height=900)
      })
      
      
      output$corrnetwork <- renderPlot({
        M() %>% correlate() %>% 
          network_plot(min_cor = input$minCor)
      },height=900,width=900)
      
    }
  )
  
}