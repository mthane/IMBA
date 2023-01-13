
mod_randomForest_ui <- function(id){
  ns <- NS(id)
  tagList(
    
    fluidRow(
      
      column(4,
             wellPanel(
               numericInput(ns("kfolds"), "K (fold)", 10),
               numericInput(ns("repeated"), "Repitition", 3),
               sliderInput(ns("ntree"),"Number of Trees",100,5000,500),
               radioButtons(ns("sampling"),label = "Sampling",choices = c("up","down")),
               actionButton(ns("train_model"),"Train model")
             )
             
      ),
      column(8,
             wellPanel(
               fluidRow(
                 column(6,               
                        plotlyOutput(ns("plot_cm"))
                 ),
                 column(6,
                        verbatimTextOutput(ns("printSummary"))
                 )
               ),
               plotlyOutput(ns("plot_vimp")),
               numericInput(ns("showN_vimp"),"Number of Variables",5,60,value = 15)
               
             )
      )
    )
    
  )
  
}
mod_randomForest_server <- function(id, data,variables) {
  moduleServer(id,
               function(input, output, session) {
                 ns <- session$ns
                 # training the model by assigning sales column 
                 # as target variable and rest other column 
                 # as independent varaible 
                 model <- eventReactive(input$train_model,{
                   n <- 2
                   i=1
                   withProgress(message = 'Making model' ,{
                     
                     incProgress(1/n, detail = "Preprocessing data")
                     target = "group_condition"
                     my_data = data()%>%select(all_of(variables()))
                     my_data$group_condition <- as.factor(data()$group_condition)
                     df <-my_data
                     df[is.na(df)] = 0
                     form <- as.formula(paste(target,paste(variables(),collapse=" + "),sep=" ~ "))
                     
                     if(input$sampling=="down"){
                       sampled <-  downSample(x = df,
                                              y = df$group_condition)%>%
                         na.omit()
                     }
                     if(input$sampling=="up"){
                       sampled <-  upSample(x = df,
                                            y = df$group_condition)%>%
                         na.omit()
                     }
                     
                     incProgress(1/n, detail = "Train Model...")
                     train_control <- trainControl(method = "repeatedcv",  
                                                   number = input$kfolds,
                                                   repeats = input$repeated,
                                                   savePredictions = T) 
                     
                     showNotification("Training model...")
                     model <- train(form, 
                                    data = sampled,  
                                    fitBest=TRUE,
                                    returnData = TRUE,
                                    method = "rf",
                                    ntrees = input$ntree,
                                    trControl = train_control)
                     showNotification("Training model finished.")
                   })
                   model
                 })
                 
                 
                 
                 output$plot_vimp <- renderPlotly({
                   model <- model()
                   labels <- variables()
                   
                   
                   vimp <- data.frame(imp = varImp(model)$importance,
                                      label = rownames(varImp(model)$importance))%>%
                     arrange(-Overall)%>%
                     slice(1:input$showN_vimp)
                   
                   #get real names
                   varNames <- global_vars()$SUMMARY_VARTABLE %>%
                     filter(variable %in% vimp$label)
                   varNames <- varNames[match(vimp$label, varNames$variable), ]
                   vimp$label <- unlist(varNames$name)
                   
                   colnames(vimp) = c("Importance", "Variable")
                   
                   print(vimp)
                   vimp.plot <- ggplot(vimp,aes(x=Importance,y=reorder(Variable,Importance)))+
                     #geom_point()+
                     geom_bar(stat = "identity",aes(fill=Importance))+
                     scale_fill_viridis_c()+
                     ylab("Variable")+
                     theme_classic()
                   
                   
                   p<-ggplotly(vimp.plot)%>%
                     config(
                       toImageButtonOptions = list(
                         format = "svg",
                         filename = paste0("variableImportance"),
                         width = 600,
                         height = 700
                       )
                     )
                   
                   p
                 })
                 
                 
                 output$printSummary <- renderPrint({
                   print(model())
                 })
                 # 
                 output$plot_cm <- renderPlotly({
                   
                   obs <- extractPrediction(list(model()))$obs
                   pred <- extractPrediction(list(model()))$pred
                   
                   table <- data.frame(confusionMatrix(model())$table)
                   
                   plotTable <- table %>%
                     mutate(goodbad = ifelse(table$Prediction == table$Reference, "good", "bad")) %>%
                     group_by(Reference) %>%
                     mutate(prop = Freq/sum(Freq))
                   
                   # fill alpha relative to sensitivity/specificity by proportional outcomes within reference groups (see dplyr code above as well as original confusion matrix for comparison)
                   p<-ggplot(data = plotTable, mapping = aes(x = Reference, y = Prediction, fill = prop)) +
                     geom_tile() +
                     geom_text(aes(label = round(prop,2)), vjust = .5, fontface  = "bold", alpha = 1) +
                     scale_fill_viridis_c(limits=c(0,1))+
                     
                     xlim(rev(levels(table$Reference)))+
                     theme_classic()+
                     scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+
                     scale_y_discrete(labels = function(x) str_wrap(x, width = 10))+
                     coord_fixed()+ 
                     theme(axis.text.x = element_text(angle = 90))
                   
                   ggplotly(p)%>%
                     config(
                       toImageButtonOptions = list(
                         format = "svg",
                         filename = paste0("confusionMatrix"),
                         width = 600,
                         height = 700
                       )
                     )
                   
                 })
                 
               })
}
