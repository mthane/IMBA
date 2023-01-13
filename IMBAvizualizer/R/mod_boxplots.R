

mod_boxplots_ui <- function(id) {
  ns <- NS(id)
  
    wellPanel(
      fluidRow(
        
        column(4,
               selectInput(
                 ns("boxplotsVariable"),
                 "Variable",
                 global_vars()$SUMMARY_VARIABLES
               )
               ),
        column(4,
               tags$head(tags$style(
                 HTML(
                   "
                  .form-control {
                      border-radius: 4px 4px 4px 4px;
                  }

                  #boxplotDescription {
                  font-family:  'Source Sans Pro','Helvetica Neue',Helvetica,Arial,sans-serif;
                  font-size: 14px;
                  width: 300px;
                  height: 500px;
                  max-width: 100%;
                  padding: 6px 12px;
                  white-space: pre-wrap;
                  }

                  "
                 )
               )),
               textOutput(ns("boxplotDescription"))
               ),
        column(4,
               
               checkboxInput(ns("violin"), "Violin plots", FALSE),
               
               checkboxInput(ns("rmoutliers"), "Remove Outliers", F)
               )
      ),


      fluidRow(
        plotlyOutput(
          ns("boxplotsSummary")
        ),
        
        DTOutput(ns(
          "sign_table"
        ))
      )
      

    )
  
  
  
}


mod_boxplots_server <- function(id,
                           data,
                           colorValues) {
  moduleServer(id,
               function(input, output, session) {
                 
                 
                 boxplot_variable <- reactive({
                   global_vars()$SUMMARY_VARTABLE%>%
                     filter(variable == input$boxplotsVariable)
                   
                 })
                 

                 output$boxplotDescription <- renderPrint({
                   cat(boxplot_variable()$description)
                 })
                 
                 boxplot <- reactive({
                   
                   #if(nrow(data()>0)){
                     boxplot_per_group(
                       data(),
                       boxplot_variable(),
                       colorValues()$groupColorValues,
                       colorValues()$groupNameValues,
                       remove_outliers = input$rmoutliers,
                       violin = input$violin,
                       height = input$height,
                       width = input$width
                     )
                   #   
                   # }else{
                   #   empty_plot("Create aggregated data using Upload -> Aggregated data -> Aggregate data")
                   # }
                 })
                 
                 output$boxplotsSummary <- renderPlotly({
                   
                     boxplot()
                 })
                 
                 sign_data <- reactive({
                   #req(data())
                   data = data()
                   
                   print(length(unique(data$group_condition)))
                   if(nrow(data)>0 & length(unique(data$group_condition))>1){
                     pvalues(data,
                             input$boxplotsVariable)
                   }else{
                     NULL
                   }

                   
                 })
                 
                 output$sign_table <- DT::renderDataTable({
                   req(sign_data())
                   DT::datatable(sign_data(),
                                 options = list(pageLength = 10,
                                                order = list(list(3, 'asc')))) %>%
                     formatStyle('pvalues',
                                 target = 'row',
                                 backgroundColor = styleInterval(0.05, c('#e09d99', 'white')))
                 })
                 
                 
                 sel_row <- reactive({
                   req(sign_data())
                   row = input$sign_table_rows_selected
                   if (!is.null(row)) {
                     print(sign_data() %>% slice(row))
                   }
                 })
                 
               })
  
}
