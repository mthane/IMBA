
mod_variableSelection_ui <- function(id) {
  ns <- NS(id)
  
  fluidRow(
    column(2,
           actionButton(ns("select_ts"),"Select all"),
           actionButton(ns("select_none_ts"),"Select none"),
           checkboxGroupInput(ns("TS_Variables"),
                              label = "TS-Variables", 
                              choices = global_vars()$TS_VARIABLES,
                              selected = global_vars()$TS_VARIABLES
           )
    ),
    column(2,
           actionButton(ns("select_run"),"Select all"),
           actionButton(ns("select_none_run"),"Select none"),
           checkboxGroupInput(ns("RUN_Variables"),
                              label = "RUN-Variables", 
                              choices = global_vars()$RUN_VARIABLES,
                              selected = global_vars()$RUN_VARIABLES
           )
    ),
    column(2,
           actionButton(ns("select_hc"),"Select all"),
           actionButton(ns("select_none_hc"),"Select none"),
           checkboxGroupInput(ns("HC_Variables"),
                              label = "HC-Variables", 
                              choices = global_vars()$HC_VARIABLES,
                              selected = global_vars()$HC_VARIABLES
           )
    ),
    column(2,
           actionButton(ns("select_pref"),"Select all"),
           actionButton(ns("select_none_pref"),"Select none"),
           checkboxGroupInput(ns("PREF_Variables"),
                              label = "PREF-Variables", 
                              choices = global_vars()$PREF_VARIABLES,
                              selected = global_vars()$PREF_VARIABLES
           )
    ),
    column(2,
           actionButton(ns("select_track"),"Select all"),
           actionButton(ns("select_none_track"),"Select none"),
           checkboxGroupInput(ns("TRACK_Variables"),
                              label = "TRACK-Variables", 
                              choices = global_vars()$TRACK_VARIABLES,
                              selected = global_vars()$TRACK_VARIABLES
           )
    )
    
    
  )
}


mod_variableSelection_server <- function(id){
  moduleServer(id,
               function(input, output, session) {
                 
                 
                 ## TS 
                 observeEvent(input$select_ts,{
                   
                   updateCheckboxGroupInput(
                     session, 'TS_Variables', choices = global_vars()$TS_VARIABLES,
                     selected = TS_VARIABLES
                   )
                   
                 })
                 
                 observeEvent(input$select_none_ts,{
                   
                   updateCheckboxGroupInput(
                     session, 'TS_Variables', choices = global_vars()$TS_VARIABLES,
                     selected = c()
                   )
                   
                 })
                 
                 
                 ## RUN
                 observeEvent(input$select_run,{
                   updateCheckboxGroupInput(
                     session, 'RUN_Variables', choices = global_vars()$RUN_VARIABLES,
                     selected = RUN_VARIABLES
                   )
                   
                 })
                 
                 observeEvent(input$select_none_run,{
                   updateCheckboxGroupInput(
                     session, 'RUN_Variables', choices = global_vars()$RUN_VARIABLES,
                     selected = c()
                   )
                 })
                 
                 observeEvent(input$select_hc,{
                   updateCheckboxGroupInput(
                     session, 'HC_Variables', choices = global_vars()$HC_VARIABLES,
                     selected = HC_VARIABLES
                   )
                 })
                 
                 observeEvent(input$select_none_hc,{
                   updateCheckboxGroupInput(
                     session, 'HC_Variables', choices = global_vars()$HC_VARIABLES,
                     selected = c()
                   )
                   
                 })
                 
                 observeEvent(input$select_pref,{
                   updateCheckboxGroupInput(
                     session, 'PREF_Variables', choices = global_vars()$PREF_VARIABLES,
                     selected = PREF_VARIABLES
                   )
                   
                 })
                 
                 observeEvent(input$select_none_pref,{
                   updateCheckboxGroupInput(
                     session, 'PREF_Variables', choices = global_vars()$PREF_VARIABLES,
                     selected = c()
                   )
                 })
                 
                 observeEvent(input$select_track,{
                   updateCheckboxGroupInput(
                     session, 'TRACK_Variables', choices = global_vars()$TRACK_VARIABLES,
                     selected = TRACK_VARIABLES
                   )
                   
                 })
                 
                 observeEvent(input$select_none_track,{
                   updateCheckboxGroupInput(
                     session, 'TRACK_Variables', choices = global_vars()$TRACK_VARIABLES,
                     selected = c()
                   )
                 })
                 
                 
                 variables <- reactive(c(input$TS_Variables,
                                         input$RUN_Variables,
                                         input$HC_Variables,
                                         input$PREF_Variables,
                                         input$TRACK_Variables
                 )
                 
                 )
                 variables
               }
               
               
  )
}
