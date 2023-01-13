#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @import shinydashboard
#' @import shinyWidgets
#' @import shinyFiles
#' 
#' @noRd
app_ui <- function(request) {
  tagList(
    # Leave this function for adding external resources
    golem_add_external_resources(),
    # Your application UI logic
    
    
    fluidPage(
      #dashboardHeader(title = "ILLAV"),
      
      # dashboardSidebar(
      #   sidebarMenu(
      #     id = "sidebarmenu",
      #     menuItem(
      #       "Analyze data",
      #       tabName = "analysis",
      #       icon = icon("database")
      #     ),
      #     menuItem(
      #       "Upload analysis",
      #       tabName = "upload_analysis",
      #       icon = icon("database")
      #     ),
      #     menuItem(
      #       "Aggregation",
      #       tabName = "aggregation",
      #       icon = icon("database")
      #     ),
      #     menuItem(
      #       "Single mode",
      #       tabName = "single",
      #       badgeColor = "green",
      #       icon = icon("th")
      #       
      #     ),
      #     # menuItem(
      #     #   "Aggregation Visualization",
      #     #   tabName = "aggregation_vis",
      #     #   badgeColor = "green",
      #     #   icon = icon("chart-bar")
      #     # ),
      #     menuItem(
      #       "Binning",
      #       tabName = "binning",
      #       badgeColor = "green",
      #       menuSubItem("1D-Binning", tabName = "1dbin"),
      #       menuSubItem("2D-Binning", tabName = "2dbin"),
      #       icon = icon("chart-line")
      #     )
      #     
      #   )
      # ),
      
      
        shinyjs::useShinyjs(),
        tags$style(
          type = "text/css",
          ".shiny-output-error { visibility: hidden; }",
          ".shiny-output-error:before { visibility: hidden; }"
        ),
        tags$head(tags$style(
          HTML("#dashboard{margin-bottom:2000px;}")
        )),
        
        
        tabsetPanel(
          tabPanel("Analysis Processing",
                  mod_analysis_ui("analysis")
                  ),
          tabPanel("Upload Analysis",
                  wellPanel(tabsetPanel(
                    tabPanel("Upload",
                             mod_upload_analysis_ui("upload_analysis")),
                    tabPanel(
                      "Accepted tracks",
                      mod_accepted_tracks_ui("accepted_tracks")
                    )
                  ))),
          tabPanel("Aggregation",
                  wellPanel(
                      h3("Aggregation"),
                      

                               tags$div(
                                 "Using the aggregation one can create a data set that contains behavioral variables that are calculated
                                      either for each Petri dish or each individual animal of the experiment."
                               ),
                               tags$div(
                                 "Using the grouping option one can select if it should be grouped by Petri dish or individual animal.
                                      With the filters one can filter out frames based on different time series values.
                                      When clicking on update aggregation the data set will be created. The data set can be downloaded with the Download button.
                                      Also it is possible to use an upload file by selecting the Use upload option.
                                      "
                               ),
                               
                               tabsetPanel(
                                 tabPanel("Create Aggregation",
                                          fluidRow(
                                            column(6,
                                                   
                                                   mod_aggregation_ui("aggregation")
                                                   ),
                                            column(6,
                                                   
                                                   mod_aggregationTable_ui("aggregationTable")
                                                   )
                                            
                                          )),
                                 
                                 tabPanel("Group Selection",
                                          mod_colorPicker_ui("colorPicker")),
                                 tabPanel(
                                   "Variables",
                                   mod_variableSelection_ui("variableSelection")
                                 )
                               )
                  ),
                  wellPanel(
                    
                    tabsetPanel(
                      tabPanel(
                        "Distributions",
                        fluidRow(
                          column(6,
                                 
                                 mod_boxplots_ui("boxplots")
                          ),
                          column(6,
                                 mod_uTests_ui("utests")
                          )
                        )
                      ),
                      tabPanel("Correlation",
                               fluidRow(
                                 column(6,
                                        
                                        mod_scatterMulti_ui("scatterMulti")
                                 ),
                                 column(6,
                                        mod_correlationAnalysis_ui("correlation")
                                 )
                                 
                               )
                               ),
                      tabPanel(
                        "Parallel Coordinates",
                        
                        mod_parcoord_ui("parcoord")
                      ),
                      tabPanel(
                        "Random Forest",
                        
                        mod_randomForest_ui("randomForest")
                      )
                      
                    )
                    


                  )
                  
                  
                  ),
          
          tabPanel("Individual Mode",
                  
                   
                    fluidRow(
                      column(4,
                             
                             mod_singleTrackSelection_ui("trackSelection")
                             ),
                      column(8,
                             mod_singleTimeSeries_ui("timeSeries")
                             )
                    )
               
                    
                  ),
          
          tabPanel("Binning",
                  fluidRow(
                    column(4,
                           wellPanel(mod_1dbinning_sidebar_ui("binning"))),
                    column(8,
                           
                           mod_1dbinning_plots_ui("binning"))
                    
                  ))

        
      
    )
  ))
}

#' Add external Resources to the Application
#' 
#' This function is internally used to add external 
#' resources inside the Shiny application. 
#' 
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function(){
  
  add_resource_path(
    'www', app_sys('app/www')
  )
 
  tags$head(
    favicon(),
    bundle_resources(
      path = app_sys('app/www'),
      app_title = 'ILLAV1'
    )
    # Add here other external resources
    # for example, you can add shinyalert::useShinyalert() 
  )
}

