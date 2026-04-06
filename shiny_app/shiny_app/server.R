


# Define server logic to summarize and view selected dataset ----
server <- function(input, output) {
  
  # Your application server logic 
  memory.limit(size=50000)
  options(shiny.maxRequestSize=30*1024^2)
  
  # observeEvent(input$close_session,{
  #   session$close()
  #   q(save="no")
  # })
  
  mod_analysis_server("analysis")
  
  
  upload <- mod_upload_analysis_server("upload_analysis")
  
  threshold <- mod_accepted_tracks_server("accepted_tracks",upload)
  
  
  track <- mod_singleTrackSelection_server("trackSelection",upload,threshold)
  
  mod_singleTimeSeries_server("timeSeries",track,upload)
  
  aggregationData <- mod_aggregation_server("aggregation",upload,threshold)
  groupValues <- mod_colorPicker_server("colorPicker",aggregationData)
  
  # Pooling/Filtering/Colors
  newaggregation <- reactive({
    if(is.null(aggregationData())){
      data <- data.frame()
    }else{
      
      data <- aggregationData()
      
    }
    if(!is.null(unlist(groupValues()$groupNameValues))){
      data <- data %>% 
        filter(group_condition %in% groupValues()$groupPrevValues)
      #TODO: renaming!
      data$group_condition = as.factor(data$group_condition)
      # levels(data$group_condition) <- unlist(groupValues()$groupNameValues)
      data
    }else{
      data
    }
  })
  
  # selecting variables
  selectedVariables <- mod_variableSelection_server("variableSelection")
  
  binnedData <- mod_1dbinning_sidebar_server("binning",upload)
  mod_1dbinning_plots_server("binning",upload,binnedData)
  mod_2dbinning_server("2dbinning",upload)
  
  # Distributions
  mod_aggregationTable_server("aggregationTable",newaggregation)
  mod_boxplots_server("boxplots",newaggregation,groupValues)
  mod_uTests_server("utests",newaggregation,selectedVariables)
  
  #Relationships
  mod_scatterMulti_server("scatterMulti",newaggregation,groupValues)
  mod_correlationAnalysis_server("correlation",newaggregation,selectedVariables)
  mod_pca_server("pca",newaggregation,selectedVariables)
  mod_parcoord_server("parcoord",newaggregation,selectedVariables)
  
  #Machine Learning
  mod_randomForest_server("randomForest",newaggregation,selectedVariables)
  
}